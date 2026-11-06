"""Parameter sensitivity using unchanged hstpipe coordinate/drizzle algorithms."""
from pathlib import Path
import contextlib,copy,logging,json
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .common import *
from .extraction import measures
import matplotlib.pyplot as plt


def rectified_metrics(path,lo,hi,radius,c,pixfrac):
    with fits.open(path) as h:
        hh=h['SCI'].header;wave=hh['CRVAL1']+(np.arange(h['SCI'].data.shape[1])+1-hh['CRPIX1'])*hh['CD1_1'];y=(hh['CRVAL2']+(np.arange(h['SCI'].data.shape[0])+1-hh['CRPIX2'])*hh['CD2_2'])/.0246
        rows=np.where(abs(y)<=radius)[0];cols=np.where((wave>=lo)&(wave<=hi))[0]
        for a,b in c['exclude_wavelengths']:cols=cols[~((wave[cols]>=a)&(wave[cols]<=b))]
        im=h['SCI'].data[np.ix_(rows,cols)];dq=h['DQ'].data[np.ix_(rows,cols)];err=h['ERR'].data[np.ix_(rows,cols)];valid=(dq==0)&np.isfinite(im)&np.isfinite(err);p=np.sum(np.where(valid,im,0),axis=1);cv=np.diag(np.sum(np.where(valid,err**2,0),axis=1)*max(1,pixfrac*.0246/hh['CD2_2']))
        mm,norm,ncov=measures(y[rows],p,cv,[3,5,7,9,11,15,21,41,81,140],wing_inner=c.get('wing_inner_pixels',10))
        half=hh['CD2_2']/.0246/2
        for height in [3,5,7,9,11,15,21,41,81,140]:
            aperture=np.clip((np.minimum(y[rows]+half,height/2)-np.maximum(y[rows]-half,-height/2))/(2*half),0,1)
            mm[f'ee_{height}']=float(aperture@norm) if height<=2*radius else float('nan')
        mm['integrated_flux']=float(p.sum()*hh['CD2_2']*hh['CD1_1'])
        return mm


def run(c):
    from hstpipe.stis.config import Config
    from hstpipe.stis.geometry import Grid,make_grid
    from hstpipe.stis.io import read_inputs
    from hstpipe.stis.workflow import stage3_joint_drizzle
    reds=load_json('data/catalogs/hstpipe_reductions.json',[]);measurements=load_json('products/measurements/hstpipe_psfs.json',[]);selected={};results=[]
    for r in reds:
        if r.get('status')!='COMPLETE':continue
        key=(r['aperture'],int(np.rint(r['trace_y_min']/c.get('detector_position_group_pixels',100))),r['target']=='GRW+70D5824')
        selected.setdefault(key,r)
    for red in selected.values():
        base=Path(red['directory']);conf=Config(**load_json(base/'config.json'));inputs=read_inputs(sorted(Path(conf.project_dir).glob(conf.input_glob)));meta=load_json(base/'intermediate/flt/registration/inputs.json');basegrid=Grid(**meta['grid']);baseline=base/'products'/conf.product_name('s2d_drz')
        variants=[(p,conf.spatial_scale) for p in c.get('drizzle_test_pixfrac',[.6,.8])]+[(conf.pixfrac,s) for s in c.get('drizzle_test_spatial_scale',[.0246])]
        for pixfrac,scale in variants:
            label=f'pixfrac{pixfrac:.2f}_scale{scale:.5f}';directory=base/'diagnostics/drizzle_sensitivity'/label;directory.mkdir(parents=True,exist_ok=True);cfg=copy.deepcopy(conf);cfg.output_dir=str(directory.resolve());cfg.pixfrac=pixfrac;cfg.spatial_scale=scale
            headers=[fits.getheader(cfg.path(cfg.control_dir)/(e.path.stem+'_perform.fits'),'SCI') for e in inputs];offsets=[e.p['POSTARG2']+cfg.source_offsets.get(e.path.name,0) for e in inputs];grid=make_grid(headers,cfg,offsets)
            (directory/'products').mkdir(exist_ok=True)
            output=directory/'products'/cfg.product_name('s2d_drz');signature=dict(baseline_sha256=sha(baseline),pixfrac=pixfrac,spatial_scale=scale)
            if load_json(directory/'cache.json')!=signature or not output.exists():
                background=directory/'intermediate/background';background.parent.mkdir(exist_ok=True)
                if not background.exists():background.symlink_to((base/'intermediate/background').resolve(),target_is_directory=True)
                registration=directory/'intermediate/flt/registration';registration.mkdir(parents=True,exist_ok=True)
                for e in inputs:
                    with np.load(base/'intermediate/flt/registration'/(e.path.stem+'.npz')) as d:
                        data={k:d[k] for k in d.files};data['mapping']=grid.map(data['wavelength'],data['spatial']);np.savez_compressed(registration/(e.path.stem+'.npz'),**data)
                with (directory/'drizzle.log').open('w') as stream,contextlib.redirect_stdout(stream),contextlib.redirect_stderr(stream):stage3_joint_drizzle(cfg,inputs,grid)
                save_json(directory/'cache.json',signature);save_json(directory/'config.json',cfg.serializable())
            for lo,hi in [[1420,1480],[1330,1680]]:
                radii=[r['normalization_radius'] for r in measurements if r['reduction_group']==red['group'] and r['mode']=='drizzled' and r['lo']==lo and r['hi']==hi]
                if not radii:continue
                radius=min(radii)
                try:
                    a=rectified_metrics(baseline,lo,hi,radius,c,conf.pixfrac);b=rectified_metrics(output,lo,hi,radius,c,pixfrac)
                    flux_change=b['integrated_flux']/a['integrated_flux']-1
                    results.append(dict(group=red['group'],target=red['target'],aperture=red['aperture'],band=f'{lo}-{hi}',normalization_radius=radius,pixfrac=pixfrac,spatial_scale_arcsec=scale,fwhm_change=b['fwhm']-a['fwhm'],ee11_change=b['ee_11']-a['ee_11'],wing_fraction_change=b['wing_fraction']-a['wing_fraction'],relative_integrated_flux_change=flux_change,flux_gate='PASS' if abs(flux_change)<.03 else 'WARNING',product=str(output)))
                except Exception:logging.exception('Drizzle sensitivity %s %s',red['group'],label)
    save_json('products/measurements/drizzle_sensitivity.json',results)
    if results:
        Table(rows=results).write('products/measurements/drizzle_sensitivity.ecsv',overwrite=True)
        fig,ax=plt.subplots(1,2,figsize=(12,4))
        for ap in sorted(set(r['aperture'] for r in results)):
            rr=[r for r in results if r['aperture']==ap];ax[0].scatter([r['pixfrac'] for r in rr],[r['fwhm_change']*.0246 for r in rr],s=18,label=ap);ax[1].scatter([r['spatial_scale_arcsec'] for r in rr],[r['wing_fraction_change'] for r in rr],s=18,label=ap)
        ax[0].set(xlabel='Drizzle pixfrac',ylabel='FWHM change [arcsec]');ax[1].set(xlabel='Spatial output scale [arcsec/pixel]',ylabel='Wing-fraction change at common radius');ax[0].legend();ax[1].legend();fig.tight_layout();fig.savefig('figures/drizzle_sensitivity.png',dpi=140);plt.close(fig)
    return results
