"""Empirical spatial PSF dependence on measured source detector position."""
from pathlib import Path
import json,hashlib
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from .extraction import fractional_aperture
from .library_plots import common_window_density


def aperture_ratio(profile,covariance,numerator,denominator):
    """Ratio of two linear apertures with full within-profile covariance."""
    p=np.asarray(profile,float);a=np.asarray(numerator,float);b=np.asarray(denominator,float)
    den=float(b@p);num=float(a@p)
    if den<=0 or not np.isfinite(den):raise ValueError('Nonpositive reference-aperture flux')
    grad=a/den-num*b/den**2
    return num/den,float(np.sqrt(max(grad@covariance@grad,0)))


def position_diagnostics(c):
    band=tuple(c.get('plot_position_band',c.get('plot_focus_band',[1520,1580])));band_name=f'{band[0]}-{band[1]}'
    table=Table.read('products/measurements/psf_measurements.fits')
    rows=table[(table['mode']=='hstpipe_native')&(table['quality']=='PASS')&(table['band']==band_name)]
    if not len(rows):raise ValueError('No accepted native profiles for positional diagnostic '+band_name)
    aperture_list=c.get('plot_apertures',['52X0.2','52X0.5','52X2']);scale=.0246
    core_height=c.get('plot_position_core_height',11);wing_inner=c.get('wing_inner_pixels',10);wing_radius=c.get('plot_position_wing_radius',40)
    colors=dict(zip(aperture_list,['#009E73','#D55E00','#0072B2']));markers=dict(zip(aperture_list,['s','^','o']))
    color_norm=Normalize(float(np.min(rows['trace_y'])),float(np.max(rows['trace_y'])));cm=plt.get_cmap('viridis');records=[];profiles=[];inventory=[]
    for row in rows:
        path=Path(str(row['profile_fits']))
        with fits.open(path,checksum=True) as h:
            data=h['PROFILE'].data;y=np.array(data['offset_pixel'],float);p=np.array(data['normalized'],float);cov=np.array(h['NORMCOV'].data,float)
        core=fractional_aperture(y,core_height)
        xx,density,dcov=common_window_density(y,p,cov,core_height/2)
        concentration,concentration_err=aperture_ratio(p,cov,fractional_aperture(y,3),core)
        wing_eligible=bool(row['wing_usable']) and row['normalization_radius']>=wing_radius
        wing=wing_err=asym=asym_err=np.nan
        if wing_eligible:
            outer=fractional_aperture(y,2*wing_radius);inner=fractional_aperture(y,2*wing_inner)
            wing,wing_err=aperture_ratio(p,cov,outer-inner,core)
            asym,asym_err=aperture_ratio(p,cov,np.sign(y)*outer,outer)
        item={k:str(row[k]) for k in ['root','target','aperture','date','profile_fits']}
        item.update(band=band_name,trace_y=float(row['trace_y']),normalization_radius=float(row['normalization_radius']),snr=float(row['snr']),wing_usable=wing_eligible,fwhm_arcsec=float(row['fwhm'])*scale,fwhm_error_arcsec=float(row['fwhm_error'])*scale,ee3_over_ee11=concentration,ee3_over_ee11_error=concentration_err,wing_over_core=wing,wing_over_core_error=wing_err,asymmetry=asym,asymmetry_error=asym_err)
        records.append(item);profiles.append((item,xx,density,np.sqrt(np.maximum(np.diag(dcov),0))))
        inventory.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    Path('figures').mkdir(exist_ok=True)
    with plt.rc_context({'font.size':10,'axes.titlesize':11,'savefig.facecolor':'white'}):
        fig,axes=plt.subplots(2,len(aperture_list),figsize=(15,8.5),squeeze=False)
        fig.subplots_adjust(left=.065,right=.9,top=.88,bottom=.15,hspace=.35,wspace=.27)
        for j,ap in enumerate(aperture_list):
            group=[q for q in profiles if q[0]['aperture']==ap];wing_group=[q for q in group if q[0]['wing_usable']]
            for item,x,density,err in group:
                color=cm(color_norm(item['trace_y']));axes[0,j].plot(x,density,color=color,lw=.9,alpha=.8)
                if len(group)==1:axes[0,j].fill_between(x,density-err,density+err,color=color,alpha=.16)
                if item['wing_usable']:
                    axes[1,j].plot(x,density,color=color,lw=.9,alpha=.8)
                    if len(wing_group)==1:axes[1,j].fill_between(x,density-err,density+err,color=color,alpha=.16)
            axes[0,j].set_title(f'{ap}: core | {len(group)} exposures')
            axes[1,j].set_title(f'{ap}: wings | {len(wing_group)} qualified exposures')
            if group:
                ymin=min(q[0]['trace_y'] for q in group);ymax=max(q[0]['trace_y'] for q in group)
                axes[0,j].text(.03,.95,f'Source rows {ymin:.1f}–{ymax:.1f}',transform=axes[0,j].transAxes,va='top',fontsize=8)
            for i in [0,1]:
                axes[i,j].set_xlabel('Offset from the measured trace [arcsec]');axes[i,j].set_ylabel(f'Profile / F({core_height:g}-pixel height) [arcsec⁻¹]');axes[i,j].axhline(0,color='.5',lw=.5);axes[i,j].grid(alpha=.13)
            axes[0,j].set_xlim(-.25,.25);axes[0,j].set_ylim(bottom=-.1)
            axes[1,j].set_xlim(-wing_radius*scale,wing_radius*scale);axes[1,j].set_yscale('symlog',linthresh=.002);axes[1,j].set_ylim(-.03,25)
        cb=fig.colorbar(ScalarMappable(norm=color_norm,cmap=cm),cax=fig.add_axes([.925,.22,.015,.58]));cb.set_label('Source position along slit: detector row [native pixel]')
        fig.suptitle(f'STIS G140L PSF shape versus target position along the slit | {band[0]}–{band[1]} Å',fontsize=14)
        fig.text(.065,.045,'Each curve is one accepted native exposure; color is its measured source row, not offset within the PSF.\nCommon normalization: flux in the central 11-pixel extraction height. Bottom row uses wing-qualified profiles only.\nSources, epochs and detector positions are confounded; sparse/disjoint sampling does not establish a pure position or slit-width effect.',fontsize=9)
        for suffix in ['png']:fig.savefig('figures/psf_along_slit_position_profiles.'+suffix,dpi=180)
        plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(13,8));fig.subplots_adjust(left=.09,right=.98,top=.9,bottom=.15,hspace=.36,wspace=.23)
        specs=[('fwhm_arcsec','fwhm_error_arcsec','Direct core FWHM [arcsec]',False),('ee3_over_ee11','ee3_over_ee11_error','Flux in 3 pixels / flux in 11 pixels',False),('wing_over_core','wing_over_core_error',f'Flux at {wing_inner:g}<|offset|<{wing_radius:g} pixels / F(11 pixels)',True),('asymmetry','asymmetry_error',f'(Right − left) / flux within ±{wing_radius:g} pixels',True)]
        for ax,(key,error,label,wing_only) in zip(axes.flat,specs):
            for ap in aperture_list:
                group=[q for q in records if q['aperture']==ap and (q['wing_usable'] or not wing_only) and np.isfinite(q[key])]
                ax.errorbar([q['trace_y'] for q in group],[q[key] for q in group],yerr=[q[error] for q in group],fmt=markers[ap],ms=4,color=colors[ap],lw=.7,capsize=1.5,label=f'{ap} (n={len(group)})')
            ax.set_xlabel('Target position along slit [detector row, native pixel]');ax.set_ylabel(label);ax.grid(alpha=.15);ax.legend(fontsize=8)
            if key=='asymmetry':ax.axhline(0,color='.5',lw=.6)
        fig.suptitle(f'PSF measurements versus target position | {band[0]}–{band[1]} Å | native exposures',fontsize=14)
        fig.text(.09,.04,'Error bars propagate each exposure’s covariance and background bound; they do not measure all between-star systematics.\nNo fitted trend or slit correction: stars/epochs differ and position coverage is sparse. Detector row is not a calibrated sky offset from the slit center.',fontsize=9)
        for suffix in ['png']:fig.savefig('figures/psf_along_slit_position_metrics.'+suffix,dpi=180)
        plt.close(fig)
    directory=Path('products/measurements');Table(rows=records).write(directory/'psf_along_slit_position.ecsv',overwrite=True)
    coverage={ap:dict(n_exposures=len([r for r in records if r['aperture']==ap]),n_wing_exposures=len([r for r in records if r['aperture']==ap and r['wing_usable']]),detector_rows=sorted({r['trace_y'] for r in records if r['aperture']==ap}),stars=sorted({r['target'] for r in records if r['aperture']==ap})) for ap in aperture_list}
    metadata=dict(band=band_name,source_position_coordinate='Median measured detector trace row in native pixels; not sky offset from aperture center',profile_normalization_height_pixels=core_height,wing_inner_pixels=wing_inner,wing_outer_pixels=wing_radius,coverage=coverage,input_profile_sha256=inventory,limitations='No controlled same-star scan or universal position/slit-width correction; conditional per-exposure errors. Regular apertures kept distinct.')
    (directory/'psf_along_slit_position_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata
