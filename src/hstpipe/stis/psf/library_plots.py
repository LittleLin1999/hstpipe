"""Plot existing FITS PSF templates without reducing or stacking new data."""
from pathlib import Path
import hashlib,json
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .extraction import fractional_aperture


def common_window_density(y,profile,covariance,radius,scale=.0246):
    """Finite-aperture normalization, quotient covariance, density per arcsec."""
    y=np.asarray(y,float);p=np.asarray(profile,float);cov=np.asarray(covariance,float)
    if len(y)<2 or np.any(np.diff(y)<=0):raise ValueError('Offsets must increase')
    widths=np.diff(y)
    if not np.allclose(widths,widths[0],rtol=1e-6):raise ValueError('Nonuniform spatial sampling')
    if y[0]-widths[0]/2>-radius+1e-6 or y[-1]+widths[0]/2<radius-1e-6:raise ValueError('Requested normalization window is not covered')
    aperture=fractional_aperture(y,2*radius);den=float(aperture@p)
    if den<=0 or not np.isfinite(den):raise ValueError('Invalid common-window flux')
    jac=np.eye(len(p))/den-p[:,None]*aperture[None,:]/den**2
    normalized=p/den;normalized_cov=jac@cov@jac.T
    width_arcsec=widths[0]*scale
    return y*scale,normalized/width_arcsec,normalized_cov/width_arcsec**2


def select_templates(hdul,band,sampling,apertures,radius):
    """One existing WING group per slit: favor independent stars, then repeats."""
    result=[]
    for ap in apertures:
        candidates=[h for h in hdul if h.name.startswith('PSF') and h.header.get('APERTURE')==ap and h.header.get('SAMPLING')==sampling and h.header.get('GRADE')=='WING' and h.header.get('NORMRAD',0)>=radius and h.header.get('WAVEMIN')==band[0] and h.header.get('WAVEMAX')==band[1]]
        if candidates:
            chosen=max(candidates,key=lambda h:(h.header['NSTARS'],h.header['NEXP'],h.header['EPOCH'],h.name))
            result.append(chosen)
    return result


def library_diagnostics(c):
    library=Path(c.get('plot_library','psf_library/stis_g140l_psf_library.fits'));out=Path('figures');out.mkdir(exist_ok=True)
    radius=c.get('plot_normalization_radius',50);comparison_radius=c.get('plot_comparison_radius',radius-.5);apertures=c.get('plot_apertures',['52X0.2','52X0.5','52X2'])
    bands=[tuple(b) for b in c['bands'] if b[1]-b[0]<100];focus=tuple(c.get('plot_focus_band',[1520,1580]));scale=.0246
    colors=dict(zip(apertures,['#009E73','#D55E00','#0072B2']));styles=dict(zip(apertures,[':','--','-']))
    modes=[('hstpipe_native','Native detector'),('drizzled','hstpipe drizzle')];records=[];points=[]
    with fits.open(library,checksum=True) as hd:
        idx={r['extension']:r for r in hd['INDEX'].data}
        targets={r['root']:r['target'] for r in hd['PROVENANCE'].data}
        def draw(ax,band,mode,detail=False,core=False,collect=False):
            selected=select_templates(hd,band,mode,apertures,radius)
            for h in selected:
                n=h.name[3:];x,density,cov=common_window_density(h.data['offset_pixel'],h.data['profile'],hd['PCOV'+n].data,comparison_radius)
                err=np.sqrt(np.maximum(np.diag(cov),0));ap=h.header['APERTURE'];label=ap
                target_names=sorted({targets.get(root,root) for root in str(idx[h.name]['roots']).split(',')})
                if detail:label+=f" | {h.header['EPOCH']} | "+(', '.join(target_names) if len(target_names)==1 else f"{len(target_names)} stars")+f", {h.header['NEXP']} exp."
                ax.plot(x,density,color=colors[ap],ls=styles[ap],lw=1.65,label=label)
                ax.fill_between(x,density-err,density+err,color=colors[ap],alpha=.13,lw=0)
                if collect:
                    row=dict(extension=h.name,aperture=ap,band=f'{band[0]}-{band[1]}',sampling=mode,epoch=h.header['EPOCH'],n_stars=h.header['NSTARS'],n_exposures=h.header['NEXP'],trace_y_min=h.header['YMIN'],trace_y_max=h.header['YMAX'],normalization_radius=h.header['NORMRAD'],comparison_radius=comparison_radius,grade=h.header['GRADE'],roots=idx[h.name]['roots'],targets=','.join(target_names))
                    records.append(row)
                    points.extend(dict(extension=h.name,aperture=ap,band=row['band'],sampling=mode,offset_arcsec=float(xx),density_arcsec_inverse=float(pp),prediction_sigma_arcsec_inverse=float(ee)) for xx,pp,ee in zip(x,density,err))
            ax.axhline(0,color='.5',lw=.55)
            ax.set_xlim((-.3,.3) if core else (-comparison_radius*scale,comparison_radius*scale))
            if core:ax.set_ylim(bottom=-.1)
            else:
                ax.set_yscale('symlog',linthresh=.002);ax.set_ylim(-.03,20)
            ax.set_xlabel('Along-slit offset [arcsec]');ax.set_ylabel('Normalized profile [arcsec⁻¹]')
            ax.grid(alpha=.13);ax.legend(fontsize=7,loc='upper right',framealpha=.9)
            missing=[ap for ap in apertures if ap not in [h.header['APERTURE'] for h in selected]]
            if missing:ax.text(.02,.03,'No matching wing-grade template: '+', '.join(missing),transform=ax.transAxes,fontsize=7,color='.35')
        with plt.rc_context({'font.size':10,'axes.titlesize':11,'savefig.facecolor':'white'}):
            fig,axes=plt.subplots(2,2,figsize=(13,8.2));fig.subplots_adjust(left=.075,right=.985,bottom=.13,top=.91,hspace=.39,wspace=.22)
            for j,(mode,label) in enumerate(modes):
                draw(axes[j,0],focus,mode,detail=True,core=True)
                draw(axes[j,1],focus,mode,core=False)
                axes[j,0].set_title(label+': core');axes[j,1].set_title(label+': core and wings (signed log scale)')
            fig.suptitle(f'STIS G140L spatial PSF by slit width | {focus[0]}–{focus[1]} Å',fontsize=15)
            fig.text(.075,.035,f'Flux normalized within ±{comparison_radius*scale:.3f} arcsec (±{comparison_radius:g} native pixels); shading: pointwise ±1σ predictive uncertainty.\nStars/epochs/positions differ; rows select independently. No slit correction; single-star scatter is unmeasured.',fontsize=9)
            for suffix in ['png']:fig.savefig(out/('psf_slit_width_offset.'+suffix),dpi=180)
            plt.close(fig)
            for mode,label in modes:
                ncols=3;nrows=(len(bands)+ncols-1)//ncols
                fig,axes=plt.subplots(nrows,ncols,figsize=(15,4*nrows+1),squeeze=False);fig.subplots_adjust(left=.065,right=.985,bottom=.13,top=.9,hspace=.39,wspace=.25)
                for ax,band in zip(axes.flat,bands):
                    draw(ax,band,mode,collect=True);ax.set_title(f'{band[0]}–{band[1]} Å')
                for ax in list(axes.flat)[len(bands):]:ax.set_visible(False)
                fig.suptitle(label+': wavelength-dependent PSF profiles by slit width',fontsize=15)
                fig.text(.065,.035,f'Common normalization: flux within ±{comparison_radius*scale:.3f} arcsec; signed log intensity; shading: ±1σ predictive uncertainty.\nExisting WING groups only; no new stack. Epoch/position differ; single-star scatter is unmeasured. See psf_slit_plot_selection.ecsv.',fontsize=9)
                for suffix in ['png']:fig.savefig(out/(f'psf_slit_width_wavelength_{mode}.'+suffix),dpi=180)
                plt.close(fig)
    tables=Path('products/measurements');tables.mkdir(parents=True,exist_ok=True)
    Table(rows=records).write(tables/'psf_slit_plot_selection.ecsv',overwrite=True)
    Table(rows=points).write(tables/'psf_slit_plot_profiles.ecsv',overwrite=True)
    meta=dict(library=str(library),library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),minimum_template_radius_native_pixels=radius,comparison_radius_native_pixels=comparison_radius,comparison_radius_arcsec=comparison_radius*scale,selection='Existing WING groups with at least the requested support; renormalized to an exact common radius; maximum independent stars, then exposures, then epoch. No newly combined slit templates.',uncertainty='Pointwise predictive covariance transformed through finite-aperture normalization; between-star scatter unmeasured for single-star groups.',confounders='Epoch, detector position, source spectrum, centering; slit-width transfer remains unvalidated.',n_selected_templates=len(records))
    (tables/'psf_slit_plot_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta
