"""Read-only STIS acquisition diagnostics with nominal, target-anchored slit overlays."""
from pathlib import Path
import contextlib,io,json,os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import ImageNormalize,AsinhStretch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from .io import digest


def slit_sky_vertices(center,pa_deg,length_arcsec,width_arcsec):
    theta=np.deg2rad(pa_deg);along=np.array([-1,-1,1,1,-1])*length_arcsec/2;cross=np.array([-1,1,1,-1,-1])*width_arcsec/2
    return center.spherical_offsets_by((along*np.sin(theta)+cross*np.cos(theta))*u.arcsec,(along*np.cos(theta)-cross*np.sin(theta))*u.arcsec)


def target_anchored_wcs(header,target_pixel,target_coord):
    w=WCS(header,naxis=2);w.wcs.crpix=np.asarray(target_pixel)+1;w.wcs.crval=[target_coord.ra.deg,target_coord.dec.deg]
    return w


def target_offset_geometry(header,target_pixel,shape):
    """Pixel-edge extent in detector-axis arcseconds about the acquired target."""
    scale=np.linalg.norm(WCS(header,naxis=2).pixel_scale_matrix,axis=0)*3600
    ny,nx=shape;target=np.asarray(target_pixel,float)
    extent=[(-.5-target[0])*scale[0],(nx-.5-target[0])*scale[0],(-.5-target[1])*scale[1],(ny-.5-target[1])*scale[1]]
    return scale,extent


def plot_acquisitions(raw_paths,science_paths,output,cutout_size_arcsec=5.):
    from stistools import tastis
    if not np.isfinite(cutout_size_arcsec) or cutout_size_arcsec<=0:raise ValueError('Positive angular cutout size required')
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    science=[]
    for path in science_paths:
        with fits.open(path) as h:
            p=h[0].header;s=h['SCI',1].header
            science.append(dict(root=p['ROOTNAME'].strip(),acqname=p['ACQNAME'].lower(),aperture=p['APERTURE'],propaper=p.get('PROPAPER',p['APERTURE']),ra_aper=s['RA_APER'],dec_aper=s['DEC_APER'],pa_aper=s['PA_APER'],postarg1=p.get('POSTARG1',0),postarg2=p.get('POSTARG2',0),date=s['DATE-OBS']))
    fig,axes=plt.subplots(len(raw_paths),2,figsize=(12,6*len(raw_paths)),squeeze=False);fig.subplots_adjust(left=.075,right=.93,bottom=.12,top=.91,hspace=.3,wspace=.2)
    reports=[]
    for i,raw in enumerate(raw_paths):
        raw=Path(raw).resolve();spt=raw.with_name(raw.name.replace('_raw.fits','_spt.fits'));before=digest(raw)
        keys=tastis._read_keywords(str(raw),str(spt),spt.exists());tastis._calculate_slews(keys)
        buf=io.StringIO();old=Path.cwd()
        try:
            os.chdir(raw.parent)
            with contextlib.redirect_stdout(buf):tastis.tastis(raw.name,update=False)
        finally:os.chdir(old)
        (out.parent/(raw.stem+'_tastis.txt')).write_text(buf.getvalue())
        assert digest(raw)==before
        matched=[s for s in science if s['acqname'][:8]==raw.name[:8].lower()]
        if not matched:raise ValueError('No science exposures match ACQNAME for '+str(raw))
        slit=matched[0];length,width=map(float,slit['aperture'].upper().split('X'));origin=None;panels=[]
        with fits.open(raw) as h:
            primary=h[0].header;origin=SkyCoord(primary['RA_TARG']*u.deg,primary['DEC_TARG']*u.deg)
            for j,ver in enumerate([1,2]):
                d=h['SCI',ver];header=d.header;data=d.data.astype(float)
                target=np.array([keys['coarse1' if ver==1 else 'fine1'],keys['coarse2' if ver==1 else 'fine2']])-1
                w=target_anchored_wcs(header,target,origin);center=SkyCoord(slit['ra_aper']*u.deg,slit['dec_aper']*u.deg)
                vertices=slit_sky_vertices(center,slit['pa_aper'],length,width);xp,yp=w.world_to_pixel(vertices)
                interior=data[3:-3,3:-3];norm=ImageNormalize(vmin=np.percentile(interior,5),vmax=np.percentile(interior,99.95),stretch=AsinhStretch(.06))
                scale,extent=target_offset_geometry(header,target,data.shape)
                display_polygon=(np.column_stack([xp,yp])-target)*scale
                ax=axes[i,j];ax.set_facecolor('#dedede')
                im=ax.imshow(data,origin='lower',cmap='gray_r',norm=norm,interpolation='nearest',extent=extent)
                ax.add_patch(Polygon(display_polygon,closed=True,edgecolor='#009eaa',facecolor='#009eaa',alpha=.18,lw=1.8))
                ax.plot(display_polygon[:,0],display_polygon[:,1],color='#008b9a',lw=1.8,label='Nominal '+slit['propaper']+' slit')
                ax.plot(0,0,'+',color='#d55e00',ms=13,mew=1.8,label='Onboard acquired target')
                # Directions retain the CCD image orientation; axes are angular
                # offsets from the target, not absolute RA/Dec or north-up.
                half=cutout_size_arcsec/2;anchor_offset=np.array([-.74*half,.62*half]);anchor=target+anchor_offset/scale;sky=w.pixel_to_world(*anchor)
                for label,pa in [('N',0),('E',90)]:
                    xx,yy=w.world_to_pixel(sky.directional_offset_by(pa*u.deg,.5*u.arcsec));end=(np.array([xx,yy])-target)*scale
                    ax.annotate('',xy=end,xytext=anchor_offset,arrowprops=dict(arrowstyle='->',color='#333333',lw=1.4));ax.text(*end,label,color='#333333',ha='center',va='bottom',fontsize=10)
                pixscale=np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600
                barx=-.86*half;bary=-.84*half
                ax.plot([barx,barx+.5],[bary,bary],color='#333333',lw=2);ax.text(barx+.25,bary+.1,'0.5 arcsec',ha='center',fontsize=9)
                ax.set(xlim=(-half,half),ylim=(-half,half),xlabel='CCD x offset from target [arcsec]',ylabel='CCD y offset from target [arcsec]',aspect='equal')
                ax.set_xticks(np.arange(np.ceil(-half),half+1e-8,1.));ax.set_yticks(np.arange(np.ceil(-half),half+1e-8,1.))
                phase='Before coarse slew' if ver==1 else 'After coarse slew / before fine slew'
                ax.set_title(f'{raw.stem.replace("_raw","")} | {phase}\nF28X50LP | {header["DATE-OBS"]} | slit PA={slit["pa_aper"]%360:.3f}°',fontsize=11)
                ax.legend(loc='lower right',fontsize=8,framealpha=.85)
                panels.append(dict(sci_extver=ver,target_pixel_zero_based=target.tolist(),pixel_scale_arcsec=pixscale,axis_scales_arcsec=scale.tolist(),image_extent_target_arcsec=extent,display_limits_arcsec=[-cutout_size_arcsec/2,cutout_size_arcsec/2],slit_polygon_target_arcsec=display_polygon.tolist(),slit_polygon_pixels=np.column_stack([xp,yp]).tolist()))
            report=dict(raw=str(raw),raw_sha256=before,spt=str(spt),target=primary['TARGNAME'],science=matched,panels=panels,tastis_report=buf.getvalue(),coarse_slew_pixels=[keys['a1coarse_pix'],keys['a2coarse_pix']],fine_slew_pixels=[keys['a1fine_pix'],keys['a2fine_pix']])
            reports.append(report)
    fig.suptitle(f'J1022+0841: target-centered {cutout_size_arcsec:g}″ × {cutout_size_arcsec:g}″ acquisition cutouts',fontsize=15)
    fig.text(.075,.035,'Slit overlay uses science RA_APER/DEC_APER and PA_APER; acquisition WCS is anchored to the onboard target centroid.\nThis is a nominal target-anchored footprint, not a direct science-slit image or a measurement of achieved centering.\nThe 52-arcsec slit extends beyond these cutouts; gray margins have no acquisition coverage. CCD orientation is retained (see N/E arrows).',fontsize=9)
    for suffix in ['png']:fig.savefig(out.with_suffix('.'+suffix),dpi=180)
    plt.close(fig)
    metadata=dict(cutout_size_arcsec=cutout_size_arcsec,display_coordinates='CCD-axis angular offsets from onboard acquired target; target at (0,0)',method='Nominal target-anchored slit footprint; no absolute astrometric registration or distortion refinement',references=['https://hst-docs.stsci.edu/stisdhb/chapter-5-stis-data-analysis/5-2-evaluating-target-acquisitions-and-guiding','https://hst-docs.stsci.edu/stisdhb/chapter-2-stis-data-structure/2-6-stis-coordinate-systems'],acquisitions=reports)
    out.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n');return metadata
