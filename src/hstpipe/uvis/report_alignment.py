"""Per-exposure alignment previews; original pixels and a fixed HST reference."""
import base64
import html
import io
import json
from pathlib import Path
import numpy as np
from ..common.paths import ReductionPaths
from ..common.provenance import identity
from ..common.target_runner import save


def _tangent(ra,dec):
    from astropy.wcs import WCS
    w=WCS(naxis=2);w.wcs.crpix=[1.,1.];w.wcs.crval=[ra,dec]
    w.wcs.cdelt=[1/3600,1/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    return w


def _panel(path,row,world,frame,halfwidth):
    """Project native pixel cell boundaries, not interpolated detector data."""
    from astropy.io import fits
    from .astrometry.wcs import detector_wcs
    with fits.open(path) as hd:
        version=int(row['extver']);h=hd['SCI',version];w=detector_wcs(hd,version)
        predicted=w.all_world2pix([world],0)[0]
        measured=frame.wcs_world2pix(w.all_pix2world([[float(row['x']),float(row['y'])]],0),0)[0]
        # Local rotation may require a larger native box than the displayed sky box.
        scale=float(np.linalg.svd(w.pixel_scale_matrix,compute_uv=False).min()*3600)
        radius=int(np.ceil(halfwidth/scale*1.5))+2
        x,y=np.rint(predicted).astype(int);ny,nx=h.shape
        x0,x1=max(0,x-radius),min(nx,x+radius+1);y0,y1=max(0,y-radius),min(ny,y+radius+1)
        if x0>=x1 or y0>=y1:raise ValueError('Reference is outside this SCI coverage')
        cut=np.asarray(h.data[y0:y1,x0:x1],float).copy()
        rate=str(h.header.get('BUNIT','')).upper() in ('ELECTRONS/S','ELECTRON/S','COUNTS/S')
        if not rate:cut/=float(hd[0].header['EXPTIME'])
        yy,xx=np.meshgrid(np.arange(y0,y1+1)-.5,np.arange(x0,x1+1)-.5,indexing='ij')
        cells=frame.wcs_world2pix(w.all_pix2world(np.column_stack((xx.ravel(),yy.ravel())),0),0)
        return dict(data=cut,x=cells[:,0].reshape(xx.shape),y=cells[:,1].reshape(yy.shape),
                    offset_mas=(measured*1000).tolist(),extver=version,ccdchip=h.header.get('CCDCHIP'),wcsname=h.header.get('WCSNAME',''))


def exposure_cutouts(paths,record):
    """Read saved FLCs/catalogs only; validate their checkpoint before plotting."""
    from astropy.io import fits
    from astropy.table import Table
    from .astrometry.wcs import detector_wcs
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}))
    entries=rel.get('detector_catalogs',[])
    raw={Path(i['path']).name:Path(i['path']) for i in record.get('inputs',[])}
    expected={Path(i['path']).name:i for i in a.get('tweakback',{}).get('inputs',a.get('stages',{}).get('provisional_drizzle',{}).get('inputs',[])) if isinstance(i,dict)}
    inputs=[];cats=[];dependencies=[]
    for item in entries:
        name=item['exposure'];before=raw.get(name);after=paths.area('intermediate','working')/name;catpath=Path(item['path'])
        if before is None or not all(p.is_file() for p in (before,after,catpath)):continue
        # A later absolute or unrelated working WCS must not masquerade as the
        # accepted relative checkpoint. Missing historical files produce no panel.
        if name in expected and identity(after)!=expected[name]:raise ValueError('Relative working checkpoint has changed: '+name)
        cat=Table.read(catpath,format='ascii.ecsv')
        if 'common_source_id' not in cat.colnames:continue
        inputs.append((before,after));cats.append(cat);dependencies.extend((before,after,catpath))
    if len(inputs)<2:return dict(status='unavailable',figures=[])
    cache=paths.cache();receipt=paths.file('intermediate','alignment','exposure_alignment_preview')
    key=cache.key(dependencies+[Path(__file__)],dict(relative_status=rel.get('status'),maximum_sources=3))
    if cache.valid('exposure_alignment_preview',key):return json.loads(receipt.read_text())
    by_source={}
    for i,cat in enumerate(cats):
        for row in cat:
            if 'centroid_valid' in cat.colnames and not bool(row['centroid_valid']):continue
            if not np.isfinite([row['x'],row['y'],row['centroid_error_pix']]).all():continue
            by_source.setdefault(str(row['common_source_id']),{})[i]=row
    chosen=sorted((sid for sid,rows in by_source.items() if len(rows)>=2),
                  key=lambda sid:(not sid.startswith('G'),-len(by_source[sid]),
                                  np.median([float(r['centroid_error_pix']) for r in by_source[sid].values()]),sid))[:3]
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from astropy.visualization import ImageNormalize,AsinhStretch
    figures=[]
    for index,sid in enumerate(chosen):
        rows=by_source[sid];anchor=min(rows,key=lambda i:float(rows[i]['centroid_error_pix']))
        row=rows[anchor]
        with fits.open(inputs[anchor][1]) as hd:
            world=detector_wcs(hd,int(row['extver'])).all_pix2world([[float(row['x']),float(row['y'])]],0)[0]
        frame=_tangent(*world);offsets=[]
        for i,row in rows.items():
            for path in inputs[i]:
                with fits.open(path) as hd:
                    sky=detector_wcs(hd,int(row['extver'])).all_pix2world([[float(row['x']),float(row['y'])]],0)
                    offsets.append(frame.wcs_world2pix(sky,0)[0])
        halfwidth=max(.32,float(np.max(np.abs(offsets)))+.24)
        panels={};measurements=[]
        for i,row in rows.items():
            for phase,path in enumerate(inputs[i]):panels[i,phase]=_panel(path,row,world,frame,halfwidth)
        samples=np.concatenate([p['data'][np.isfinite(p['data'])] for p in panels.values()])
        lo,hi=np.percentile(samples,[10,99.5]);norm=ImageNormalize(vmin=float(lo),vmax=float(max(hi,lo+1e-9)),stretch=AsinhStretch(.1))
        fig,axes=plt.subplots(len(inputs),2,figsize=(6.5,2.5*len(inputs)),squeeze=False,layout='constrained')
        for i,(before,after) in enumerate(inputs):
            item=dict(exposure=before.name,source_id=sid,status='measured' if i in rows else 'not_measured')
            for phase in (0,1):
                ax=axes[i,phase];panel=panels.get((i,phase))
                if panel:
                    ax.pcolormesh(panel['x'],panel['y'],panel['data'],cmap='gray',norm=norm,shading='flat',rasterized=True)
                    ax.plot(0,0,'+',color='#00d9ff',markersize=11,markeredgewidth=1.2)
                    distance=float(np.linalg.norm(panel['offset_mas']))
                    ax.text(.97,.96,f'offset {distance:.1f} mas',transform=ax.transAxes,ha='right',va='top',color='white',fontsize=8,bbox=dict(facecolor='black',alpha=.65,pad=2))
                    item['before' if phase==0 else 'after']={k:panel[k] for k in ('offset_mas','extver','ccdchip','wcsname')}
                else:ax.text(.5,.5,'No reliable source measurement',transform=ax.transAxes,ha='center',wrap=True,fontsize=8)
                ax.set(xlim=(-halfwidth,halfwidth),ylim=(-halfwidth,halfwidth),aspect='equal')
                ax.tick_params(labelsize=7)
                if i==0:ax.set_title('Before: original FLC WCS' if phase==0 else 'After: relative WCS',fontsize=10)
                if phase==0:ax.set_ylabel(before.stem.removesuffix('_flc')+'\nNorth offset (arcsec)',fontsize=8)
                if i==len(inputs)-1:ax.set_xlabel('East offset (arcsec)',fontsize=8)
            measurements.append(item)
        output=paths.file('diagnostics','alignment',f'exposure_reference_source{index+1}','jpg');output.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(output,dpi=110,pil_kwargs=dict(quality=85,optimize=True));plt.close(fig)
        figures.append(dict(source_id=sid,path=str(output),reference_exposure=inputs[anchor][0].name,reference_ra=float(world[0]),reference_dec=float(world[1]),halfwidth_arcsec=halfwidth,measurements=measurements))
    result=dict(status='complete',figures=figures,relative_status=rel.get('status'),
                reference='Fixed sky position from a measured HST source in the accepted relative frame; not an external absolute reference',
                rendering='Original detector pixels drawn through each complete HSTWCS on common east/north axes; no coaddition or new image resampling')
    save(receipt,result);cache.record('exposure_alignment_preview',key,[receipt]+[Path(f['path']) for f in figures]);return result


def previews(paths,filters,budget):
    from PIL import Image
    sections=[];used=0;count=0;esc=lambda x:html.escape(str(x),quote=True)
    for key,record in sorted(filters.items()):
        group=ReductionPaths(paths.root,paths.target,'uvis',record.get('filter',key),record.get('setup'))
        try:comparison=exposure_cutouts(group,record)
        except (OSError,ValueError,KeyError) as exc:
            sections.append('<p>'+esc(key)+': exposure cutouts unavailable — '+esc(exc)+'</p>');continue
        figures=[]
        for item in comparison['figures']:
            with Image.open(item['path']) as original:
                # Reserve equal space for up to three sources in every filter.
                limit=max(2500,(budget-4096*len(filters))//max(1,3*len(filters)))
                for width,quality in ((760,70),(650,60),(560,50),(480,42),(400,35)):
                    image=original.convert('RGB');image.thumbnail((width,int(width*1.8)))
                    stream=io.BytesIO();image.save(stream,format='JPEG',quality=quality,optimize=True)
                    if 4*((len(stream.getvalue())+2)//3)+1024<=limit:break
            encoded=base64.b64encode(stream.getvalue()).decode('ascii');sid=item['source_id']
            caption=f"{key} · {sid}. Fixed reference measured in {item['reference_exposure']}. Each row is one exposure; cyan crosses mark the same reference position."
            fragment='<figure class="alignment-pair"><img alt="'+esc(key+' '+sid+' individual exposures before and after alignment with reference crosses')+'" src="data:image/jpeg;base64,'+encoded+'"><figcaption>'+esc(caption)+'</figcaption></figure>'
            size=len(fragment.encode('utf-8'))
            if used+size+2048>budget:continue
            figures.append(fragment);used+=size;count+=1
        if not figures:continue
        accepted=comparison['relative_status']=='validated'
        state='Accepted relative WCS' if accepted else 'No validated relative correction'
        note='<p><strong>One row per exposure. Left: original WCS. Right: relative WCS.</strong> The cyan cross is the same fixed HST reference position in every panel. Compare each source center with the cross. All panels of a source share angular scale, orientation and brightness stretch; these are individual detector images, not coadds. Missing measurements are labelled.</p>'
        note+='<p class="muted">Native pixels are displayed through the complete distortion-aware WCS. The fixed HST reference checks relative registration, not absolute Gaia accuracy. Cosmic rays remain visible in these input images.</p>'
        a=record.get('scientific_assessment',{})
        if a.get('astrometry',{}).get('final_validation',{}).get('individual_absolute_status')=='review_required':note+='<p class="attention">An individual-exposure external-catalog check still needs review; it is separate from the validated relative alignment.</p>'
        section='<details class="alignment-comparison"><summary>'+esc(key)+' — '+str(len(figures))+' sources, individual exposures · '+state+'</summary>'+note+''.join(figures)+'</details>'
        used+=len(section.encode('utf-8'))-sum(len(f.encode('utf-8')) for f in figures);sections.append(section)
    return sections,count
