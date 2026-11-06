"""Configuration-specific finite-window PSF library and normalized covariance."""
from collections import defaultdict
from astropy.io import fits
from astropy.table import Table
import numpy as np
from .common import *


def eligible(m,rows,wing=False,modes=('hstpipe_native',)):
    byroot={r['root']:r for r in rows}
    return [r for r in m if r['background']==0 and r['mode'] in modes and not r['quality'] and byroot.get(r['root'],{}).get('point_source_prior',False) and (r['wing_usable'] or not wing)]


def combine(items):
    # Equal independent-star weights, with both repeat and inter-star scatter.
    stars=defaultdict(list)
    for r in items:stars[r['target']].append(r)
    profiles=[];mean_covariances=[];predictive_covariances=[];repeat_scatter=[]
    for target,group in stars.items():
        d=[np.load(r['profile_file']) for r in group];w=np.array([max(r['snr'],1)**2 for r in group],dtype=float);w/=w.sum()
        values=np.array([x['normalized'] for x in d]);mean=np.sum(w[:,None]*values,axis=0)
        measurement=np.sum([ww**2*x['normalized_covariance'] for ww,x in zip(w,d)],axis=0)
        within=np.cov(values,rowvar=False,ddof=1) if len(values)>1 else np.zeros_like(measurement)
        # Undeblended sample scatter is conservative: photon noise is not removed.
        neff=1/np.sum(w*w)
        profiles.append(mean);mean_covariances.append(measurement+within/neff)
        predictive_covariances.append(measurement+within);repeat_scatter.append(within)
    profiles=np.array(profiles);mean=np.mean(profiles,axis=0)
    between=np.cov(profiles,rowvar=False,ddof=1) if len(profiles)>1 else np.zeros_like(mean_covariances[0])
    cov_mean=np.sum(mean_covariances,axis=0)/len(profiles)**2+between/len(profiles)
    predictive=np.mean(predictive_covariances,axis=0)+between
    scatter=between+np.mean(repeat_scatter,axis=0)
    return mean,cov_mean,predictive,scatter,len(stars)


def build_library(c):
    rows=load_json('data/catalogs/observations.json');m=load_json('products/measurements/profiles.json');good=eligible(m,rows,modes=('hstpipe_native','drizzled'))
    if not good:raise RuntimeError('No hstpipe PSF profiles pass QA; refusing empty library')
    groups=defaultdict(list)
    for r in good:
        key=(r['aperture'],r['band'],str(r['date'])[:4],int(np.rint(r['trace_y']/c.get('detector_position_group_pixels',100))),r['mode'],r.get('normalization_radius',70))
        groups[key].append(r)
    hdr=fits.Header();hdr['ORIGIN']='Empirical recent STIS G140L';hdr['NORMTYPE']='FINITE_PER_GROUP';hdr['NORMUNIT']='detector pixel';hdr['WINGIN']=c.get('wing_inner_pixels',10);hdr['POSBIN']=c.get('detector_position_group_pixels',100);hdr['WINGNOTE']='Wing grade requires exposure-level wing S/N gate';hdr['EXTRAPOL']='PROHIBITED'
    hdus=[fits.PrimaryHDU(header=hdr)];index=[]
    for k,items in sorted(groups.items()):
        allroots=sorted({root for item in items for root in item.get('roots',[item['root']])})
        mean,cm,pred,scatter,nstars=combine(items);y=np.load(items[0]['profile_file'])['y'];number=len(index)+1
        table=Table(dict(offset_pixel=y,offset_arcsec=y*.0246,profile=mean,err_mean=np.sqrt(np.maximum(np.diag(cm),0)),scatter=np.sqrt(np.maximum(np.diag(scatter),0))))
        h=fits.BinTableHDU(table,name=f'PSF{number:03d}')
        h.header.update(dict(APERTURE=k[0],WAVEMIN=items[0]['lo'],WAVEMAX=items[0]['hi'],EPOCH=k[2],POSGROUP=k[3],SAMPLING=k[4],NSTARS=nstars,NEXP=len(allroots),NPROF=len(items),YMIN=min(r['trace_y'] for r in items),YMAX=max(r['trace_y'] for r in items),NORMRAD=k[5],GRADE='WING' if all(r['wing_usable'] for r in items) else 'CORE',SCATSTAT='BETWEEN_STARS' if nstars>1 else ('REPEAT_ONLY' if len(items)>1 else 'UNMEASURED')))
        hdus.extend([h,fits.ImageHDU(cm,name=f'MCOV{number:03d}'),fits.ImageHDU(pred,name=f'PCOV{number:03d}')])
        index.append(dict(extension=h.name,aperture=k[0],band=k[1],epoch=k[2],position_group=k[3],sampling=k[4],normalization_radius=k[5],stars=nstars,profiles=len(items),roots=allroots,quality=h.header['GRADE']))
        np.savez_compressed(f'products/combined_psfs/{h.name}.npz',y=y,profile=mean,covariance_mean=cm,covariance_prediction=pred,scatter=scatter)
    provenance=[{k:r.get(k,'') for k in ['root','target','program','date','aperture','trace_y','exptime','cal_ver','psf_status','processing_status']} for r in rows if r.get('processing_status')]
    used={root for item in index for root in item['roots']};catalog={r['root']:r for r in rows};manifest=load_json('data/catalogs/download_manifest.json',{})
    for r in provenance:
        r['library_used']=r['root'] in used;r['wavecorr']=catalog[r['root']].get('wavecorr','')
        r['configuration_aperture']=catalog[r['root']].get('configuration_aperture',r['aperture'])
        r['mast_obsid']=catalog[r['root']].get('obsid','')
        r['flt_uri']=next((uri for uri,p in manifest.items() if p.get('root')==r['root'] and p.get('path','').endswith('_flt.fits')),'')
    if provenance:hdus.append(fits.BinTableHDU(Table(rows=provenance),name='PROVENANCE'))
    hdus.append(fits.BinTableHDU(Table(rows=[dict(extension=r['extension'],aperture=r['aperture'],band=r['band'],epoch=r['epoch'],position_group=r['position_group'],sampling=r['sampling'],normalization_radius=r['normalization_radius'],n_stars=r['stars'],n_profiles=r['profiles'],n_exposures=len(r['roots']),roots=','.join(r['roots'])) for r in index]),name='INDEX'))
    Path('psf_library/combined').mkdir(parents=True,exist_ok=True)
    fits.HDUList(hdus).writeto('psf_library/stis_g140l_psf_library.fits',overwrite=True,checksum=True)
    import shutil
    shutil.copyfile('psf_library/stis_g140l_psf_library.fits','products/psf_library.fits')
    save_json('products/combined_psfs/index.json',index)
    return index


def interpolate_library(aperture,wavelength,epoch,position_group,sampling='hstpipe_native',normalization_radius=70):
    index=load_json('products/combined_psfs/index.json');group=[r for r in index if r['aperture']==aperture and r['epoch']==str(epoch) and r['position_group']==position_group and r['sampling']==sampling and r.get('normalization_radius',70)==normalization_radius and '-' in r['band']]
    group=[r for r in group if float(r['band'].split('-')[1])-float(r['band'].split('-')[0])<100]
    group.sort(key=lambda r:sum(map(float,r['band'].split('-')))/2)
    centers=np.array([sum(map(float,r['band'].split('-')))/2 for r in group])
    if len(centers)<2 or wavelength<centers[0] or wavelength>centers[-1]:raise ValueError('Unconstrained wavelength extrapolation prohibited')
    right=min(np.searchsorted(centers,wavelength,side='right'),len(centers)-1);left=max(0,right-1);f=(wavelength-centers[left])/(centers[right]-centers[left])
    a=np.load('products/combined_psfs/'+group[left]['extension']+'.npz');b=np.load('products/combined_psfs/'+group[right]['extension']+'.npz')
    # Unknown cross-wavelength covariance: conservative fully correlated bound.
    p=(1-f)*a['profile']+f*b['profile'];e=(1-f)*np.sqrt(np.maximum(np.diag(a['covariance_prediction']),0))+f*np.sqrt(np.maximum(np.diag(b['covariance_prediction']),0))
    return a['y'],p,e
