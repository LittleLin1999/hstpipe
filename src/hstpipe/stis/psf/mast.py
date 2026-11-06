"""Metadata-first discovery and immutable, checksummed, bounded downloads."""
from pathlib import Path
import logging,time,re
import numpy as np
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from astroquery.mast import Observations,MastMissions
from .common import *


def discover(c):
    p=Path('data/catalogs/all_recent_g140l.ecsv')
    if not p.exists():
        t=Observations.query_criteria(obs_collection='HST',instrument_name='STIS/FUV-MAMA',filters='G140L',dataRights='PUBLIC',t_min=[Time(c['date_start']).mjd,Time(c['date_end']).mjd+1])
        t.write(p)
    p=Path('data/catalogs/narrow_slit_mission.ecsv')
    if not p.exists():
        m=MastMissions(mission='hst')
        t=m.query_criteria(sci_instrume='STIS',sci_spec_1234='G140L',sci_aper_1234=c['apertures'][1:],sci_start_time='>='+c['date_start'],limit=5000)
        t.write(p)
    allobs=Table.read('data/catalogs/all_recent_g140l.ecsv')
    narrow=Table.read(p)
    names={str(r['sci_data_set_name']).lower() for r in narrow if str(r['sci_aper_1234'])=='52X0.5' and int(r['sci_pep_id']) in c.get('broad_stellar_programs',[17223]) and str(r['sci_targname']).upper().startswith(('B1.5IV','B2V','O7V','WD','GD','HZ','GRW','G191','FEIGE','AGK'))}
    previous={r['root']:r for r in load_json('data/catalogs/observations.json',[])}
    rows=[]
    for r in allobs:
        root=str(r['obs_id']).lower()
        selected=int(r['proposal_id']) in c['priority_programs']+c.get('additional_programs',[18113,17804]) or root in names
        rows.append(dict(obsid=str(r['obsid']),root=root,target=str(r['target_name']),program=int(r['proposal_id']),mjd=float(r['t_min']),date=Time(float(r['t_min']),format='mjd').isot,archive_exptime=float(r['t_exptime']),public=str(r['dataRights']),selected_for_preview=selected,status='PENDING_PREVIEW' if selected else 'NOT_PRIORITY_OR_NARROW_SLIT',reason=''))
    for r in rows:
        if r['root'] in previous:
            old=previous[r['root']];selection=r['selected_for_preview'];r.update(old);r['selected_for_preview']=selection
    save_json('data/catalogs/observations.json',rows)
    logging.info('Discovered %s G140L observations; %s previews',len(rows),sum(r['selected_for_preview'] for r in rows))
    return rows


def product_list(r):
    p=Path('data/catalogs/products')/(r['root']+'.ecsv');p.parent.mkdir(exist_ok=True)
    if not p.exists():Observations.get_product_list(r['obsid']).write(p)
    return Table.read(p)


def fetch(product,r):
    name=str(product['productFilename']);uri=str(product['dataURI'])
    base='data/raw' if name.endswith(('_raw.fits','_tag.fits')) else 'data/mast_calibrated'
    p=Path(base)/str(r['program'])/safe(r['target'])/name[:9]/name
    p.parent.mkdir(parents=True,exist_ok=True)
    manifest=load_json('data/catalogs/download_manifest.json',{})
    if uri in manifest and manifest[uri].get('status')=='COMPLETE' and p.exists():
        if p.stat().st_size==manifest[uri]['bytes'] and sha(p)==manifest[uri]['sha256']:return p
        raise RuntimeError('Immutable archive file changed: '+str(p))
    for attempt in range(3):
        try:
            tmp=p.with_suffix('.part')
            if tmp.exists():tmp.unlink()
            Observations.download_file(uri,local_path=str(tmp),cache=False)
            if not tmp.exists() or tmp.stat().st_size!=int(product['size']):raise IOError('Missing/truncated download')
            if name.endswith('.fits'):
                with fits.open(tmp,memmap=False) as h:
                    h.verify('exception')
                    for hd in h:
                        _=hd.data
                        if 'CHECKSUM' in hd.header and hd.verify_checksum()!=1:raise IOError('CHECKSUM mismatch')
            tmp.replace(p);p.chmod(0o444)
            manifest[uri]=dict(status='COMPLETE',path=str(p),url='https://mast.stsci.edu/api/v0.1/Download/file?uri='+uri,bytes=p.stat().st_size,sha256=sha(p),obsid=r['obsid'],root=r['root'],attempts=attempt+1)
            save_json('data/catalogs/download_manifest.json',manifest)
            return p
        except Exception as e:
            logging.warning('Download %s attempt %d: %s',name,attempt+1,e)
            manifest[uri]=dict(status='FAILED',error=str(e),path=str(p),attempts=attempt+1)
            save_json('data/catalogs/download_manifest.json',manifest)
            time.sleep(2**attempt)
    return None


def preview(c):
    rows=load_json('data/catalogs/observations.json')
    for r in sorted(rows,key=lambda r: ([18169,17644,17649,18229,18113,16966,17804].index(r['program']) if r['program'] in [18169,17644,17649,18229,18113,16966,17804] else 99,r['root'])):
        if not r['selected_for_preview'] or r.get('preview_done'):continue
        try:
            products=product_list(r)
            # Preserve exact product metadata; only dataset-matched small products at this stage.
            pp=[p for p in products if str(p['productFilename']).startswith(r['root']+'_')]
            x=[p for p in pp if str(p['productFilename']).endswith('_x1d.fits')]
            if not x:
                r.update(status='NO_X1D',reason='No dataset-matched X1D available; calibration needed')
                if r['program'] in c['priority_programs']:
                    raw=[p for p in pp if str(p['productFilename']).endswith('_raw.fits')]
                    if raw:
                        path=fetch(raw[0],r)
                        if path:
                            with fits.open(path) as h:
                                r['raw_header']={k:h[0].header.get(k) for k in ['INSTRUME','DETECTOR','OPT_ELEM','APERTURE','CENWAVE','OBSTYPE','OBSMODE','TARGNAME','SCLAMP','DARKCORR','X1DCORR','X2DCORR']}
                            r['paths']={'raw':str(path)}
            else:
                path=fetch(x[0],r)
                if path:
                    with fits.open(path) as h:
                        hdr=h[0].header
                        r.update({k.lower():hdr.get(k) for k in ['INSTRUME','DETECTOR','OPT_ELEM','APERTURE','PROPAPER','CENWAVE','OBSMODE','CAL_VER','POSTARG1','POSTARG2']})
                        r['x1d']=str(path)
                        sn=[]
                        for hd in h[1:]:
                            if hd.data is None or not hd.data.dtype.names:continue
                            for spec in hd.data:
                                w=spec['WAVELENGTH'];f=spec['FLUX'];e=spec['ERROR']
                                mask=(w>1330)&(w<1680)&np.isfinite(f)&(e>0)
                                sn.append(float(np.sum(f[mask])/np.sqrt(np.sum(e[mask]**2))))
                        r['preview_snr']=max(sn,default=0)
                    correct=r['instrume']=='STIS' and r['detector']=='FUV-MAMA' and r['opt_elem']=='G140L' and r['aperture'] in c['apertures'] and r['cenwave']==1425
                    if not correct:r.update(status='REJECT_CONFIGURATION',reason='FITS configuration does not match selected modes')
                    elif r['preview_snr']<c['preview_snr']:r.update(status='REJECT_PREVIEW_SNR',reason='Broad continuum preview S/N below threshold')
                    else:r.update(status='SELECTED_FOR_PIXELS',reason='Preview only; 2D QA pending')
                else:r.update(status='DOWNLOAD_FAILED',reason='X1D download failed')
            r['preview_done']=True
        except Exception as e:r.update(status='PREVIEW_ERROR',reason=str(e));logging.exception('Preview %s',r['root'])
        save_json('data/catalogs/observations.json',rows)
        logging.info('Preview %s %s SNR=%s',r['root'],r['status'],r.get('preview_snr'))
    return rows


def download(c):
    rows=load_json('data/catalogs/observations.json')
    for r in sorted(rows,key=lambda r: ([18169,17644,17649,18229,18113,16966,17804].index(r['program']) if r['program'] in [18169,17644,17649,18229,18113,16966,17804] else 99,r['root'])):
        if not r['selected_for_preview']:continue
        if r['status']!='SELECTED_FOR_PIXELS' and not (r['program']==18113 and r.get('x1d')):continue
        products=product_list(r);paths={};available=[]
        for p in products:
            name=str(p['productFilename'])
            if not name.startswith(r['root']+'_') or not name.endswith(('.fits','.txt')):continue
            suffix=name.split('_',1)[1].split('.')[0]
            if suffix in ['raw','flt','x2d','x1d','spt','wav','wsp','trl','tag','asn','jit','jif']:
                available.append(suffix);path=fetch(p,r)
                if path:paths[suffix]=str(path)
        r['paths']=paths;r['missing_products']=sorted(set(['raw','flt','x2d','x1d','spt','wav','trl'])-set(available))
        save_json('data/catalogs/observations.json',rows)
        logging.info('Downloaded %s %s',r['root'],list(paths))
