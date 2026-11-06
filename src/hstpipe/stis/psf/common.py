from pathlib import Path
import json, hashlib, logging
import numpy as np
import yaml

def config(path='config/psf.yaml'):
    """Read PSF settings through hstpipe's standard ReductionConfig schema."""
    from hstpipe.common.config import ReductionConfig
    p=Path(path).resolve();base=ReductionConfig.read(p)
    if base.instrument!='stis':raise ValueError('PSF calibration requires instrument: stis')
    settings=dict(base.options.get('psf',{}));c={}
    selection=ReductionConfig.read(p.parent/settings.pop('selection_config','selection.yaml'))
    reduction=ReductionConfig.read(p.parent/settings.pop('reduction_config','reduction.yaml'))
    c.update(selection.options.get('psf_selection',{}));c.update(settings)
    c['workdir']=str(Path(base.output).parent)
    c['reference_dir']=reduction.options['reference_dir']
    c['reduction_options']=dict(reduction.options)
    c['hstpipe_spatial_scale']=reduction.options.get('spatial_scale',.0123)
    c['hstpipe_pixfrac']=reduction.options.get('pixfrac',1.)
    return c
def save_json(path, value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    p.with_suffix(p.suffix+'.tmp').write_text(json.dumps(value,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
    p.with_suffix(p.suffix+'.tmp').replace(p)
def load_json(path, default=None):
    value=json.loads(Path(path).read_text()) if Path(path).exists() else default
    if str(path)=='data/catalogs/observations.json' and isinstance(value,list):
        q=Path('products/measurements/exposure_qa.json')
        if q.exists():
            byroot={r['root']:r for r in json.loads(q.read_text())}
            for r in value:
                if r['root'] in byroot:r.update(byroot[r['root']])
    return value

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def safe(s):
    return ''.join(c if c.isalnum() or c in '-+.' else '_' for c in str(s))
def setup():
    for p in ['data/catalogs','data/raw','data/mast_calibrated','data/recalibrated','products/individual_profiles','products/combined_psfs','products/measurements','figures','reports','logs']:
        Path(p).mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.StreamHandler(),logging.FileHandler('logs/pipeline.log')])


def aperture_configuration(header):
    """Requested pseudo-aperture, kept distinct from the hardware APERTURE."""
    return str(header.get('PROPAPER') or header['APERTURE']).strip()


def catalog_apertures(c):
    from astropy.io import fits
    rows=load_json('data/catalogs/observations.json')
    for r in rows:
        path=r.get('paths',{}).get('flt') or r.get('paths',{}).get('raw')
        if not path:continue
        h=fits.getheader(path,0)
        r.update(physical_aperture=h.get('APERTURE'),propaper=h.get('PROPAPER'),configuration_aperture=aperture_configuration(h))
    save_json('data/catalogs/observations.json',rows)
    return rows
