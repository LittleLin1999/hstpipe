import json
import shutil
import numpy as np
from astropy.io import fits
from astropy.table import Table
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.common.target_runner import save
from hstpipe.uvis.reporting import generate_report
from hstpipe.uvis.report_alignment import exposure_cutouts


def dataset(flc,tmp_path):
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F621M')
    raw=[flc(name='a',size=128),flc(name='b',size=128)];working=[];catalogs=[]
    for i,p in enumerate(raw):
        with fits.open(p,mode='update') as hd:
            if i:hd['SCI',7].header['CRPIX1']+=8
            yy,xx=np.indices((128,128));hd['SCI',7].data[:]=sum(100*np.exp(-((xx-x)**2+(yy-y)**2)/4) for x,y in [(44,64),(84,64),(64,94)])+1
        q=paths.exposure_file('intermediate','working',p.stem.removesuffix('_flc'),'flc')
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
        if i:
            with fits.open(q,mode='update') as hd:hd['SCI',7].header['CRPIX1']-=8
        working.append(q)
        cat=Table(dict(x=[44.,84.,64.],y=[64.,64.,94.],extver=[7]*3,
                       common_source_id=['G0','G1','G2'],centroid_error_pix=[.08]*3,centroid_valid=[True]*3))
        dest=paths.file('intermediate','alignment',p.stem+'_sources','ecsv');dest.parent.mkdir(parents=True,exist_ok=True)
        cat.write(dest,format='ascii.ecsv');catalogs.append(dict(exposure=p.name,path=str(dest)))
    record=dict(filter='F621M',status='complete',inputs=[identity(p) for p in raw],
                scientific_assessment=dict(alignment=dict(status='validated',detector_catalogs=catalogs),tweakback=dict(inputs=[identity(p) for p in working])))
    return paths,record,raw+working


def test_native_exposures_use_one_fixed_reference_and_preserve_original_data(flc,tmp_path):
    paths,r,files=dataset(flc,tmp_path);before=[identity(p) for p in files]
    result=exposure_cutouts(paths,r)
    assert len(result['figures'])==3
    for figure in result['figures']:
        a,b=figure['measurements']
        assert np.linalg.norm(a['before']['offset_mas'])<.001
        assert np.linalg.norm(b['before']['offset_mas'])>300
        assert np.linalg.norm(b['after']['offset_mas'])<.001
        assert b['after']['extver']==7 and b['after']['ccdchip']==2
    assert [identity(p) for p in files]==before
    assert exposure_cutouts(paths,r)==result


def test_report_contains_individual_exposure_reference_crosses_and_honest_status(flc,tmp_path):
    paths,r,files=dataset(flc,tmp_path)
    state=dict(target='TEST',filters={'F621M':r});save(paths.directory/'manifest.json',state)
    output=generate_report(paths.directory);text=output.read_text()
    assert 'Before vs after alignment' in text and 'One row per exposure.' in text
    assert '3 sources, individual exposures' in text
    assert 'individual detector images, not coadds' in text and 'cyan cross' in text
    assert text.count('class="alignment-pair"')==3
    assert 'data:image/jpeg;base64,' in text and output.stat().st_size<1_000_000
    assert json.loads((paths.directory/'manifest.json').read_text())['filters']==state['filters']
    r['scientific_assessment']['alignment']['status']='unvalidated'
    save(paths.directory/'manifest.json',state)
    assert 'No validated relative correction' in generate_report(paths.directory).read_text()


def test_stale_relative_working_wcs_is_not_presented_as_accepted(flc,tmp_path):
    paths,r,files=dataset(flc,tmp_path)
    with fits.open(files[-1],mode='update') as hd:hd['SCI',7].header['CRVAL1']+=.01
    save(paths.directory/'manifest.json',dict(target='TEST',filters={'F621M':r}))
    text=generate_report(paths.directory).read_text()
    assert 'checkpoint has changed' in text and 'class="alignment-pair"' not in text
