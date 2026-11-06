import numpy as np
import pytest
from astropy.io import fits
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry.relative_fit import quality,match

@pytest.fixture
def options():return settings({})['astrometry']['relative']

def grid(n=30):return np.random.default_rng(7).uniform(0,1000,(n,2))

def test_initial_offsets_are_distinct_from_fit_residuals(options):
    a=grid();q=quality(a,a+[4,2],np.full(len(a),.1),options)
    assert q['status']=='refinement_supported' and q['initial_rms_pix']>4
    assert q['postfit_scatter_pix']<1e-8 and q['loo_rms_pix']<1e-8
    assert quality(a,a,np.full(len(a),.1),options)['status']=='validated'

def test_histogram_finds_offset_outside_narrow_tolerance(options):
    a=grid();i,j=match(a,a+[12,-7],options,True)
    assert len(i)==len(a) and np.array_equal(i,j)

def test_geometry_counts_and_degeneracy(options):
    a=grid(4)
    assert quality(a,a,np.ones(4)*.1,options,'shift')['status']=='validated'
    assert quality(a,a,np.ones(4)*.1,options,'rshift')['status']=='validated'
    assert quality(a,a,np.ones(4)*.1,options,'rscale')['status']=='validated'
    a[:,1]=0
    assert quality(a,a,np.ones(4)*.1,options,'rshift')['status']=='validated'
    assert quality(a,a,np.ones(4)*.1,options,'general')['status']=='degenerate_distribution'

def test_unreliable_fits_cannot_be_accepted(options):
    a=grid();b=a+np.random.default_rng(9).normal(0,4,a.shape)
    assert quality(a,b,np.ones(len(a))*.1,options)['status']=='unreliable_fit'
    b[0,0]=np.nan
    assert quality(a,b,np.ones(len(a))*.1,options)['status']=='nonfinite'

def test_same_visit_unresolved_registration_blocks_combination(flc,tmp_path,monkeypatch,options):
    from hstpipe.uvis.astrometry import relative
    from hstpipe.common.paths import ReductionPaths
    from astropy.table import Table
    files=[flc('a'),flc('b')]
    cat=Table(names=['x','y','ra','dec','flux','extver'],dtype=['f8']*5+['i8'])
    monkeypatch.setattr(relative,'detect_sources',lambda *a:cat)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True)
    r=relative.validate_relative(files,options,paths)
    assert not r['may_combine']

def test_absolute_and_tweakback_never_run_with_unvalidated_relative(monkeypatch):
    from hstpipe.uvis.astrometry import workflow
    monkeypatch.setattr(workflow,'absolute_alignment',lambda *a:pytest.fail('Absolute calibration must remain blocked'))
    products,report=workflow.finalize(None,[],[],None,{}, {},dict(status='unvalidated'))
    assert products==[] and report['status']=='relative_astrometry_unvalidated'

def test_segmentation_rejects_isolated_hot_pixel(flc,options):
    from hstpipe.uvis.astrometry.sources import sources
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:hd['SCI',7].data[30,30]=1e6
    c=sources(p,options['detection'],True,6)
    assert len(c)>0 and np.min(np.hypot(c['x']-30,c['y']-30))>5

def test_validation_cannot_be_disabled():
    with pytest.raises(ValueError,match='mandatory'):settings({'astrometry':{'relative':{'validate':False}}})

@pytest.mark.parametrize('accept',[True,False])
def test_refinement_transaction_requires_image_quality(flc,tmp_path,monkeypatch,options,accept):
    from astropy.table import Table
    from astropy.wcs import WCS
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.common.provenance import identity
    from hstpipe.uvis.astrometry import relative_engine as engine
    files=[flc('first',size=256),flc('second',size=256)]
    with fits.open(files[1],mode='update') as hd:hd['SCI',7].header['CRPIX1']+=5
    before=[identity(p) for p in files]
    def catalog(p,*args,**kwargs):
        xy=np.array([[20,20],[100,20],[230,20],[20,100],[230,100],[20,230],[100,230],[230,230]])
        with fits.open(p) as hd:world=WCS(hd['SCI'].header,hd).all_pix2world(xy,0)
        return Table(dict(x=xy[:,0],y=xy[:,1],ra=world[:,0],dec=world[:,1],flux=np.full(8,100),extver=np.full(8,7),centroid_error_pix=np.full(8,.1),width_pix=np.ones(8)))
    monkeypatch.setattr(engine,'sources',catalog)
    def tweak(inputs,*args,**kwargs):
        with fits.open(inputs[0]) as hd:value=hd['SCI',7].header['CRPIX1']
        with fits.open(inputs[1],mode='update') as hd:hd['SCI',7].header['CRPIX1']=value
    monkeypatch.setattr(engine,'tweak',tweak)
    def mosaic(inputs,prefix,*a):
        dest=prefix.with_suffix('.fits');fits.PrimaryHDU().writeto(dest,overwrite=True);return dest
    monkeypatch.setattr(engine,'internal_mosaic',mosaic)
    monkeypatch.setattr(engine,'image_quality',lambda *a,**kw:dict(status='validated' if accept else 'image_quality_unvalidated'))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    for a,s in [('intermediate','alignment'),('intermediate','drizzle'),('diagnostics','alignment')]:paths.area(a,s,create=True)
    options=dict(options,methods=['existing_wcs','tweakreg'])
    report=engine.refine_relative(files,paths,options)
    assert report['may_combine']==accept
    assert identity(files[0])==before[0]
    with fits.open(files[1]) as hd:assert hd['SCI',7].header['CRPIX1']==(128 if accept else 133)
    if not accept:assert [identity(p) for p in files]==before

@pytest.mark.parametrize('after_width,accepted',[(1.8,True),(3.,False)])
def test_no_cr_coadd_shape_comparison(tmp_path,options,after_width,accepted):
    from astropy.table import Table
    from astropy.wcs import WCS
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry.relative_engine import image_quality
    w=WCS(naxis=2);w.wcs.crpix=[128,128];w.wcs.crval=[10,20];w.wcs.cdelt=[-.04/3600,.04/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    xy=np.array([[50,50],[200,50],[50,200],[200,200]])
    yy,xx=np.indices((256,256));files=[]
    for name,width in [('before',2.),('after',after_width)]:
        data=sum(np.exp(-((xx-x)**2+(yy-y)**2)/(2*width**2))/width**2 for x,y in xy)
        p=tmp_path/(name+'.fits');fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(data,header=w.to_header(),name='SCI'),fits.ImageHDU(np.ones(data.shape),name='WHT')]).writeto(p);files.append(p)
    world=w.all_pix2world(xy,0);catalog=Table(dict(ra=world[:,0],dec=world[:,1]))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('diagnostics','alignment',create=True)
    q=image_quality(*files,catalog,options,paths,'comparison')
    assert (q['status']=='validated')==accepted and q['no_cr_rejection']

@pytest.mark.parametrize('line,valid', [('0 0 0 1 .1 .1',True),('nan 0 0 1 0 0',False),('0 0 1 1 .1 .1',False)])
def test_native_custom_reference_and_transform_checks(flc,tmp_path,monkeypatch,options,line,valid):
    from astropy.table import Table
    from hstpipe.uvis.astrometry.relative_engine import tweak
    from drizzlepac import tweakreg
    p=flc();work=tmp_path/'work';work.mkdir()
    catalog=Table(dict(x=[20,30,40,50],y=[20,40,60,80],ra=[10.]*4,dec=[20.]*4,flux=[100.]*4,extver=[7]*4))
    from hstpipe.uvis.astrometry.centroids import QUALITY_VERSION
    catalog['centroid_valid']=np.ones(len(catalog),bool);catalog.meta['centroid_quality']=dict(version=QUALITY_VERSION)
    reference=catalog.copy();reference['ra']=[11.]*4
    def native(files,**kw):
        assert np.allclose(np.loadtxt(kw['refcat'])[:,0],11)
        assert 'refimage' not in kw and kw['refxyunits']=='degrees'
        (work/'relative_shifts.txt').write_text(str(p)+' '+line+'\n')
    monkeypatch.setattr(tweakreg,'TweakReg',native)
    if valid:tweak([p],[catalog],work,options,'shift',True,reference)
    else:
        with pytest.raises(ValueError):tweak([p],[catalog],work,options,'shift',True,reference)


def test_better_catalogs_can_validate_original_without_any_correction(flc,tmp_path,monkeypatch,options):
    from astropy.table import Table
    from astropy.wcs import WCS
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.common.provenance import identity
    from hstpipe.uvis.astrometry import relative_engine as engine
    files=[flc('one',size=256),flc('two',size=256)];before=[identity(p) for p in files]
    def catalog(p,options,extended=False,width=3.5,**kwargs):
        xy=np.array([[20,20],[100,20],[230,20],[20,100],[230,100],[20,230],[100,230],[230,230]])
        if not extended:xy=xy[:1]
        with fits.open(p) as hd:world=WCS(hd['SCI'].header,hd).all_pix2world(xy,0)
        return Table(dict(x=xy[:,0],y=xy[:,1],ra=world[:,0],dec=world[:,1],flux=np.full(len(xy),100),extver=np.full(len(xy),7),centroid_error_pix=np.full(len(xy),.1)))
    monkeypatch.setattr(engine,'sources',catalog)
    monkeypatch.setattr(engine,'tweak',lambda *a,**k:pytest.fail('Noisy refinement must not replace a validated original'))
    monkeypatch.setattr(engine,'internal_mosaic',lambda *a,**k:pytest.fail('Original solution needs no candidate mosaic'))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    for a,s in [('intermediate','alignment'),('intermediate','drizzle'),('diagnostics','alignment')]:paths.area(a,s,create=True)
    report=engine.refine_relative(files,paths,dict(options,methods=['improved_catalogs']))
    assert report['status']=='validated' and report['method']=='existing_wcs_improved_catalogs'
    assert report['pairs'][0]['matched']==8
    assert [identity(p) for p in files]==before


def test_requested_user_initialization_remains_supported(flc,tmp_path,monkeypatch):
    from hstpipe.uvis import pipeline
    raw=tmp_path/'raw';raw.mkdir();p=flc();p.rename(raw/p.name)
    def user(*args,**kw):
        assert args[3]['mode']=='user'
        raise RuntimeError('requested initialization reached')
    monkeypatch.setattr(pipeline,'register',user)
    outputs=pipeline.reduce_uvis(dict(input=str(raw),output=str(tmp_path/'redux'),instrument='uvis',options={'alignment':{'mode':'user','transformations':{'test':{'dx':0,'dy':0}}}}))
    import json
    state=json.loads((outputs[0]/'manifest.json').read_text())['filters']['F275W']
    assert state['status']=='failed' and 'requested initialization reached' in state['error']
    assert not state.get('products')


@pytest.mark.parametrize('bit,fresh',[(16,False),(256,False),(2048,False),(8192,False),(4096,True)])
def test_segment_centroid_inside_masked_core_is_rejected(flc,options,bit,fresh):
    from hstpipe.uvis.astrometry.sources import sources
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:hd['DQ',7].data[127:130,127:130]=bit;hd[0].header['HPCRRUN']=fresh
    catalog=sources(p,options['detection'],True,6)
    assert not np.any(np.hypot(catalog['x']-128,catalog['y']-128)<4)


@pytest.mark.parametrize('case',['noise','weight_deficit'])
def test_image_quality_cannot_validate_noise_or_masked_source_cores(tmp_path,options,case):
    from astropy.table import Table
    from astropy.wcs import WCS
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry.relative_engine import image_quality
    w=WCS(naxis=2);w.wcs.crpix=[128,128];w.wcs.crval=[10,20];w.wcs.cdelt=[-.04/3600,.04/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    xy=np.array([[50,50],[200,50],[50,200],[200,200]])
    yy,xx=np.indices((256,256));data=np.random.default_rng(42).normal(0,1,(256,256));weight=np.ones(data.shape)
    if case=='weight_deficit':
        for x,y in xy:
            data+=100*np.exp(-((xx-x)**2+(yy-y)**2)/8)
            weight[(xx-x)**2+(yy-y)**2<9]=.1
    files=[]
    for name in ('before','after'):
        p=tmp_path/(name+'.fits');fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(data,header=w.to_header(),name='SCI'),fits.ImageHDU(weight,name='WHT')]).writeto(p);files.append(p)
    world=w.all_pix2world(xy,0);catalog=Table(dict(ra=world[:,0],dec=world[:,1]))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('diagnostics','alignment',create=True)
    qa=image_quality(*files,catalog,options,paths,'quality')
    assert qa['status']=='image_quality_unvalidated' and qa['before']['sources']==0
