import json
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.target_runner import save
from hstpipe.uvis.reporting import generate_report
from hstpipe.uvis.reporting.decisions import decision


def test_report_is_a_brief_partial_success_summary_with_closed_details(tmp_path):
    p=ReductionPaths(tmp_path,'TEST','uvis','F275W');product=p.product('drc');product.parent.mkdir(parents=True);product.write_bytes(b'preserved')
    ready=dict(filter='F275W',status='complete',products=[{'path':str(product)}],scientific_assessment=dict(status='science_ready',relative_astrometry=dict(status='validated'),astrometry=dict(status='validated',reference='gaia_dr3',final_validation=dict(rms_arcsec=.01))))
    blocked=dict(filter='F621M',status='relative_astrometry_unvalidated',scientific_assessment=dict(relative_astrometry=dict(status='unvalidated',original=dict(pairs=[dict(postfit_scatter_pix=1.5,scatter_limit_pix=.4)]))))
    original=dict(target='TEST',filters={'F275W':ready,'F621M':blocked});save(p.directory/'manifest.json',original)
    report=generate_report(p.directory);text=report.read_text()
    assert '1 of 2 filters science-ready' in text
    assert '1.50 pixel scatter' in text and '0.40 pixel limit' in text
    assert 'Applied corrections and current WCS' in text and 'Technical details' in text
    assert '<details open' not in text
    assert 'No relative correction accepted; original WCS retained.' in text
    assert json.loads((p.directory/'manifest.json').read_text())['filters']==original['filters']
    assert product.read_bytes()==b'preserved' and report.stat().st_size<1_000_000


def test_unverified_low_rms_or_local_solution_not_called_ready():
    for status in ('unvalidated','local_only'):
        record=dict(status='relative_astrometry_unvalidated',scientific_assessment=dict(relative_astrometry=dict(status=status,postfit_scatter_pix=.01)))
        assert decision(record)[0]=='Blocked'


def test_failed_aperture_and_source_are_visible_before_collapsed_details(tmp_path):
    p=ReductionPaths(tmp_path,'TEST','uvis','F336W')
    source=dict(source_id='G123',status='review_required',output_x=40.5,output_y=60.25,
        failures=['significant_persistent_source_flux_loss'],observations=[],
        apertures=[dict(radius_pix=2.5,status='review_required',loss_fraction=.066,loss_significance=5.79),
                   dict(radius_pix=5.,status='consistent',loss_fraction=.049,loss_significance=1.94)])
    record=dict(filter='F336W',status='image_quality_review_required',products=[],scientific_assessment=dict(
        status='review_required',relative_astrometry=dict(status='validated'),astrometry=dict(status='validated'),
        unvalidated_final='/example/TEST_F336W_unvalidated_final.fits',
        source_validation=dict(status='review_required',validated_sources=9,failed_sources=1,sources=[source])))
    save(p.directory/'manifest.json',dict(target='TEST',filters={'F336W':record}))
    report=generate_report(p.directory);text=report.read_text();visible=text.split('<details',1)[0]
    assert 'Human review required' in visible and 'G123' in visible
    assert 'TEST_F336W_unvalidated_final.fits' in visible and '41.50, 61.25' in visible
    assert '2.5-pixel radius: 6.60% lower flux; 5.79σ' in visible
    assert 'agreeing within uncertainty: 5.0 pixels' in visible
    assert 'Final / reference flux (5-pixel radius)' not in text
    assert 'Aperture radius (pixels)' in text
    assert json.loads((p.directory/'manifest.json').read_text())['filters']['F336W']==record
    assert report.stat().st_size<1_000_000


def test_review_figure_is_embedded_and_saved_pixels_unchanged(tmp_path):
    import numpy as np
    from astropy.io import fits
    from hstpipe.uvis.reporting.report import _review_figures
    y,x=np.indices((77,77));star=np.exp(-((x-38.2)**2+(y-37.8)**2)/4)
    hd=[fits.PrimaryHDU()]
    for name,offset in [('INPUT',1.),('INPUT',2.),('MEDIAN',0.),('FINAL',1.),('WHT',0.)]:
        a=np.ones_like(star)*100 if name=='WHT' else star+offset
        h=fits.ImageHDU(a,name=name);h.header['SOURCEID']='Gtest';h.header['EXPOSURE']='exposure_flc.fits';hd.append(h)
    path=tmp_path/'cutouts.fits';fits.HDUList(hd).writeto(path);original=path.read_bytes()
    source=dict(source_id='Gtest',status='review_required',output_x=38.2,output_y=37.8,apertures=[dict(radius_pix=2.5,status='review_required')])
    records={'F336W':dict(scientific_assessment=dict(source_validation=dict(cutouts=str(path),sources=[source])))}
    figures,count=_review_figures(records,220_000);fragment=figures['F336W','Gtest']
    assert count==1 and 'data:image/jpeg;base64,' in fragment
    assert '2 of 2 individual exposures' in fragment and 'not raw detector panels' in fragment
    assert 'same brightness stretch' in fragment and 'not an error-significance map' in fragment
    assert len(fragment.encode())<=220_000 and path.read_bytes()==original
    limited,n=_review_figures(records,1)
    assert n==0 and 'budget exhausted' in limited['F336W','Gtest']


def test_review_missing_cutouts_remain_explicit():
    from hstpipe.uvis.reporting.report import _review_figures
    records={'F336W':dict(scientific_assessment=dict(source_validation=dict(cutouts='/missing/qa.fits',sources=[dict(source_id='Gtest',status='review_required')])))}
    figures,count=_review_figures(records,220_000)
    assert count==0 and 'Review figure unavailable' in figures['F336W','Gtest']
