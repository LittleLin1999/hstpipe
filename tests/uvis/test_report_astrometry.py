"""Human-readable fit membership, proposed/applied offsets and WCS provenance."""
import json
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.target_runner import save
from hstpipe.common.provenance import identity
from hstpipe.uvis.reporting import generate_report


@pytest.mark.parametrize('applied',[False,True])
def test_report_lists_measured_sources_and_distinguishes_applied_WCS(flc,tmp_path,applied):
 raw=flc('a');paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F621M')
 work=paths.exposure_file('intermediate','alignment','a','flc');paths.copy(raw,work)
 with fits.open(work,mode='update') as hd:hd['SCI',7].header['WCSNAME']='HP-RELATIVE' if applied else 'ORIGINAL'
 before=[identity(raw),identity(work)]
 rows=[dict(source_id=sid,exposure=work.name,classification=kind,catalog_ra=ra,catalog_dec=20.,status='confirmed_common',reason='Measured native HST centroid',extver=7,x=x,y=64.,centroid_error_pix=.08) for sid,kind,ra,x in [('G123','gaia_stellar_candidate',10.,30.),('L456','legacy_point_candidate',10.001,60.)]]
 paths.area('intermediate','alignment',create=True)
 Table(rows=rows).write(paths.file('intermediate','alignment','common_sources','ecsv'),format='ascii.ecsv')
 f=dict(exposure=work.name,source_ids=['G123','L456'],matched=2,geometry='shift',shift=[.2,-.1],dx_uncertainty_pix=.08,dy_uncertainty_pix=.09,rotation_measured=False,validation_scope='source_constrained')
 joint=dict(may_combine=True,corrections_applied=applied,fits=[dict(f,exposure='reference_flc.fits'),f],measurement_inlier_mask=[True]*4,distinct_sources=2,residual_rms_pix=.03,degrees_of_freedom=2,reason='Relative fit accepted')
 relative=dict(status='validated',joint_models=joint,wcs_solutions=[dict(exposure=work.name,chips=[dict(extver=7,primary='ORIGINAL')])])
 record=dict(filter='F621M',status='astrometry_unvalidated',inputs=[{'path':str(raw)}],observing_configurations=[dict(ROOTNAME='a',RA_TARG=10.001,DEC_TARG=20.)],scientific_assessment=dict(stages=dict(cross_filter_alignment=dict(status='anchor')),relative_astrometry=relative,astrometry=dict(status='unvalidated'),multi_chip=[dict(exposure=work.name,rootname='a',chips=[dict(extver=7,ccdchip=2,wcsname='HP-RELATIVE' if applied else 'ORIGINAL')])]))
 save(paths.directory/'manifest.json',dict(target='TEST',filters={'F621M':record}))
 report=generate_report(paths.directory);text=report.read_text()
 assert 'Objects used for astrometry' in text and 'G123' in text and 'L456' in text
 assert 'Gaia stellar candidate' in text and 'near science target' in text
 assert '30.000 / 64.000' in text and '0.2000 ± 0.0800' in text
 assert 'No validated absolute WCS' in text and 'Fixed; not measured' in text
 assert 'Which WCS is in use?' in text and 'Numerical working WCS' in text
 assert ('Relative correction applied' if applied else 'Not applied; existing WCS retained') in text
 assert ('0.2000 / -0.1000' if applied else '0 / 0 (retained)') in text
 assert 'NOT_RUN — anchor' in text and 'FAIL — anchor' not in text
 assert '<svg' in text and len(report.read_bytes())<1_000_000
 assert identity(raw)==before[0] and identity(work)==before[1]


def test_failed_report_does_not_invent_sources_or_corrected_WCS(tmp_path):
 paths=ReductionPaths(tmp_path,'T','uvis','F275W')
 save(paths.directory/'manifest.json',dict(target='T',filters={'F275W':dict(filter='F275W',status='relative_astrometry_unvalidated',scientific_assessment=dict(alignment=dict(status='unvalidated')))}))
 text=generate_report(paths.directory).read_text()
 assert 'No saved object-level fit identities' in text
 assert 'No accepted correction is implied' in text and 'No validated absolute WCS' in text
 assert 'HP-RELATIVE' not in text and 'Gaia-aligned</strong>' not in text


def test_accepted_existing_wcs_is_not_hidden_by_rejected_joint_proposal(tmp_path):
 from hstpipe.uvis.report_astrometry import context,render
 paths=ReductionPaths(tmp_path,'T','uvis','F336W');paths.area('intermediate','alignment',create=True)
 cats={}
 for name in ('r_flc.fits','a_flc.fits'):
  p=paths.file('intermediate','alignment',name.removesuffix('.fits')+'_identified','ecsv')
  Table(dict(common_source_id=['Gstar','HST-object'])).write(p,format='ascii.ecsv');cats[name]=str(p)
 Table(dict(reference_index=[0,1],source_index=[0,1])).write(paths.file('intermediate','alignment','original_a_flc_matches','ecsv'),format='ascii.ecsv')
 Table(rows=[dict(source_id='Gstar',exposure='a_flc.fits',classification='gaia_stellar_candidate',status='confirmed_common',extver=1,x=20.,y=30.,centroid_error_pix=.08)]).write(paths.file('intermediate','alignment','common_sources','ecsv'),format='ascii.ecsv')
 fit=dict(exposure='a_flc.fits',source_ids=['Gstar','HST-object'],shift=[.2,-.1],geometry='shift',matched=2)
 joint=dict(status='unvalidated',may_combine=False,corrections_applied=False,fits=[dict(fit,exposure='r_flc.fits'),fit])
 relative=dict(status='validated',method='existing_wcs',joint_models=joint,detector_catalogs=[dict(exposure=n,path=p) for n,p in cats.items()],pairs=[dict(status='validated',reference='r_flc.fits',exposure='a_flc.fits',fit_inlier_mask=[True,False])])
 record=dict(scientific_assessment=dict(relative_astrometry=relative));ctx=context(paths,record);text=render(ctx,record)
 assert ctx['used'][('a_flc.fits','Gstar')] and not ctx['used'][('a_flc.fits','HST-object')]
 assert 'Used in accepted relative fit / frame validation' in text
 assert 'joint correction below was proposed and rejected' in text
 assert '0 / 0 (retained)' in text and 'Not applied; existing WCS retained' in text
