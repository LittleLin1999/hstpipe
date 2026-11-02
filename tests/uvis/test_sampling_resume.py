"""Final sampling changes preserve validated WCS; no network or drizzle jobs."""
from contextlib import contextmanager
import json
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.common.target_runner import save
from hstpipe.uvis.astrometry import workflow
from hstpipe.uvis.config import settings

@pytest.mark.parametrize('tamper',[False,True])
def test_sampling_resume_reuses_only_intact_absolute_checkpoint(flc,tmp_path,monkeypatch,tamper):
    p=flc();paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    paths.area('intermediate','alignment',create=True);paths.area('intermediate','drizzle',create=True)
    with fits.open(p,mode='update') as hd:hd['SCI',7].header['WCSNAME']='ACCEPTED'
    catalog=paths.file('intermediate','alignment','reference','ecsv');Table(dict(ra=[1.],dec=[2.])).write(catalog,format='ascii.ecsv')
    solution=paths.file('intermediate','alignment','absolute');save(solution,dict(status='validated',corrected_wcs_name='ACCEPTED',reference_catalog=str(catalog),reference_mode='cross_filter'))
    save(paths.file('intermediate','drizzle','provisional_inputs'),[identity(p)])
    paths.cache().record('uvis_absolute','accepted',[p,catalog,solution])
    if tamper:catalog.write_text(catalog.read_text()+'# changed\n')
    def forbidden(*a,**k):raise AssertionError('Must not re-query or refit absolute WCS')
    monkeypatch.setattr(workflow,'absolute_alignment',forbidden)
    class FinalStageReached(Exception):pass
    @contextmanager
    def scratch(*a,**k):raise FinalStageReached;yield
    monkeypatch.setattr(ReductionPaths,'scratch',scratch)
    cfg=settings({});before=identity(p)
    with pytest.raises(ValueError if tamper else FinalStageReached):
        workflow.finalize(p,[p],[],paths,cfg['astrometry'],cfg['drizzle'],dict(status='validated'),reuse_absolute=True)
    assert identity(p)==before
    assert cfg['drizzle']['final_pixfrac']==.8 and cfg['drizzle']['driz_sep_pixfrac']==1.
