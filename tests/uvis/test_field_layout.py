"""Explicit field membership, immutable input identity, and layout compatibility."""
import json
from pathlib import Path
import pytest
from astropy.io import fits
from hstpipe.common.config import ReductionConfig
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.common.target_runner import run_uvis,save
from hstpipe.uvis.grouping import group_exposures
from hstpipe.uvis.inventory import discover
from hstpipe.uvis.reporting import generate_report


def test_explicit_field_groups_different_targets_and_preserves_headers(flc,tmp_path):
    raw=tmp_path/'raw';raw.mkdir()
    for i,name in enumerate(('VISIT-A','VISIT-B',None)):
        p=flc(str(i))
        with fits.open(p,mode='update') as h:
            if name is None:del h[0].header['TARGNAME']
            else:h[0].header['TARGNAME']=name
        p.rename(raw/p.name)
    before={str(p):identity(p) for p in raw.glob('*.fits')}
    rows=discover(raw)
    with pytest.raises(ValueError,match='--field'):group_exposures(rows)
    grouped=group_exposures(rows,field='GOODS-S')
    assert list(grouped)==[('GOODS-S','F275W',None)]
    assert len(next(iter(grouped.values())))==3
    cfg=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',field='GOODS-S')
    plans=run_uvis(cfg,None,dry_run=True)
    assert len(plans)==1 and plans[0]['output']==str(tmp_path/'redux/uvis/GOODS-S')
    assert not Path(cfg.output).exists()
    seen=[]
    def engine(config,paths,exposures):
        seen.extend(e.metadata['TARGNAME'] for e in exposures)
        return []
    [directory]=run_uvis(cfg,engine,_retry_anchors=False)
    state=json.loads((directory/'manifest.json').read_text())
    assert state['field']=='GOODS-S' and state['grouping']['mode']=='explicit_field'
    assert set(seen)=={'VISIT-A','VISIT-B',None}
    assert {str(p):identity(p) for p in raw.glob('*.fits')}==before
    assert (directory/'index.html').is_file()


def test_explicit_field_still_separates_disconnected_footprints(flc,tmp_path):
    flc('a',ra=10);flc('b',ra=12)
    rows=discover(tmp_path)
    groups=group_exposures(rows,field='DEEP')
    assert len(groups)==2 and all(key[2].startswith('field-') for key in groups)
    assert set(groups)==set(group_exposures(list(reversed(rows)),field='DEEP'))


@pytest.mark.parametrize('legacy',[False,True])
def test_report_and_cleanup_accept_both_layouts(tmp_path,legacy):
    from hstpipe.common.uvis_cleanup import cleanup_uvis
    from hstpipe.cli import main
    directory=tmp_path/('T/uvis' if legacy else 'uvis/T')
    save(directory/'manifest.json',dict(schema_version=2,instrument='uvis',target='T',filters={'F275W':dict(filter='F275W',status='incomplete')}))
    paths=ReductionPaths.from_directory(directory,'F275W',target='T')
    assert paths.directory==directory
    assert ReductionPaths(tmp_path,'T','uvis','F275W').directory==directory
    main(['uvis','report','--directory',str(directory)])
    assert (directory/'index.html').exists()
    assert cleanup_uvis(tmp_path)['files']==[]
    assert cleanup_uvis(directory)['files']==[]


def test_independent_legacy_and_current_trees_are_not_overwritten(tmp_path):
    (tmp_path/'T/uvis').mkdir(parents=True)
    (tmp_path/'uvis/T').mkdir(parents=True)
    with pytest.raises(ValueError,match='Ambiguous'):
        ReductionPaths(tmp_path,'T','uvis','F275W').directory


def test_config_rejects_unsafe_or_conflicting_field(tmp_path):
    for field in ('../elsewhere','bad/name','bad name'):
        with pytest.raises(ValueError):ReductionConfig(str(tmp_path/'raw'),str(tmp_path/'redux'),'uvis',field=field)
    with pytest.raises(ValueError,match='must match'):
        ReductionConfig(str(tmp_path/'raw'),str(tmp_path/'redux'),'uvis',target='A',field='B')


def test_publication_view_cannot_be_overwritten_by_reduction(tmp_path):
    current=tmp_path/'uvis/T'
    save(current/'publication.json',dict(kind='published_results'))
    with pytest.raises(ValueError,match='published-results view'):
        ReductionPaths(tmp_path,'T','uvis','F275W').directory
