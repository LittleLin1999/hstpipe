import json
from PIL import Image
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.target_runner import save
from hstpipe.uvis.reporting import generate_report


def record(accepted):
    rows=[dict(source_id=f'star{i}',status='measured',fwhm_pix=12.,elongation=4.) for i in range(3)]
    after=[dict(r,fwhm_pix=3.,elongation=1.1) for r in rows]
    joint=dict(status='validated' if accepted else 'unvalidated',may_combine=accepted,corrections_applied=accepted,
               image_quality=dict(sources_by_phase=dict(before=rows,after=after)))
    return dict(filter='F621M',status='relative_astrometry_unvalidated',
                scientific_assessment=dict(alignment=dict(status='validated' if accepted else 'unvalidated',joint_models=joint)))


def test_report_contains_three_identified_offline_before_after_pairs(tmp_path):
    p=ReductionPaths(tmp_path,'TEST','uvis','F621M')
    plot=p.file('diagnostics','alignment','joint_relative_coadd_cutouts','png');plot.parent.mkdir(parents=True)
    Image.new('RGB',(1800,300),'grey').save(plot)
    original=dict(target='TEST',filters={'F621M':record(True)})
    save(p.directory/'manifest.json',original)
    output=generate_report(p.directory);text=output.read_text()
    assert 'Before vs after alignment' in text and 'Left: before. Right: after.' in text
    assert '3 source pairs' in text and 'Accepted relative correction' in text
    assert text.count('class="alignment-pair"')==3
    assert all(f'F621M star{i} before and after relative alignment' in text for i in range(3))
    assert 'no-CR diagnostic coadds, not the final DRC' in text
    assert 'data:image/jpeg;base64,' in text and output.stat().st_size<1_000_000
    state=json.loads((p.directory/'manifest.json').read_text())
    assert state['filters']==original['filters'] and state['report']['alignment_source_pairs']==3


def test_rejected_candidate_is_not_labelled_applied_and_old_unpaired_plots_skipped(tmp_path):
    p=ReductionPaths(tmp_path,'TEST','uvis','F621M')
    plot=p.file('diagnostics','alignment','joint_relative_coadd_cutouts','png');plot.parent.mkdir(parents=True)
    Image.new('RGB',(1800,300),'grey').save(plot)
    state=dict(target='TEST',filters={'F621M':record(False)})
    save(p.directory/'manifest.json',state)
    text=generate_report(p.directory).read_text()
    assert 'Candidate correction — not applied' in text and 'Accepted relative correction' not in text
    state['filters']['F621M']['scientific_assessment']['alignment']['joint_models']['image_quality']={}
    save(p.directory/'manifest.json',state)
    assert 'class="alignment-pair"' not in generate_report(p.directory).read_text()
