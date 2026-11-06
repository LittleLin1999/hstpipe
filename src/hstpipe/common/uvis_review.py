"""Human-facing UVIS review adapter; scientific approval stays in source_qa."""
import json
from pathlib import Path
from .paths import ReductionPaths
from .provenance import identity


def review_uvis(directory, filter_key, *, approve=False, reviewer=None, reason=None):
    directory = Path(directory).expanduser().resolve()
    state = json.loads((directory / 'manifest.json').read_text())
    record = state['filters'][filter_key]
    paths = ReductionPaths.from_directory(directory, record['filter'], record.get('setup'), target=state['target'])
    if paths.directory != directory or paths.group_key != filter_key:
        raise ValueError('Use the actual target reduction directory and filter/setup key')
    assessment = record.get('scientific_assessment', {})
    image = Path(assessment.get('unvalidated_final', ''))
    if not image.is_absolute() or not image.resolve().is_relative_to(directory / 'intermediate'):
        raise ValueError('No local review image recorded; use the actual reduction, not a merged report view')
    sources = [s for s in assessment.get('source_validation', {}).get('sources', [])
               if s.get('status') == 'review_required']
    if not sources:
        raise ValueError('No flagged source assessment to review')
    previous = assessment.get('source_review', {}).get('input_image')
    candidates = ([previous] if previous else []) + record.get('artifacts', []) + record.get('intermediate_files', [])
    recorded = [v for v in candidates if v.get('path') == str(image)]
    if not recorded or any(v != recorded[0] for v in recorded):
        raise ValueError('Missing or inconsistent recorded image identity; recover its QA provenance before approval')
    if identity(image) != recorded[0]:
        raise ValueError('Review FITS changed since its recorded assessment; review and validate the new image first')
    report=paths.report()
    legacy=paths.area('products')/f'{paths.target}_uvis_report.html'
    if not report.is_file() and legacy.is_file():report=legacy
    result = dict(target=state['target'], filter=filter_key, status=record['status'],
                  image=str(image), report=str(report),
                  sources=[{k:s.get(k) for k in ('source_id','failures','output_x','output_y','apertures')} for s in sources],
                  approved=False)
    if approve:
        if not reviewer or not reviewer.strip() or not reason or not reason.strip():
            raise ValueError('Approval requires --reviewer and a scientific --reason')
        from ..uvis import approve_source_review
        product = approve_source_review(directory, filter_key, reviewer=reviewer,
            reason=reason, source_ids=[s['source_id'] for s in sources],
            expected_sha256=recorded[0]['sha256'])
        result.update(approved=True, product=str(product), report=str(paths.report()))
    return result
