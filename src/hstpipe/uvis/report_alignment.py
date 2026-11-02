"""Compact paired-source previews from saved no-CR alignment diagnostics."""
import base64
import html
import io
from pathlib import Path
from ..common.paths import ReductionPaths


def previews(paths,filters,budget):
    """Reflow saved paired panels; never rerun astrometry or read full FITS data.

    Only the current paired-panel format, identified by per-source QA records,
    is eligible. Older unpaired plots must not be labelled before/after pairs.
    """
    from PIL import Image
    sections=[];used=0;count=0
    esc=lambda x:html.escape(str(x),quote=True)
    for key,record in sorted(filters.items()):
        a=record.get('scientific_assessment',{})
        rel=a.get('relative_astrometry',a.get('alignment',{}))
        joint=rel.get('joint_models',{})
        qa=joint.get('image_quality',{})
        phases=qa.get('sources_by_phase',{})
        after={r['source_id']:r for r in phases.get('after',[]) if r.get('status')=='measured'}
        sources=[r for r in phases.get('before',[]) if r.get('status')=='measured' and r['source_id'] in after][:3]
        group=ReductionPaths(paths.root,paths.target,'uvis',record.get('filter',key),record.get('setup'))
        imagefile=group.file('diagnostics','alignment','joint_relative_coadd_cutouts','png')
        if not sources or not imagefile.is_file():continue
        accepted=bool(joint.get('corrections_applied') and joint.get('may_combine'))
        state='Accepted relative correction' if accepted else 'Candidate correction — not applied'
        if joint.get('method')=='joint_existing_wcs':state='Original WCS retained; proposed correction shown for comparison'
        figures=[]
        try:
            with Image.open(imagefile) as strip:
                strip=strip.convert('RGB')
                for index,source in enumerate(sources):
                    # image_quality writes equal-width before/after panel pairs.
                    # Crop layout only: scientific intensities are not altered.
                    pair=strip.crop((round(index*strip.width/len(sources)),0,
                                     round((index+1)*strip.width/len(sources)),strip.height))
                    pair.thumbnail((760,440))
                    stream=io.BytesIO();pair.save(stream,format='JPEG',quality=80,optimize=True)
                    encoded=base64.b64encode(stream.getvalue()).decode('ascii')
                    sid=source['source_id'];new=after[sid]
                    caption=f"{key} · {sid}. Moment width: {source['fwhm_pix']:.2f} → {new['fwhm_pix']:.2f} pixels. Elongation: {source['elongation']:.2f} → {new['elongation']:.2f}."
                    fragment='<figure class="alignment-pair"><img alt="'+esc(key+' '+sid+' before and after relative alignment')+'" src="data:image/jpeg;base64,'+encoded+'"><figcaption>'+esc(caption)+'</figcaption></figure>'
                    size=len(fragment.encode('utf-8'))
                    if used+size+2048>budget:continue
                    figures.append(fragment);used+=size;count+=1
        except (OSError,ValueError):continue
        if not figures:continue
        note='<p><strong>Left: before. Right: after.</strong> Same source and angular scale; the cutouts include the original displaced images. Both are no-CR diagnostic coadds, not the final DRC. Display stretches may differ; assess positions and shape rather than brightness.</p>'
        warning=a.get('astrometry',{}).get('final_validation',{}).get('individual_absolute_status')
        if warning=='review_required':note+='<p class="attention">Final combined-image checks passed, but an individual-exposure external-catalog check still needs review. This is separate from relative alignment.</p>'
        section='<details class="alignment-comparison"><summary>'+esc(key)+' — '+str(len(figures))+' source pairs · '+esc(state)+'</summary>'+note+''.join(figures)+'</details>'
        used+=len(section.encode('utf-8'))-sum(len(f.encode('utf-8')) for f in figures)
        sections.append(section)
    return sections,count
