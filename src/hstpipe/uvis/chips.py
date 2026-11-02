"""Chip provenance and diagnostics; detector geometry always comes from FITS."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack

def geometry(image):
    with fits.open(image) as hd:
        return [dict(extver=int(h.header.get('EXTVER',1)),ccdchip=h.header.get('CCDCHIP'),shape=list(h.shape),
                     err_extver=int(h.header.get('EXTVER',1)),dq_extver=int(h.header.get('EXTVER',1)),
                     wcsname=h.header.get('WCSNAME'),processed=True,excluded=False)
                for h in hd if h.name=='SCI']

def source_products(image,catalog,methods,paths):
    """Materialize per-chip QC tables even when detector catalogs are reused."""
    result=[];audits=[]
    for method in methods:
        if method.get('quality_file') and Path(method['quality_file']).exists():
            a=Table.read(method['quality_file'],format='ascii.ecsv');a.meta={};audits.append(a)
    allrows=vstack(audits,join_type='outer') if audits else catalog.copy()
    root=fits.getheader(image).get('ROOTNAME',Path(image).stem.removesuffix('_flc'))
    for chip in geometry(image):
        v=chip['extver'];rows=allrows[allrows['extver']==v].copy();rows['ccdchip']=np.full(len(rows),chip['ccdchip'] if chip['ccdchip'] is not None else -1)
        rows.meta={'extver':v,'ccdchip':chip['ccdchip'],'detector_rows_may_overlap':True,
                   'provenance':'Native HST positions. Rejected detector candidates retained; valid unique sources are counted after cross-method deduplication.'}
        output=paths.exposure_file('intermediate','alignment',root,f'sci{v}_sources','ecsv')
        rows.write(output,format='ascii.ecsv',overwrite=True)
        counts=[]
        for method,audit in zip(methods,audits):
            use=audit['extver']==v;valid=np.asarray(audit['selected_for_astrometry'],bool) if 'selected_for_astrometry' in audit.colnames else np.asarray(audit['centroid_valid'],bool)
            counts.append(dict(method=method['method'],kernel_fwhm_pix=method['kernel_fwhm_pix'],detected=int(use.sum()),valid=int((use&valid).sum()),rejected=int((use&~valid).sum())))
        good=catalog[catalog['extver']==v]
        recovered=np.asarray(good['recovered_from_inherited_cr'],bool) if 'recovered_from_inherited_cr' in good.colnames else np.zeros(len(good),bool)
        chip.update(recovered_valid_bright_sources=int(np.count_nonzero(recovered&(np.asarray(good['snr'])>=20))),bright_snr_definition=20,methods=counts,detected_method_rows=len(rows),valid_unique=int(np.count_nonzero(catalog['extver']==v)),catalog=str(output))
        result.append(chip)
    return result

def match_provenance(reference,source,ri,ii,delta,assessment,options):
    """Track paired native SCI IDs and flag contradictory chip residual evidence."""
    ri=np.asarray(ri,int);ii=np.asarray(ii,int);delta=np.asarray(delta,float)
    keep=np.asarray(assessment.get('fit_inlier_mask',np.ones(len(ri),bool)),bool)
    shift=np.asarray(assessment.get('shift',[0.,0.]),float)
    details=[]
    for side,cat,index in (('reference',reference,ri),('exposure',source,ii)):
        for version in np.unique(cat['extver']):
            use=(np.asarray(cat['extver'])[index]==version)&keep
            d=delta[use]
            fitted=np.asarray(assessment.get('residuals_pix',delta-shift),float)
            res=fitted[use]
            ccd=int(cat['ccdchip'][index[np.flatnonzero(use)[0]]]) if use.any() and 'ccdchip' in cat.colnames else None
            details.append(dict(side=side,extver=int(version),ccdchip=ccd,candidates=int(np.count_nonzero(np.asarray(cat['extver'])[index]==version)),
                robust_inliers=int(use.sum()),offset_pix=np.median(d,axis=0).tolist() if len(d) else None,
                postfit_rms_pix=float(np.sqrt(np.mean(np.sum(res**2,axis=1)))) if len(d) else None,
                independently_constrained=int(use.sum())>=2))
    assessment['chip_matches']=details
    contradictory=[r for r in details if r['independently_constrained'] and r['postfit_rms_pix']>options['max_scatter_pix']]
    if contradictory and assessment['status'] in ('validated','refinement_supported'):
        assessment['status']='inconsistent_chip_residuals';assessment['chip_warning']='At least one chip contradicts the shared exposure-level fit; no independent chip shifts applied'
    return details

def report_summary(inputs,rows,alignment):
    catalogs={r['exposure']:r['path'] for r in alignment.get('detector_catalogs',[])}
    records=[]
    for image,row in zip(inputs,rows):
        record=dict(exposure=Path(image).name,rootname=row.metadata['ROOTNAME'],chips=geometry(image),all_sci_included=True)
        path=catalogs.get(Path(image).name)
        cat=Table.read(path,format='ascii.ecsv') if path and Path(path).exists() else None
        for chip in record['chips']:
            chip['relative_catalog_sources']=int(np.count_nonzero(cat['extver']==chip['extver'])) if cat is not None else None
            chip['distortion_model']='stwcs.wcsutil.HSTWCS; no detector resampling'
        records.append(record)
    return records

def coverage_report(product,inputs):
    """Sample expected chip footprints against native drizzle coverage, not ERR."""
    from astropy.wcs import WCS
    from .astrometry.wcs import detector_wcs
    from .fits import image_data
    records=[]
    with fits.open(product) as out:
        frame=WCS(out['SCI'].header,out);weight=out['WHT'].data;ny,nx=weight.shape
        for image in inputs:
            with fits.open(image) as hd:
                for h in hd:
                    if h.name!='SCI':continue
                    v=h.header.get('EXTVER',1);yy,xx=np.meshgrid(np.linspace(.1*h.shape[0],.9*h.shape[0],5),np.linspace(.1*h.shape[1],.9*h.shape[1],5))
                    xy=np.column_stack([xx.ravel(),yy.ravel()]);dq=image_data(hd['DQ',v]);ok=dq[xy[:,1].astype(int),xy[:,0].astype(int)]==0
                    xy=xy[ok];world=detector_wcs(hd,v).all_pix2world(xy,0);positions=frame.all_world2pix(world,0)
                    finite=np.isfinite(positions).all(axis=1);positions=np.rint(positions[finite]).astype(int)
                    inside=(positions[:,0]>=0)&(positions[:,0]<nx)&(positions[:,1]>=0)&(positions[:,1]<ny);p=positions[inside]
                    covered=int(np.count_nonzero(weight[p[:,1],p[:,0]]>0))
                    records.append(dict(exposure=Path(image).name,extver=int(v),ccdchip=h.header.get('CCDCHIP'),input_included=True,
                        usable_samples=int(len(xy)),covered_samples=covered,status='coverage_present' if covered else 'coverage_not_verified',
                        caveat='WHT sampling verifies expected footprint coverage, not unique per-chip flux attribution'))
    return records
