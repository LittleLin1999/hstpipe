"""HAP's maintained matcher/fitter, with independent policy and diagnostics."""
import numpy as np
from astropy.table import Table
from tweakwcs.matchutils import XYXYMatch
from tweakwcs.linearfit import iter_linear_fit

# Mathematical minima. Scientific validation is assessed separately below.
MINIMUM={'shift':1,'rshift':2,'rscale':2,'general':3}
PARAMETERS={'shift':2,'rshift':3,'rscale':4,'general':6}


def match(reference,image,options,histogram=False,guess=None):
    a=Table(dict(TPx=reference[:,0],TPy=reference[:,1]));b=Table(dict(TPx=image[:,0],TPy=image[:,1]))
    kw=dict(searchrad=options['search_radius_pix'],tolerance=options['match_tolerance_pix'],separation=options['separation_pix'],use2dhist=histogram)
    if guess is not None:kw.update(xoffset=float(guess[0]),yoffset=float(guess[1]))
    i,j=XYXYMatch(**kw)(a,b,tp_pscale=1.,tp_units='pixels')
    return np.asarray(i),np.asarray(j)

def catalog_match(reference,image,refcat,imgcat,options,histogram=False,guess=None,matcher=None):
    """Common IDs select pairs; coordinate arrays contain HST measurements only."""
    if 'common_source_id' in refcat.colnames and 'common_source_id' in imgcat.colnames:
        a=np.asarray(refcat['common_source_id']).astype(str);b=np.asarray(imgcat['common_source_id']).astype(str)
        if len(np.unique(a))!=len(a) or len(np.unique(b))!=len(b):raise ValueError('Ambiguous duplicate HST source identities')
        lookup={k:i for i,k in enumerate(b)};i=[i for i,k in enumerate(a) if k in lookup]
        return np.asarray(i,int),np.asarray([lookup[a[k]] for k in i],int)
    return (matcher or match)(reference,image,options,histogram,guess)



def _parameter_errors(image,error,matrix,center,geometry,residual):
    """Weighted local information matrix; measurement errors never become zero."""
    q=(image-center)@matrix.T;n=len(q);jac=np.zeros((n,2,PARAMETERS[geometry]))
    jac[:,0,0]=1.;jac[:,1,1]=1.
    if geometry in ('rshift','rscale'):
        jac[:,0,2]=-q[:,1];jac[:,1,2]=q[:,0]
        if geometry=='rscale':jac[:,:,3]=q
    elif geometry=='general':
        x=image-center
        jac[:,0,2:4]=x;jac[:,1,4:6]=x
    a=(jac/error[:,None,None]).reshape(-1,PARAMETERS[geometry])
    if np.linalg.matrix_rank(a)<a.shape[1]:raise ValueError('Rank-deficient transformation')
    cov=np.linalg.inv(a.T@a);dof=2*n-a.shape[1]
    chi2=float(np.sum((residual/error[:,None])**2))
    cov*=max(1.,chi2/dof) if dof>0 else 1.
    correlation=cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
    return cov,correlation,chi2,dof



def _native_fit(reference,image,weights,geometry,nclip):
    """Native fit, with SciPy's proper-rotation constraint for reflection ties.

    A two-point covariance is rank one. Installed tweakwcs can select an
    improper reflection in that tie. Rotation.align_vectors with an infinite
    normal-vector weight constrains the maintained solver to the TAN plane
    and a proper rotation; it does not add an astronomical source constraint.
    """
    f=iter_linear_fit(reference,image,wxy=weights,fitgeom=geometry,nclip=nclip,sigma=(3.,'rmse'))
    if geometry not in ('rshift','rscale') or np.linalg.det(np.asarray(f['matrix'],float))>0:return f
    from scipy.spatial.transform import Rotation
    keep=np.asarray(f['fitmask']);a=reference[keep];b=image[keep];w=weights[keep]
    am=np.average(a,axis=0,weights=w);bm=np.average(b,axis=0,weights=w)
    av=np.vstack(([0.,0.,1.],np.column_stack((a-am,np.zeros(len(a))))))
    bv=np.vstack(([0.,0.,1.],np.column_stack((b-bm,np.zeros(len(b))))))
    rotation,_=Rotation.align_vectors(av,bv,weights=np.r_[np.inf,w]);matrix=rotation.as_matrix()[:2,:2]
    scale=1.
    if geometry=='rscale':
        scale=float(np.sum(w[:,None]*(a-am)*((b-bm)@matrix.T))/np.sum(w[:,None]*(b-bm)**2))
        if scale<=0:raise ValueError('Nonphysical similarity scale')
    matrix*=scale;center=np.asarray(f['center'],float);shift=am-((bm-center)@matrix.T+center)
    f.update(matrix=matrix,shift=shift,**{'<scale>':scale})
    f['proper_rotation_solver']='scipy.spatial.transform.Rotation.align_vectors (fixed TAN normal)'
    return f

def quality(reference,image,error,options,geometry='shift',identities_verified=False):
    """Native fit with rank/uncertainty checks, separate from independent QA."""
    reference=np.asarray(reference,float);image=np.asarray(image,float);errors=np.asarray(error,float)
    n=len(image);r=dict(status='insufficient_matches',matched=n,geometry=geometry,mathematical_minimum=MINIMUM[geometry])
    if n<MINIMUM[geometry]:return r
    if not np.isfinite(reference).all() or not np.isfinite(image).all():return dict(r,status='nonfinite')
    if not np.isfinite(errors).all() or (errors<=0).any():return dict(r,status='nonfinite_centroid_errors')
    if min(len(np.unique(image,axis=0)),len(np.unique(reference,axis=0)))<MINIMUM[geometry]:return dict(r,status='degenerate_distribution')
    if geometry=='general' and np.linalg.matrix_rank(image-image.mean(axis=0))<2:return dict(r,status='degenerate_distribution')
    weights=1/errors**2
    try:
        f=_native_fit(reference,image,weights,geometry,nclip=3)
        keep=np.asarray(f['fitmask']);matrix=np.asarray(f['matrix'],float);center=np.asarray(f['center'],float);shift=np.asarray(f['shift'],float)
        if keep.sum()<MINIMUM[geometry] or not np.isfinite(np.r_[matrix.ravel(),center,shift]).all():return dict(r,status='unreliable_fit')
        predicted=(image-center)@matrix.T+center+shift;res=reference[keep]-predicted[keep]
        cov,corr,chi2,dof=_parameter_errors(image[keep],errors[keep],matrix,center,geometry,res)
    except (ValueError,np.linalg.LinAlgError):return dict(r,status='degenerate_distribution')
    rms=lambda a:float(np.sqrt(np.mean(np.sum(np.asarray(a,float)**2,axis=1))))
    loo=[];shifts=[]
    if keep.sum()>MINIMUM[geometry]:
        for k in np.flatnonzero(keep):
            use=keep.copy();use[k]=False
            try:
                g=_native_fit(reference[use],image[use],weights[use],geometry,nclip=0)
                c=np.asarray(g['center'],float);m=np.asarray(g['matrix'],float);s=np.asarray(g['shift'],float)
                loo.append(reference[k]-((image[k]-c)@m.T+c+s))
                shifts.append((center-c)@m.T+c+s-center)
            except (ValueError,np.linalg.LinAlgError):loo=[];shifts=[];break
    delta=reference-image;offset=np.median(delta,axis=0);err=errors[keep]
    expected=float(np.sqrt(np.mean(err**2)));limit=min(options['max_scatter_pix'],max(options['scatter_floor_pix'],options['centroid_sigma']*expected))
    fitrms=rms(res);actual=rms(delta[keep]);w=weights[keep];effective=float(w.sum()**2/np.sum(w**2))
    formal=float(np.sqrt(np.linalg.eigvalsh(cov[:2,:2]).max()))
    uncertainty=max(formal,float(np.max(np.std(shifts,axis=0))) if shifts else formal)
    rot_error=float(np.rad2deg(np.sqrt(cov[2,2]))) if geometry in ('rshift','rscale') else None
    scale=float(f['<scale>']);scale_error=float(scale*np.sqrt(cov[3,3])) if geometry=='rscale' else None
    span=float(np.linalg.norm(np.ptp(image[keep],axis=0)))
    rot_constrained=rot_error is not None and options['centroid_sigma']*rot_error<=options.get('max_rotation_uncertainty_deg',options['max_rotation_deg'])
    r.update(robust_inliers=int(keep.sum()),fit_inlier_mask=keep.tolist(),robust_scatter_pix=float(1.4826*np.linalg.norm(np.median(np.abs(res-np.median(res,axis=0)),axis=0))),
        effective_matches=effective,max_weight_fraction=float(w.max()/w.sum()),initial_offset_pix=offset.tolist(),initial_rms_pix=rms(delta),postfit_scatter_pix=fitrms,
        weighted_rms_pix=float(np.sqrt(np.average(np.sum(res**2,axis=1),weights=w))),current_rms_pix=actual,centroid_precision_pix=expected,scatter_limit_pix=limit,
        loo_rms_pix=rms(loo) if loo else None,independent_validation='leave_one_source_out' if loo else 'not_available',transform_uncertainty_pix=uncertainty,
        parameter_covariance=cov.tolist(),parameter_correlations=corr.tolist(),parameter_order=['dx','dy']+(['theta_rad'] if geometry in ('rshift','rscale') else [])+(['log_scale'] if geometry=='rscale' else [])+(['m11','m12','m21','m22'] if geometry=='general' else []),
        dx_uncertainty_pix=float(np.sqrt(cov[0,0])),dy_uncertainty_pix=float(np.sqrt(cov[1,1])),rotation_uncertainty_deg=rot_error,scale_uncertainty=scale_error,
        rotation_measured=geometry in ('rshift','rscale'),rotation_constrained=rot_constrained,chi2=chi2,degrees_of_freedom=dof,
        rejected=int((~keep).sum()),outlier_fraction=float((~keep).mean()),baseline_pix=span,center=center.tolist(),shift=shift.tolist(),matrix=matrix.tolist(),
        rotation_deg=float(np.rad2deg(np.arctan2(matrix[1,0],matrix[0,0]))),rotation_convention='Counterclockwise mapping of measured image TAN coordinates into reference TAN coordinates',scale=scale,residuals_pix=(reference-predicted).tolist(),
        validation_scope='local_shift_only' if keep.sum()==1 else ('limited_redundancy' if not loo else 'source_constrained'),
        residual_trend_pix=float(np.max(np.abs(np.nan_to_num(np.corrcoef(np.column_stack([image[keep],res]).T)[:2,2:])))) if len(res)>2 and np.std(res)>1e-10 else None)
    if f.get('proper_rotation_solver'):r['proper_rotation_solver']=f['proper_rotation_solver']
    good=r['outlier_fraction']<=options['max_outlier_fraction'] and fitrms<=limit and uncertainty<=options['max_transform_error_pix']
    if loo:good&=r['loo_rms_pix']<=limit*1.5
    if geometry in ('rshift','rscale'):good&=rot_constrained
    if keep.sum()>10:good&=(r['residual_trend_pix'] or 0)<=options.get('max_residual_trend',.6)
    good&=np.linalg.norm(shift)<=options['max_shift_pix'] and abs(r['rotation_deg'])<=options['max_rotation_deg'] and abs(scale-1)<=options['max_scale_change']
    r['status']='validated' if good and actual<=options['max_current_rms_pix'] and np.linalg.norm(offset)<=options['max_offset_pix'] else ('refinement_supported' if good else 'unreliable_fit')
    if good and not loo and not identities_verified:r['status']='candidate_fit_only'
    # These are source-fit assessments. A one-source result is always local;
    # exact/minimal fits are never labelled independent full-field validation.
    return r


def select_model(reference,image,error,options):
    """Prefer shift unless measured rotation/scale and predictive QA justify it."""
    requested=options.get('fit_geometry','auto');models=['shift','rshift','rscale'] if requested=='auto' else ['shift']+([requested] if requested!='shift' else [])
    candidates={g:quality(reference,image,error,options,g,identities_verified=True) for g in models};chosen=candidates['shift'];reasons={}
    sigma=options.get('model_significance_sigma',options['centroid_sigma'])
    for geometry in models[1:]:
        fit=candidates[geometry]
        if fit['status'] not in ('validated','refinement_supported'):
            reasons[geometry]='Rank, measurement precision, physical bounds or residual validation failed';continue
        rotation=abs(fit['rotation_deg'])/max(fit['rotation_uncertainty_deg'],1e-15)
        scale=abs(fit['scale']-1)/max(fit['scale_uncertainty'] or np.inf,1e-15) if geometry=='rscale' else 0.
        significant=rotation>=sigma if geometry=='rshift' else scale>=sigma
        gain=chosen.get('chi2',np.inf)-fit['chi2'];added=PARAMETERS[geometry]-PARAMETERS[chosen['geometry']]
        # Compare the same objects: clipping cannot improve model evidence by
        # silently changing the training sample.
        same_sample=fit.get('fit_inlier_mask')==chosen.get('fit_inlier_mask')
        improved=chosen.get('postfit_scatter_pix',np.inf)-fit['postfit_scatter_pix']>=options['min_improvement_pix']
        if significant and same_sample and gain>=sigma**2*added and improved:
            chosen=fit;reasons[geometry]='Supported parameter significance, weighted residual improvement and prediction checks'
        else:reasons[geometry]='Extra parameter not significantly supported, sample changed, or improvement below configured requirement'
    return dict(selected=chosen,candidates=candidates,reasons=reasons)


def joint_models(xy,catalogs,options):
    """Weighted common-frame adjustment of verified HST measurements only.

    Native tweakwcs fits seed and select the per-exposure rigid/similarity
    models. SciPy's maintained least_squares solves their shared source
    coordinates simultaneously, with exposure zero fixing the frame gauge.
    Missing measurements add no constraints and receive no invented positions.
    """
    from scipy.optimize import least_squares
    ids=[np.asarray(c['common_source_id']).astype(str) for c in catalogs]
    if any(len(set(a))!=len(a) for a in ids):return dict(status='ambiguous_source_ids')
    counts={sid:sum(sid in a for a in ids) for sid in set(np.concatenate(ids))}
    common=sorted(sid for sid,n in counts.items() if n>=2)
    selected=[np.isin(a,common) for a in ids];xy=[np.asarray(a,float)[s] for a,s in zip(xy,selected)]
    errors=[np.asarray(c['centroid_error_pix'],float)[s] for c,s in zip(catalogs,selected)];ids=[a[s] for a,s in zip(ids,selected)]
    if not common or any(not len(a) for a in xy):return dict(status='disconnected',reason='An exposure has no common measured object')
    if any(not np.isfinite(a).all() for a in xy+errors) or any((e<=0).any() for e in errors):return dict(status='nonfinite')
    if any(len(xy[i])<1 for i in range(len(xy))):return dict(status='disconnected')
    lookup={s:k for k,s in enumerate(common)};source_indices=[np.array([lookup[s] for s in a]) for a in ids]
    pivots=[np.average(a,axis=0,weights=1/e**2) for a,e in zip(xy,errors)]
    consensus_errors={sid:errors[0][k] for k,sid in enumerate(ids[0])}
    consensus={sid:xy[0][k] for k,sid in enumerate(ids[0])};seed={0:dict(center=pivots[0].tolist(),shift=[0.,0.],matrix=np.eye(2).tolist(),geometry='shift')};comparisons={};pending=set(range(1,len(xy)))
    while pending:
        choices=[]
        for i in pending:
            use=np.array([s in consensus for s in ids[i]])
            if not use.any():continue
            ref=np.array([consensus[s] for s in ids[i][use]])
            pair_errors=np.hypot([consensus_errors[s] for s in ids[i][use]],errors[i][use])
            model=select_model(ref,xy[i][use],pair_errors,options)
            if 'matrix' not in model['selected']:continue
            choices.append((int(use.sum()),-float(np.mean(errors[i][use])),-i,i,model))
        if not choices:return dict(status='disconnected',reason='No full-rank HST overlap path to anchored exposure')
        *_,i,model=max(choices,key=lambda v:v[:3]);f=model['selected'];seed[i]=f;comparisons[i]=model
        mapped=(xy[i]-f['center'])@np.asarray(f['matrix']).T+f['center']+np.asarray(f['shift'])
        for j,(sid,pos) in enumerate(zip(ids[i],mapped)):
            if sid not in consensus:
                consensus[sid]=pos;consensus_errors[sid]=np.hypot(errors[i][j],f['transform_uncertainty_pix'])
        pending.remove(i)
    index=[(i,j) for i in range(len(xy)) for j in range(len(xy[i]))]
    sigma=np.asarray([errors[i][j] for i,j in index]);observed=np.asarray([xy[i][j] for i,j in index])
    source_row=np.asarray([source_indices[i][j] for i,j in index]);exposure_row=np.asarray([i for i,j in index])
    def solve(models,clip=False):
        starts={};parameters=[]
        for i in range(1,len(xy)):
            starts[i]=len(parameters);f=seed[i];center=np.asarray(f['center']);matrix=np.asarray(f['matrix']);p=pivots[i]
            shift=(p-center)@matrix.T+center+np.asarray(f['shift'])-p
            parameters.extend(shift)
            if models[i] in ('rshift','rscale'):parameters.append(float(np.arctan2(matrix[1,0],matrix[0,0])))
            if models[i]=='rscale':parameters.append(float(np.log(np.sqrt(np.linalg.det(matrix)))))
        source_start=len(parameters);parameters.extend(np.asarray([consensus[s] for s in common]).ravel())
        def unpack(v,i):
            if i==0:return np.zeros(2),np.eye(2)
            k=starts[i];shift=v[k:k+2];theta=v[k+2] if models[i]!='shift' else 0.
            scale=np.exp(v[k+3]) if models[i]=='rscale' else 1.
            return shift,scale*np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        def residual(v):
            latent=v[source_start:].reshape(-1,2);mapped=observed.copy()
            for i in range(len(xy)):
                use=exposure_row==i;shift,matrix=unpack(v,i)
                mapped[use]=(observed[use]-pivots[i])@matrix.T+pivots[i]+shift
            return ((mapped-latent[source_row])/sigma[:,None]).ravel()
        robust=least_squares(residual,np.asarray(parameters),loss='soft_l1',f_scale=options['centroid_sigma'],jac='3-point',x_scale='jac')
        keep=np.linalg.norm(residual(robust.x).reshape(-1,2),axis=1)<=options['centroid_sigma'] if clip else np.ones(len(index),bool)
        # A minimally constrained observation cannot be clipped into a valid fit.
        for i in range(len(xy)):
            if np.count_nonzero(keep&(exposure_row==i))<MINIMUM[models[i]]:return dict(status='underconstrained_after_rejection')
        mask=np.repeat(keep,2)
        final=least_squares(lambda v:residual(v)[mask],robust.x,jac='3-point',x_scale='jac')
        jac=final.jac;rank=np.linalg.matrix_rank(jac)
        if not final.success or rank<len(final.x) or not np.isfinite(final.x).all():return dict(status='degenerate_joint_fit')
        dof=len(final.fun)-len(final.x);chi2=float(final.fun@final.fun)
        covariance=np.linalg.inv(jac.T@jac)*max(1.,chi2/dof) if dof>0 else np.linalg.inv(jac.T@jac)
        delta=residual(final.x).reshape(-1,2)*sigma[:,None];fits=[]
        for i in range(len(xy)):
            use=exposure_row==i;inside=use&keep;n=int(use.sum());shift,matrix=unpack(final.x,i)
            geom=models[i];k=starts.get(i,0);size=PARAMETERS[geom]
            cov=covariance[k:k+size,k:k+size] if i else np.zeros((2,2))
            uncertainty=float(np.sqrt(np.linalg.eigvalsh(cov[:2,:2]).max()))
            rotation_error=float(np.rad2deg(np.sqrt(cov[2,2]))) if i and geom!='shift' else None
            r=dict(geometry=geom,center=pivots[i].tolist(),shift=shift.tolist(),matrix=matrix.tolist(),rotation_deg=float(np.rad2deg(np.arctan2(matrix[1,0],matrix[0,0]))) if geom!='shift' else None,
                rotation_measured=bool(i and geom!='shift'),rotation_uncertainty_deg=rotation_error,rotation_constrained=bool(rotation_error is not None and options['centroid_sigma']*rotation_error<=options.get('max_rotation_uncertainty_deg',options['max_rotation_deg'])),
                scale=float(np.sqrt(np.linalg.det(matrix))),scale_uncertainty=float(np.sqrt(np.linalg.det(matrix))*np.sqrt(cov[3,3])) if i and geom=='rscale' else None,
                dx_uncertainty_pix=float(np.sqrt(cov[0,0])),dy_uncertainty_pix=float(np.sqrt(cov[1,1])),transform_uncertainty_pix=uncertainty,parameter_covariance=cov.tolist(),
                parameter_correlations=np.divide(cov,np.sqrt(np.outer(np.diag(cov),np.diag(cov))),out=np.zeros_like(cov),where=np.outer(np.diag(cov),np.diag(cov))>0).tolist(),
                matched=n,robust_inliers=int(inside.sum()),rejected=int((use&~keep).sum()),outlier_fraction=float(np.count_nonzero(use&~keep)/n),
                postfit_scatter_pix=float(np.sqrt(np.mean(np.sum(delta[inside]**2,axis=1)))),robust_scatter_pix=float(1.4826*np.linalg.norm(np.median(np.abs(delta[inside]-np.median(delta[inside],axis=0)),axis=0))),
                initial_rms_pix=float(np.sqrt(np.mean(np.sum((observed[use]-np.array([consensus[common[s]] for s in source_row[use]]))**2,axis=1)))),
                centroid_precision_pix=float(np.sqrt(np.mean(sigma[inside]**2))),baseline_pix=float(np.linalg.norm(np.ptp(xy[i],axis=0))),source_ids=ids[i].tolist(),
                validation_scope='anchor_frame' if i==0 else ('local_shift_only' if n==1 else ('limited_redundancy' if n==MINIMUM[geom] else 'source_constrained')),
                independent_validation='No independent catalog or held-out objects in joint fit; model stability and no-CR image QA assessed separately',
                rotation_fixed_reason='Anchored frame' if i==0 else ('Single source: retain calibrated rotation' if n==1 else ('Rotation not statistically justified' if geom=='shift' else None)))
            fits.append(r)
        return dict(status='fit_supported',fits=fits,chi2=chi2,degrees_of_freedom=dof,measurement_inlier_mask=keep.tolist(),measurement_residuals_pix=delta.tolist(),
                    residual_rms_pix=float(np.sqrt(np.mean(np.sum(delta[keep]**2,axis=1)))),residual_definition='Individual HST measurements versus shared fitted source positions; not independent validation',source_positions=final.x[source_start:].reshape(-1,2).tolist())
    baseline=solve(['shift']*len(xy));models=['shift']+[seed[i]['geometry'] for i in range(1,len(xy))]
    auto=solve(models) if any(m!='shift' for m in models) else None;chosen=baseline;reason='Shift is the simplest supported model; no significant rotation/scale improvement'
    if auto and auto['status']=='fit_supported':
        extra=sum(PARAMETERS[m]-2 for m in models);sig=options.get('model_significance_sigma',options['centroid_sigma'])
        rotation_ok=all(not f['rotation_measured'] or (f['rotation_constrained'] and (abs(f['rotation_deg'])>=sig*f['rotation_uncertainty_deg'] if f['geometry']=='rshift' else abs(f['scale']-1)>=sig*f['scale_uncertainty'])) for f in auto['fits'])
        if baseline['status']=='fit_supported' and rotation_ok and baseline['chi2']-auto['chi2']>=sig**2*extra and baseline['residual_rms_pix']-auto['residual_rms_pix']>=options['min_improvement_pix'] and baseline['measurement_inlier_mask']==auto['measurement_inlier_mask']:
            chosen=auto;reason='Measured parameter significance and joint weighted residual improvement support extra degrees of freedom'
    chosen_models=models if chosen is auto else ['shift']*len(xy)
    chosen=solve(chosen_models,clip=True) if chosen.get('status')=='fit_supported' else dict(chosen)
    chosen.update(reference_index=0,distinct_sources=len(common),source_ids=common,measurement_count=len(index),external_coordinates_used_for_fit=False,
        model_selection_reason=reason,shift_only=baseline if chosen is not baseline else {k:v for k,v in baseline.items() if k not in ('shift_only','rotation_enabled','model_comparisons')},rotation_enabled=auto,
        model_comparisons={str(i):v for i,v in comparisons.items()})
    return chosen
