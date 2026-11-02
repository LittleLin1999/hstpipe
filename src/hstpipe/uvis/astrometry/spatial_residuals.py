"""Small first-order residual diagnostic; never a distortion/WCS correction.

Weighted planes are evaluated only at matched positions (the extrema of their
convex hull), not at detector corners. Known centroid errors and leave-one-out
sensitivity bound the amplitude uncertainty. Correlation alone cannot reject.
"""
import numpy as np


def spatial_residuals(position,residual,error,options):
    p=np.asarray(position,float);r=np.asarray(residual,float);e=np.asarray(error,float);n=len(p)
    out=dict(status='underconstrained',matched=n,correlations=[],max_abs_correlation=None,
        significant=False,exceeds_accuracy=False,amplitude_pix=None,amplitude_mas=None,
        amplitude_uncertainty_pix=None,scope='Matched-source convex hull only; no detector extrapolation')
    if n<2:return out
    for i in range(2):
        for j in range(2):
            value=float(np.corrcoef(p[:,i],r[:,j])[0,1]) if np.std(p[:,i])>0 and np.std(r[:,j])>1e-12 else None
            out['correlations'].append(dict(position='xy'[i],residual='xy'[j],matched=n,pearson=value))
    finite=[abs(c['pearson']) for c in out['correlations'] if c['pearson'] is not None]
    out['max_abs_correlation']=max(finite) if finite else None
    center=np.average(p,axis=0,weights=1/e**2);extent=np.ptp(p,axis=0)
    out.update(center_pix=center.tolist(),extent_pix=extent.tolist(),baseline_pix=float(np.linalg.norm(extent)))
    if n<=3 or np.any(extent==0):return out
    # Normalize the design to avoid numerical dependence on detector dimensions.
    q=(p-center)/extent;design=np.column_stack([np.ones(n),q]);a=design/e[:,None]
    if np.linalg.matrix_rank(a)<3:return out
    def fit(use):
        aw=a[use];coef=np.linalg.lstsq(aw,r[use]/e[use,None],rcond=None)[0]
        trend=q@coef[1:];amplitudes=np.linalg.norm(trend,axis=1)
        return coef,trend,amplitudes
    use=np.ones(n,bool);coef,trend,amplitudes=fit(use);peak=int(np.argmax(amplitudes));amplitude=float(amplitudes[peak])
    noise=r-design@coef;chi=np.sum((noise/e[:,None])**2,axis=0)
    cov=np.linalg.inv(a.T@a)
    # Offset is removed from the reported variation; propagate slope covariance.
    slope_variance=float(q[peak]@cov[1:,1:]@q[peak])
    formal=float(np.sqrt(max(0,slope_variance)*np.maximum(1,chi/(n-3)).sum()))
    loo=[];influence=[]
    for i in range(n):
        subset=use.copy();subset[i]=False
        if np.linalg.matrix_rank(a[subset])<3:continue
        c,t,amp=fit(subset);loo.append(float(amp.max()));influence.append(float(np.max(np.linalg.norm(t-trend,axis=1))))
    stability=max(influence,default=0.);uncertainty=max(formal,stability)
    sigma=options['centroid_sigma'];accuracy=min(options['max_scatter_pix'],max(options['scatter_floor_pix'],sigma*float(np.sqrt(np.mean(e**2)))))
    significant=amplitude>sigma*uncertainty
    exceeds=significant and amplitude-sigma*uncertainty>accuracy
    slopes=coef[1:]/extent[:,None];scale=options.get('pixel_scale_arcsec');factor=1000*scale if scale is not None else None
    out.update(status='exceeds_accuracy' if exceeds else ('significant_within_accuracy' if significant else 'consistent_with_measurement_noise'),
        significant=bool(significant),exceeds_accuracy=bool(exceeds),amplitude_pix=amplitude,
        amplitude_mas=amplitude*factor if factor is not None else None,amplitude_uncertainty_pix=uncertainty,
        amplitude_uncertainty_mas=uncertainty*factor if factor is not None else None,
        systematic_rms_pix=float(np.sqrt(np.mean(np.sum(trend**2,axis=1)))),accuracy_limit_pix=accuracy,
        confidence_sigma=sigma,amplitude_lower_bound_pix=max(0.,amplitude-sigma*uncertainty),
        covariance_uncertainty_pix=formal,leave_one_out_max_change_pix=stability,leave_one_out_count=len(loo),
        dominant_source_index=int(np.argmax(influence)) if len(influence)==n else None,
        fractional_max_source_influence=stability/amplitude if amplitude>0 else None,
        gradient=slopes.T.tolist(),rotation_like_deg=float(np.rad2deg((slopes[0,1]-slopes[1,0])/2)),
        scale_like_fraction=float((slopes[0,0]+slopes[1,1])/2),
        diagnostic='WARNING: coherent residual variation exceeds accuracy' if exceeds else ('WARNING: measurable spatial pattern; within accuracy' if significant else 'No significant spatial trend at centroid precision'),
        uncertainty_method='Maximum of centroid-error covariance and largest leave-one-source-out trend change; not independent validation')
    return out
