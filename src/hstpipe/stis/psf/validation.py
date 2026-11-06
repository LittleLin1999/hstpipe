from collections import defaultdict
import numpy as np
from scipy.linalg import pinvh
from scipy.signal import fftconvolve
from astropy.table import Table
from .common import *
from .psf import eligible,combine
from .extraction import fractional_aperture
import matplotlib.pyplot as plt


from functools import lru_cache

@lru_cache(maxsize=1024)
def _restricted_metrics(profile_file,radius):
    from .extraction import measures
    d=np.load(profile_file);keep=abs(d['y'])<=radius
    m,_,_=measures(d['y'][keep],d['profile'][keep],d['covariance'][np.ix_(keep,keep)],[11,21])
    return m

def delta_metrics(a,b):
    radius=min(a.get('normalization_radius',70),b.get('normalization_radius',70))
    if a.get('normalization_radius',70)!=b.get('normalization_radius',70):
        from .extraction import measures
        def restricted(item):
            m=_restricted_metrics(item['profile_file'],radius)
            return dict(item,**m)
        a,b=restricted(a),restricted(b)
    out={'comparison_normalization_radius':radius}
    for k in ['fwhm','wing_fraction','asymmetry','ee_11','ee_21']:
        out[k+'_difference']=a[k]-b[k]
        ek=k+'_error' if k!='wing_fraction' else 'wing_error'
        if ek in a and ek in b:
            e=np.hypot(a[ek],b[ek]);out[k+'_error']=float(e);out[k+'_z_measurement_only']=(a[k]-b[k])/e if e else None
    return out


def extension_direction(y,p,scale):
    kernel=np.exp(-abs(y)/(scale/.0246));kernel/=kernel.sum();e=fftconvolve(p,kernel,mode='same');e/=e.sum()
    return e-p


def extension_design(y,p,scale):
    direction=extension_direction(y,p,scale)
    centroid=np.gradient(p);centroid-=p*centroid.sum()
    # Normalizing a source plus a residual constant background couples source
    # amplitude to that background: dP/db = 1 - N*P, not a constant alone.
    background=np.ones(len(p))-len(p)*p
    slope=(y-np.mean(y))/max(abs(y))
    return np.column_stack([direction,centroid,background,slope])


def extension_fit(y,data,p,cov,scale):
    x=extension_design(y,p,scale)
    ridge=max(np.max(np.diag(cov))*1e-7,1e-14);w=pinvh(cov+np.eye(len(p))*ridge,rtol=1e-8)
    mc=pinvh(x.T@w@x,rtol=1e-10);coef=mc@x.T@w@(data-p)
    residual=data-p-x@coef;chi=float(residual@w@residual)
    source=1-len(p)*coef[2]
    if source<=0:return float('nan'),float('nan'),chi
    fraction=coef[0]/source;gradient=np.zeros(4);gradient[0]=1/source;gradient[2]=len(p)*coef[0]/source**2
    return float(fraction),float(np.sqrt(max(gradient@mc@gradient,0))),chi


def validate(c):
    from .drizzle_sensitivity import run as drizzle_test
    drizzle_test(c)
    rows=load_json('data/catalogs/observations.json');m=load_json('products/measurements/profiles.json');good=eligible(m,rows);results={};epoch=[];rect=[];background=[];slit=[];null=[]
    grw=[r for r in good if r['target']=='GRW+70D5824' and r['program'] in [17644,18169]]
    for a in grw:
        if a['program']!=18169:continue
        for b in grw:
            if b['program']==17644 and a['band']==b['band']:
                epoch.append(dict(root_a=a['root'],root_b=b['root'],band=a['band'],trace_y_a=a['trace_y'],trace_y_b=b['trace_y'],**delta_metrics(a,b)))
    for a in good:
        for b in m:
            if a['root'] not in b.get('roots',[b['root']]) or b['band']!=a['band']:continue
            if b['mode'] in ['x2d','drizzled','aligned'] and b['background']==0 and not b['quality']:rect.append(dict(root=a['root'],band=a['band'],comparison='native_minus_'+b['mode'],**delta_metrics(a,b)))
            if b['mode']=='unaligned' and b['background']==0 and not b['quality']:rect.append(dict(root=a['root'],band=a['band'],comparison='aligned_minus_native_column_sum',**delta_metrics(a,b)))
            if b['mode']=='aligned' and b['background']>0 and a['wing_usable']:
                base=next((x for x in m if x['root']==b['root'] and x['band']==b['band'] and x['mode']=='aligned' and x['background']==0),None)
                if base:
                    radius=a.get('normalization_radius',70)
                    ba=dict(base,**_restricted_metrics(base['profile_file'],radius));bb=dict(b,**_restricted_metrics(b['profile_file'],radius))
                    ba['normalization_radius']=radius;bb['normalization_radius']=radius
                    background.append(dict(root=a['root'],band=a['band'],background=b['background'],method='same_pixels_local_model_variation',**delta_metrics(ba,bb)))
    groups=defaultdict(list)
    for r in good:groups[(r['aperture'],r['band'],int(np.rint(r['trace_y']/c.get('detector_position_group_pixels',100))),r.get('normalization_radius',70))].append(r)
    residuals={}
    for key,items in groups.items():
        for star in sorted(set(r['target'] for r in items)):
            train=[r for r in items if r['target']!=star];test=[r for r in items if r['target']==star]
            if not train:continue
            p,cm,pred,scatter,ns=combine(train)
            for r in test:
                d=np.load(r['profile_file']);cov=d['normalized_covariance']+pred;y=d['y'];res=d['normalized']-p
                residuals.setdefault(key,[]).append(res)
                f,e,chi=extension_fit(y,d['normalized'],p,cov,.2)
                null.append(dict(root=r['root'],star=star,aperture=key[0],band=key[1],position_group=key[2],training_stars=ns,normalization_radius=r.get('normalization_radius',70),wing_usable=r['wing_usable'],extended_fraction=f,error=e,measurement_plus_template_z=f/e if e else None,chi2=chi,epoch_transfer_included=True))
    # Slit comparisons are diagnostic only unless the same star and overlapping epoch/position exist.
    for a in good:
        if a['aperture']!='52X0.5' or not a['wing_usable']:continue
        b=[r for r in good if r['aperture']=='52X2' and r['band']==a['band'] and r['wing_usable']]
        for r in b:slit.append(dict(narrow_root=a['root'],wide_root=r['root'],band=a['band'],same_star=a['target']==r['target'],trace_separation=abs(a['trace_y']-r['trace_y']),**delta_metrics(a,r)))
    results.update(epoch=epoch,rectification=rect,background=background,slit=slit,leave_one_star_out=null)
    for key,data in results.items():
        save_json(f'products/measurements/{key}.json',data)
        if data:Table(rows=data).write(f'products/measurements/{key}.ecsv',overwrite=True)
    # Empirical null samples are small; simulation FPR is conditional on the measured library.
    simulations=[];rng=np.random.default_rng(c['seed'])
    for key,items in groups.items():
        if key[0]!='52X2' or key[1]!='1330-1680':continue
        p,cm,pred,scatter,ns=combine(items);d=np.load(items[0]['profile_file']);y=d['y'];residualpool=residuals.get(key,[])
        for snr in [30,100,300,1000]:
            photo_variance=np.maximum(abs(p),1e-5)/snr**2
            jac=np.eye(len(p))-p[:,None]
            noise=(jac*photo_variance[None,:])@jac.T
            cov=noise+pred
            for scale in c['extended_scales_arcsec']:
                direction=extension_direction(y,p,scale)
                # Precompute the GLS linear estimator once per scale/SNR.
                x=extension_design(y,p,scale);w=pinvh(cov+np.eye(len(p))*max(np.max(np.diag(cov))*1e-7,1e-14),rtol=1e-8);coefcov=pinvh(x.T@w@x,rtol=1e-10);estimator_matrix=coefcov@x.T@w
                def score(draws):
                    pars=draws@estimator_matrix.T;source=1-len(p)*pars[:,2]
                    return pars[:,0]/np.maximum(source,1e-6)
                n=c['simulation_trials'];noise_draw=rng.normal(size=(n,len(y)))*np.sqrt(photo_variance)
                noise_draw-=noise_draw.sum(axis=1)[:,None]*p[None,:]
                if residualpool:noise_draw+=np.array(residualpool)[rng.integers(0,len(residualpool),n)]
                else:
                    ev,vec=np.linalg.eigh(pred);noise_draw+=rng.normal(size=(n,len(y)))@(vec*np.sqrt(np.maximum(ev,0))).T
                scores=score(noise_draw);threshold=float(np.quantile(scores,.99));hold=rng.normal(size=(n,len(y)))*np.sqrt(photo_variance)
                hold-=hold.sum(axis=1)[:,None]*p[None,:]
                if residualpool:hold+=np.array(residualpool)[rng.integers(0,len(residualpool),n)]
                else:hold+=rng.normal(size=(n,len(y)))@(vec*np.sqrt(np.maximum(ev,0))).T
                nullscores=score(hold);fpr=float(np.mean(nullscores>threshold))
                for fraction in c['extended_fractions']:
                    truth=p+fraction*direction
                    injected_noise=rng.normal(size=(n,len(y)))*np.sqrt(np.maximum(truth,1e-5))/snr
                    injected_noise-=injected_noise.sum(axis=1)[:,None]*truth[None,:]
                    if residualpool:injected_noise+=np.array(residualpool)[rng.integers(0,len(residualpool),n)]
                    else:injected_noise+=rng.normal(size=(n,len(y)))@(vec*np.sqrt(np.maximum(ev,0))).T
                    recovered=score(fraction*direction+injected_noise);complete=float(np.mean(recovered>threshold))
                    simulations.append(dict(aperture=key[0],position_group=key[2],snr=snr,scale_arcsec=scale,injected_fraction=fraction,bias=float(np.mean(recovered)-fraction),scatter=float(np.std(recovered)),completeness=complete,false_positive_rate=fpr,threshold_fraction=threshold,n_trials=n,independent_training_stars=ns,null_residuals=len(residualpool),noise_model='Gaussian count-limit approximation plus empirical residual resampling',conditional=True))
    save_json('products/measurements/injection_recovery.json',simulations)
    if simulations:Table(rows=simulations).write('products/measurements/injection_recovery.ecsv',overwrite=True)
    make_plots(results,simulations)
    return results


def make_plots(results,sim):
    specs=[('epoch','fwhm_difference','FWHM difference [pixel]'),('rectification','fwhm_difference','FWHM aligned − comparison [pixel]'),('background','wing_fraction_difference','Wing-fraction change'),('slit','ee_11_difference','EE(11) narrow − wide'),('leave_one_star_out','extended_fraction','Apparent extended fraction, 0.2 arcsec')]
    for key,metric,label in specs:
        rows=results[key]
        if not rows:continue
        fig,ax=plt.subplots(figsize=(10,5));x=np.arange(len(rows));ax.scatter(x,[r[metric] for r in rows],s=10)
        if key=='leave_one_star_out':ax.errorbar(x,[r[metric] for r in rows],yerr=[r['error'] for r in rows],fmt='none',lw=.5)
        ax.axhline(0,color='k',lw=.7);ax.set(xlabel='Comparison index (see corresponding ECSV table)',ylabel=label,title=key.replace('_',' '));fig.tight_layout();fig.savefig(f'figures/{key}.png',dpi=140);plt.close(fig)
    if sim:
        fig,ax=plt.subplots(1,2,figsize=(11,4))
        for sn in [30,100,300,1000]:
            g=[r for r in sim if r['snr']==sn and r['scale_arcsec']==.2]
            ax[0].plot([r['injected_fraction'] for r in g],[r['completeness'] for r in g],'o-',label=f'S/N={sn}')
        ax[0].set(xlabel='Injected extended fraction',ylabel='Completeness',ylim=(0,1),title='Scale 0.2 arcsec; conditional simulation');ax[0].legend()
        ax[1].scatter([r['snr'] for r in sim],[r['false_positive_rate'] for r in sim],s=8);ax[1].set(xlabel='Integrated continuum S/N',ylabel='Null false-positive rate',title='Independent noise draws, empirical template scatter')
        fig.tight_layout();fig.savefig('figures/injection_recovery.png',dpi=140);plt.close(fig)
