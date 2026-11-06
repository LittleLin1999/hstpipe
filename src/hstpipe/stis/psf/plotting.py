import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np
from .common import *


def detector_plot(r,s,valid,trace,wave,c):
    fig,ax=plt.subplots(figsize=(10,5))
    v=np.nanpercentile(abs(s),99.8);ax.imshow(s,origin='lower',aspect='auto',norm=SymLogNorm(linthresh=max(v/100,.01),vmin=-v,vmax=v),cmap='magma')
    x=np.arange(len(trace));ax.plot(x,trace,color='cyan',lw=1,label='Measured/calibrated trace')
    for d in [-160,-100,100,160]:ax.plot(x,trace+d,'g--',lw=.6)
    # Mark only severe DQ at a modest sampling density for readability.
    yb,xb=np.where(~valid);ax.scatter(xb[::max(1,len(xb)//3000)],yb[::max(1,len(yb)//3000)],s=1,c='red',alpha=.4,label='Serious DQ')
    ax.set(xlabel='Detector column [pixel]',ylabel='Detector row [pixel]',title=f'{r["root"]} {r["target"]} {r["aperture"]} {r["date"]}');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(f'figures/{r["root"]}_detector.png',dpi=120);plt.close(fig)


def diagnostics(measurements):
    m=[r for r in measurements if r['mode']=='aligned' and r['background']==0]
    for band in sorted(set(r['band'] for r in m)):
        group=[r for r in m if r['band']==band]
        fig,axes=plt.subplots(1,3,figsize=(15,4))
        for r in group:
            d=np.load(r['profile_file']);y=d['y'];p=d['normalized'];e=np.sqrt(np.maximum(np.diag(d['normalized_covariance']),0));label=f'{r["root"]} {r["target"]} {r["aperture"]}'
            axes[0].plot(y,p,lw=.7,label=label);axes[1].plot(y,p,lw=.7);axes[1].fill_between(y,p-e,p+e,alpha=.07)
            heights=[3,5,7,9,11,15,21,41,81,140];axes[2].errorbar(np.array(heights)/2,[r[f'ee_{h}'] for h in heights],yerr=[r[f'ee_{h}_error'] for h in heights],lw=.5,markersize=2)
        axes[0].set(xlim=(-12,12),xlabel='Along-slit offset [native pixel]',ylabel='Finite-window flux fraction / pixel',title=band+' Å core')
        axes[1].set(xlim=(-70,70),yscale='symlog',ylim=(-.0005,.5),xlabel='Along-slit offset [native pixel]',ylabel='Normalized profile',title='Wings; shaded ±1σ measurement')
        axes[2].set(xlabel='Aperture half-height [native pixel]',ylabel='EE / flux within ±70 pixels',title='Finite-window encircled energy')
        axes[0].legend(fontsize=4,loc='upper left');fig.tight_layout();fig.savefig(f'figures/profiles_{band}.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for star in sorted(set(r['target'] for r in m)):
        g=[r for r in m if r['target']==star and r['hi']-r['lo']<100]
        ax.scatter([(r['lo']+r['hi'])/2 for r in g],[r['fwhm']*.0246 for r in g],s=10,label=star)
    ax.set(xlabel='Observed wavelength [Å]',ylabel='Direct FWHM [arcsec]',title='Profile FWHM; detector-space overlap-aligned');ax.legend(fontsize=6);fig.tight_layout();fig.savefig('figures/fwhm_wavelength.png',dpi=140);plt.close(fig)
