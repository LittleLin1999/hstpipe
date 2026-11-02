import numpy as np
import pytest
from astropy.io import fits
from hstpipe.uvis.astrometry.centroids import pixel_masks,validate_centroids
from astropy.table import Table

@pytest.mark.parametrize('version',[1,2])
def test_inherited_cr_core_recovered_with_other_flags_preserved(tmp_path,version):
 y,x=np.indices((51,51));s=5+1000*np.exp(-((x-25)**2+(y-25)**2)/(2*1.3**2));d=np.zeros(s.shape,'u2');d[22:29,22:29]=4096
 hd=fits.HDUList([fits.PrimaryHDU()])
 for v in [1,2]:
  for name,data in [('SCI',s.astype('f4')),('ERR',np.ones(s.shape,'f4')),('DQ',d.copy())]:
   h=fits.ImageHDU(data,name=name);h.header['EXTVER']=v;hd.append(h)
 original=hd['DQ',version].data.copy();dq,base,mask=pixel_masks(hd,version,{'edge_pixels':3})
 assert not base[25,25] and not mask[25,25]
 assert np.array_equal(dq,original)
 p=tmp_path/'test.fits';hd.writeto(p)
 cat=Table(dict(x=[25.],y=[25.],extver=[version],snr=[100.]))
 accepted=validate_centroids(p,cat,{'edge_pixels':3})
 assert len(accepted)==1
 for bit in [1<<k for k in range(16) if k!=12]:
  hd['DQ',version].data[25,25]=4096|bit
  _,b,m=pixel_masks(hd,version,{'edge_pixels':0})
  assert b[25,25] and m[25,25]
  if bit!=256:assert not m[25,26]  # ordinary detector flags do not grow

@pytest.mark.parametrize('run_cr,expected',[(True,4096),(False,0)])
def test_split_stage_resetbits(monkeypatch,run_cr,expected):
 from drizzlepac import astrodrizzle
 from hstpipe.uvis.drizzle import drizzle
 seen={}
 monkeypatch.setattr(astrodrizzle,'AstroDrizzle',lambda inputs,**kw:seen.update(kw))
 drizzle(['working_a_flc.fits','working_b_flc.fits'],'test',{'driz_cr':run_cr})
 assert seen['resetbits']==expected
