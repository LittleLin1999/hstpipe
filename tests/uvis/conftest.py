from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

@pytest.fixture
def flc(tmp_path):
    def make(name='test',chips=(7,),size=128,ra=10,sub=True,filt='F275W'):
        p=tmp_path/f'{name}_flc.fits';hdr=fits.Header()
        for k,v in dict(INSTRUME='WFC3',DETECTOR='UVIS',OBSTYPE='IMAGING',TARGNAME='TEST',FILTER=filt,EXPTIME=100.,ROOTNAME=name,
                        RA_TARG=ra,DEC_TARG=20.,EXPSTART=58000.,DATE_OBS='2017-09-04',SUBARRAY=sub,CCDAMP='C',APERTURE='UVIS2-C1K1C-SUB',FGSLOCK='FINE',EXPFLAG='NORMAL').items():hdr[k]=v
        for k in ('BIASCORR','BLEVCORR','DARKCORR','FLATCORR','PCTECORR','PHOTCORR','DQICORR'):hdr[k]='COMPLETE'
        hd=[fits.PrimaryHDU(header=hdr)];rng=np.random.default_rng(11)
        for ver in chips:
            w=WCS(naxis=2);w.wcs.crpix=[size/2,size/2];w.wcs.cdelt=[-0.04/3600,0.04/3600];w.wcs.crval=[ra,20];w.wcs.ctype=['RA---TAN','DEC--TAN']
            from astropy.wcs import Sip
            w.wcs.ctype=['RA---TAN-SIP','DEC--TAN-SIP'];w.sip=Sip(np.zeros((3,3)),np.zeros((3,3)),None,None,w.wcs.crpix)
            yy,xx=np.indices((size,size));sci=rng.normal(5,1,(size,size))+100*np.exp(-((xx-size/2)**2+(yy-size/2)**2)/12)
            for label,data in [('DQ',np.zeros((size,size),'u2')),('SCI',sci.astype('f4')),('ERR',np.ones((size,size),'f4'))]:
                h=fits.ImageHDU(data,name=label);h.header['EXTVER']=ver
                if label=='SCI':
                    h.header.update(w.to_header(relax=True));h.header['BUNIT']='ELECTRONS';h.header['CCDCHIP']=2;h.header['PHOTFLAM']=1e-18;h.header['PHOTPLAM']=2750.
                hd.append(h)
        fits.HDUList(hd).writeto(p);return p
    return make


@pytest.fixture(autouse=True)
def no_live_identification_queries(monkeypatch):
    from hstpipe.uvis.astrometry import identification
    from hstpipe.uvis.astrometry.catalogs import CatalogQueryError
    def offline(*a,**kw):raise CatalogQueryError('Unit tests use mocked external-identification catalogs')
    monkeypatch.setattr(identification,'query_gaia_dr3',offline)
    monkeypatch.setattr(identification,'query_legacy_dr10',offline)


@pytest.fixture(autouse=True)
def fresh_catalog_failure_cache():
    from hstpipe.uvis.astrometry.catalogs import _FAILED_REQUESTS
    _FAILED_REQUESTS.clear()
    yield
    _FAILED_REQUESTS.clear()
