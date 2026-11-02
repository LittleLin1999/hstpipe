"""Independent direct Astroquery/pyvo providers. No image fitting or resampling."""
from datetime import datetime,timezone
from functools import partial
from pathlib import Path
import hashlib,json,multiprocessing,time
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord,Distance,CartesianRepresentation
from astropy.io import fits
from astropy.table import Table,MaskedColumn
from astropy.time import Time
from astropy.wcs import WCS
from ...common.cache import StageCache
from ...common.target_runner import save
from ...common.provenance import versions

TAP_URL='https://datalab.noirlab.edu/tap'
GAIA_TAP_URLS=('https://gea.esac.esa.int/tap-server/tap','https://gaia.ari.uni-heidelberg.de/tap')
_FAILED_REQUESTS={}
GAIA_COLUMNS=('source_id ra dec ra_error dec_error pmra pmdec pmra_error pmdec_error parallax parallax_error '
              'ref_epoch phot_g_mean_mag ruwe astrometric_params_solved visibility_periods_used duplicated_source '
              'radial_velocity ra_pmra_corr dec_pmdec_corr').split()
LEGACY_REQUIRED='ls_id ra dec type brick_primary maskbits flux_r flux_z'.split()
LEGACY_OPTIONAL='release photsys ra_ivar dec_ivar flux_ivar_r flux_ivar_z nobs_r nobs_z fracflux_r fracflux_z fracmasked_r fracmasked_z fracin_r fracin_z shape_r ref_epoch ref_cat ref_id pmra pmdec pmra_ivar pmdec_ivar mjd_min mjd_max flux_g flux_ivar_g nobs_g rchisq_g rchisq_r rchisq_z'.split()

class CatalogError(RuntimeError):pass
class CatalogQueryError(CatalogError):pass
class CatalogSchemaError(CatalogError):pass
class CatalogOfflineError(CatalogQueryError):pass


def _epoch(value):
    if isinstance(value,(int,float,np.integer,np.floating)):return Time(value,format='jyear',scale='tcb')
    return Time(value).tcb


def _numbers(ra,dec,radius):
    ra=float(ra);dec=float(dec);radius=float(u.Quantity(radius,u.deg).to_value(u.deg))
    if not np.isfinite([ra,dec,radius]).all() or not -90<=dec<=90 or not 0<radius<180:raise ValueError('Invalid sky cone')
    return ra%360,dec,radius


def gaia_adql(ra,dec,radius):
    ra,dec,radius=_numbers(ra,dec,radius)
    return (f"SELECT {', '.join(GAIA_COLUMNS)} FROM gaiadr3.gaia_source WHERE "
            f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.12f},{dec:.12f},{radius:.12f}))")


def legacy_adql(ra,dec,radius,columns):
    ra,dec,radius=_numbers(ra,dec,radius)
    allowed=set(LEGACY_REQUIRED+LEGACY_OPTIONAL)
    if not set(columns)<=allowed:raise ValueError('Unsupported Tractor columns')
    return (f"SELECT {', '.join(columns)} FROM ls_dr10.tractor WHERE "
            f"'t' = Q3C_RADIAL_QUERY(ra,dec,{ra:.12f},{dec:.12f},{radius:.12f}) AND brick_primary = 1")


def _request(kind,query,tap_url):
    if kind=='gaia_dr3':
        from astroquery.gaia import GaiaClass
        from urllib.parse import urlsplit
        endpoint=GAIA_TAP_URLS[0] if tap_url==TAP_URL else tap_url
        parts=urlsplit(endpoint);context=parts.path.removesuffix('/tap').strip('/')
        client=GaiaClass(gaia_tap_server=f'{parts.scheme}://{parts.netloc}/',tap_server_context=context,show_server_messages=False)
        return client.launch_job_async(query,output_format='votable',dump_to_file=False).get_results()
    import pyvo
    return pyvo.dal.TAPService(tap_url).run_sync(query).to_table()


def _worker(kind,query,tap_url,connection,timeout):
    import socket
    socket.setdefaulttimeout(min(15.,max(1.,timeout/3)))
    try:connection.send((True,_request(kind,query,tap_url)))
    except Exception as exc:connection.send((False,f'{type(exc).__name__}: {exc}'))
    finally:connection.close()


def request(kind,query,*,tap_url=TAP_URL,timeout=60,retries=2,backoff=1):
    """Bound wall-clock network time, without replacing either supported client.

    A disposable process allows cancellation even when a client's underlying HTTP
    read cannot be interrupted. timeout=None is useful for injected unit clients.
    """
    failure_key=(kind,query,tap_url)
    if failure_key in _FAILED_REQUESTS:raise CatalogQueryError('Identical request already failed in this reduction: '+_FAILED_REQUESTS[failure_key])
    for attempt in range(retries+1):
        try:
            if timeout is None:return _request(kind,query,tap_url)
            import sys
            context=multiprocessing.get_context('fork' if sys.platform!='darwin' and 'fork' in multiprocessing.get_all_start_methods() else 'spawn')
            read,write=context.Pipe(duplex=False)
            process=context.Process(target=_worker,args=(kind,query,tap_url,write,timeout));process.start();write.close()
            try:
                if not read.poll(timeout):raise TimeoutError(f'{kind} query exceeded {timeout}s')
                ok,result=read.recv()
                if not ok:raise CatalogQueryError(result)
                return result
            finally:
                read.close()
                if process.is_alive():process.terminate()
                process.join()
        except (KeyboardInterrupt,SystemExit):raise
        except Exception as exc:
            if attempt==retries:
                _FAILED_REQUESTS[failure_key]=f'{kind}: {exc} ({attempt+1} attempts)'
                raise CatalogQueryError(_FAILED_REQUESTS[failure_key]) from exc
            time.sleep(backoff*2**attempt)


def _column(table,name,default=np.nan,dtype=float):
    if name not in table.colnames:return np.full(len(table),default,dtype=dtype)
    try:column=np.ma.asarray(table[name],dtype=float if np.issubdtype(np.dtype(dtype),np.floating) else None)
    except (ValueError,TypeError) as exc:raise CatalogSchemaError(f'Non-numerical or malformed {name}') from exc
    values=column.filled(default) if np.ma.isMaskedArray(column) else column
    if dtype is bool and np.asarray(values).dtype.kind in ('U','S','O'):
        text=np.char.lower(np.asarray(values).astype(str))
        if not np.isin(text,['t','true','1','f','false','0']).all():raise CatalogSchemaError(f'Malformed boolean {name}')
        return np.isin(text,['t','true','1'])
    try:return np.asarray(values,dtype=dtype)
    except (ValueError,TypeError) as exc:raise CatalogSchemaError(f'Non-numerical or malformed {name}') from exc


def _coordinates(table):
    for name in ('ra','dec'):
        if name not in table.colnames:raise CatalogSchemaError(f'Missing {name}')
    ra=_column(table,'ra');dec=_column(table,'dec')
    for name in ('ra','dec'):
        if table[name].unit is not None:
            try:converted=u.Quantity(table[name]).to_value(u.deg)
            except u.UnitConversionError as exc:raise CatalogSchemaError(f'Invalid angular units for {name}') from exc
            if name=='ra':ra=np.asarray(converted,dtype=float)
            else:dec=np.asarray(converted,dtype=float)
    if ra.ndim!=1 or dec.ndim!=1:raise CatalogSchemaError('Coordinates must be scalar numerical columns')
    good=np.isfinite(ra)&np.isfinite(dec)&(ra>=0)&(ra<360)&(abs(dec)<=90)
    good&=~np.ma.getmaskarray(table['ra'])&~np.ma.getmaskarray(table['dec'])
    if len(table) and not good.any():raise CatalogSchemaError('No valid numerical ICRS coordinates in response')
    return ra,dec,good



def _verify_cone(table,ra,dec,radius):
    x,y,good=_coordinates(table)
    if good.any():
        center=SkyCoord(ra*u.deg,dec*u.deg);sky=SkyCoord(x[good]*u.deg,y[good]*u.deg)
        if np.any(center.separation(sky).deg>radius+max(1e-8,radius*1e-6)):
            raise CatalogSchemaError('Returned coordinates are inconsistent with the requested sky cone')


def _isolated(ra,dec,minimum):
    if len(ra)<2:return np.ones(len(ra),bool)
    sky=SkyCoord(ra*u.deg,dec*u.deg)
    _,distance,_=sky.match_to_catalog_sky(sky,nthneighbor=2)
    return distance.to_value(u.arcsec)>=minimum


def standardize_gaia(raw,observation_epoch=None,quality=None):
    quality=dict(quality or {});ra,dec,good=_coordinates(raw)
    for key in ('source_id','ref_epoch','ra_error','dec_error'):
        if key not in raw.colnames:raise CatalogSchemaError(f'Missing Gaia {key}')
    pmra=_column(raw,'pmra');pmdec=_column(raw,'pmdec');epoch=_column(raw,'ref_epoch')
    has_pm=np.isfinite(pmra)&np.isfinite(pmdec)&np.isfinite(epoch)
    good&=np.isfinite(_column(raw,'ra_error'))&np.isfinite(_column(raw,'dec_error'))&(_column(raw,'ra_error')>0)&(_column(raw,'dec_error')>0)
    if observation_epoch is not None and quality.get('apply_proper_motion',True) and not quality.get('allow_missing_pm',False):
        good&=np.isfinite(_column(raw,'pmra_error'))&np.isfinite(_column(raw,'pmdec_error'))&(_column(raw,'pmra_error')>=0)&(_column(raw,'pmdec_error')>=0)
    good&=np.isfinite(epoch)&(_column(raw,'ruwe')<=quality.get('max_ruwe',1.4))
    good&=_column(raw,'visibility_periods_used',0)>=quality.get('min_visibility_periods',8)
    good&=np.isin(_column(raw,'astrometric_params_solved',0,int),[31,95])
    good&=~_column(raw,'duplicated_source',False,bool)
    if 'duplicated_source' in raw.colnames:good&=~np.ma.getmaskarray(raw['duplicated_source'])
    good&=_column(raw,'phot_g_mean_mag')<=quality.get('max_mag',21)
    if observation_epoch is not None and quality.get('apply_proper_motion',True) and not quality.get('allow_missing_pm',False):good&=has_pm
    good&=~np.ma.getmaskarray(raw['source_id'])
    ids=np.asarray(raw['source_id']).astype(str);_,first=np.unique(ids,return_index=True)
    unique=np.zeros(len(raw),bool);unique[first]=True;good&=unique
    selected=raw[good];ra=ra[good].copy();dec=dec[good].copy();epoch=epoch[good];pmra=pmra[good];pmdec=pmdec[good];has_pm=has_pm[good]
    native_ra=ra.copy();native_dec=dec.copy()
    target=_epoch(observation_epoch) if observation_epoch is not None else None
    position_epoch=epoch.copy();methods=np.full(len(selected),'not_propagated',dtype='U64')
    dt=np.zeros(len(selected));parallax=_column(selected,'parallax');pe=_column(selected,'parallax_error')
    distance_ok=np.isfinite(parallax)&np.isfinite(pe)&(parallax>0)&(pe>0)&(parallax/np.maximum(pe,1e-30)>=quality.get('min_parallax_snr',5))
    if target is not None and quality.get('apply_proper_motion',True):
        for with_distance in (False,True):
            use=has_pm&(distance_ok==with_distance)
            if not use.any():continue
            kwargs=dict(ra=ra[use]*u.deg,dec=dec[use]*u.deg,pm_ra_cosdec=pmra[use]*u.mas/u.yr,
                        pm_dec=pmdec[use]*u.mas/u.yr,obstime=Time(epoch[use],format='jyear',scale='tcb'))
            if with_distance:
                kwargs['distance']=Distance(parallax=parallax[use]*u.mas)
                rv=_column(selected,'radial_velocity')[use]
                kwargs['radial_velocity']=np.where(np.isfinite(rv),rv,0)*u.km/u.s
            moved=SkyCoord(**kwargs).apply_space_motion(new_obstime=target)
            ra[use]=moved.ra.deg;dec[use]=moved.dec.deg;position_epoch[use]=target.jyear;dt[use]=target.jyear-epoch[use]
            methods[use]=np.where(np.isfinite(rv),'space_motion_distance_rv','space_motion_distance_rv_assumed_zero') if with_distance else 'angular_space_motion_no_distance'
    out=Table();out['source_id']=np.asarray(selected['source_id']).astype(str)
    out['ra']=ra*u.deg;out['dec']=dec*u.deg
    for axis,pm in (('ra','pmra'),('dec','pmdec')):
        error=_column(selected,axis+'_error');pmerror=_column(selected,pm+'_error')
        rho=_column(selected,axis+'_'+pm+'_corr',0)
        variance=error**2+(dt*np.where(np.isfinite(pmerror),pmerror,0))**2+2*dt*rho*error*np.where(np.isfinite(pmerror),pmerror,0)
        out[axis+'_error_mas']=np.sqrt(np.maximum(variance,0))*u.mas
    out['ra_reference']=native_ra*u.deg;out['dec_reference']=native_dec*u.deg
    for key in ('pmra_error','pmdec_error','ra_pmra_corr','dec_pmdec_corr','radial_velocity','astrometric_params_solved','visibility_periods_used','duplicated_source'):
        out[key]=_column(selected,key,False if key=='duplicated_source' else np.nan,bool if key=='duplicated_source' else float)
    out['ra_error_reference']=_column(selected,'ra_error')*u.mas;out['dec_error_reference']=_column(selected,'dec_error')*u.mas
    out['reference_epoch']=epoch;out['position_epoch']=position_epoch;out['pm_applied']=methods!='not_propagated';out['propagation_method']=methods
    out['pmra']=pmra*u.mas/u.yr;out['pmdec']=pmdec*u.mas/u.yr
    for key in ('phot_g_mean_mag','ruwe','parallax','parallax_error'):out[key]=_column(selected,key)
    out['catalog']=np.full(len(out),'gaia');out['release']=np.full(len(out),'DR3');out['survey_region']=np.full(len(out),'all_sky')
    isolated=_isolated(ra,dec,quality.get('isolation_arcsec',.5));out=out[isolated]
    out.meta.update(epoch_scale='TCB',raw_count=len(raw),quality_rejected=int((~good).sum()),isolation_rejected=int((~isolated).sum()),quality_cuts=quality,
        observation_epoch=target.isot if target is not None else None,uncertainty_method='position/PM covariance in each tangent axis; no full perspective-error covariance')
    return out


def standardize_legacy(raw,quality=None):
    quality=dict(quality or {});ra,dec,good=_coordinates(raw)
    for name in LEGACY_REQUIRED:
        if name not in raw.colnames:raise CatalogSchemaError(f'Missing Tractor {name}')
    good&=_column(raw,'brick_primary',False,bool)
    good&=~np.ma.getmaskarray(raw['type'])&~np.ma.getmaskarray(raw['ls_id'])&~np.ma.getmaskarray(raw['maskbits'])
    kind=np.char.strip(np.asarray(raw['type']).astype(str));good&=np.isin(kind,quality.get('types',['PSF','REX']))
    mask=_column(raw,'maskbits',0,np.int64)
    linked=np.char.startswith(np.asarray(raw['ref_cat']).astype(str),'G') if 'ref_cat' in raw.colnames else np.zeros(len(raw),bool)
    identity_prior=bool(quality.get('allow_masked_stellar_priors',False))&(kind=='PSF')&linked&(_column(raw,'ref_id',0)>0)
    mask_ok=(mask&quality.get('reject_maskbits',2+4+8+16+32+64+256+512+1024+2048+4096))==0
    good&=mask_ok|identity_prior
    snrs=[];coverage=[]
    for band in ('r','z'):
        flux=_column(raw,'flux_'+band);ivar=_column(raw,'flux_ivar_'+band)
        snrs.append(np.where(np.isfinite(ivar)&(ivar>0),flux*np.sqrt(np.maximum(ivar,0)),np.nan))
        if 'nobs_'+band in raw.colnames:coverage.append(_column(raw,'nobs_'+band,0)>0)
    snr=np.nanmax(np.array(snrs),axis=0) if len(raw) else np.array([])
    good&=np.isfinite(snr)&(snr>=quality.get('min_snr',10))
    if coverage:good&=np.any(coverage,axis=0)
    for band in ('r','z'):
        if 'fracflux_'+band in raw.colnames:good&=_column(raw,'fracflux_'+band,0)<=quality.get('max_fracflux',.5)
        if 'fracmasked_'+band in raw.colnames:good&=(_column(raw,'fracmasked_'+band,0)<=quality.get('max_fracmasked',.2))|identity_prior
    if 'shape_r' in raw.colnames:good&=(kind=='PSF')|(_column(raw,'shape_r')<=quality.get('max_shape_arcsec',1.0))
    ids=np.asarray(raw['ls_id']).astype(str);_,first=np.unique(ids,return_index=True);unique=np.zeros(len(raw),bool);unique[first]=True;good&=unique
    selected=raw[good];out=Table();out['source_id']=ids[good];out['ra']=ra[good]*u.deg;out['dec']=dec[good]*u.deg
    out['ra_reference']=ra[good]*u.deg;out['dec_reference']=dec[good]*u.deg
    out['identification_prior_only']=identity_prior[good]
    out['prior_quality_warning']=np.where(identity_prior[good],'Gaia-linked PSF retained for identification despite survey masking; HST centroids still require validation','')
    out['type']=kind[good];out['maskbits']=mask[good];out['snr']=snr[good]
    system=np.char.upper(np.char.strip(np.asarray(np.ma.asarray(selected['photsys']).filled('?')).astype(str))) if 'photsys' in selected.colnames else np.full(len(selected),'?')
    codes=_column(selected,'release',-1,int)
    system=np.where((system=='?')&np.isin(codes,[10000,10002]),'S',system)
    system=np.where((system=='?')&(codes==9011),'N',system)
    out['catalog']=np.full(len(out),'legacy');out['release']=np.where(system=='N','DR9-North',np.where(system=='S','DR10-South','DR10-family-unresolved'))
    out['survey_region']=np.where(system=='N','north',np.where(system=='S','south','unknown'))
    out['release_code']=codes
    for axis in ('ra','dec'):
        ivar=_column(selected,axis+'_ivar');error=np.where(ivar>0,1/np.sqrt(np.maximum(ivar,1e-300))*3600000,np.nan)
        if axis=='ra':error*=np.cos(np.deg2rad(dec[good]))
        out[axis+'_error_mas']=error*u.mas
    native_epoch=_column(selected,'ref_epoch');valid_epoch=np.isfinite(native_epoch)&(native_epoch>1900)&(native_epoch<2100)
    out['reference_epoch']=MaskedColumn(native_epoch,mask=~valid_epoch);out['position_epoch']=MaskedColumn(native_epoch,mask=~valid_epoch)
    for key in ('pmra','pmdec','pmra_ivar','pmdec_ivar','mjd_min','mjd_max'):out[key]=_column(selected,key)
    out['ref_cat']=np.asarray(selected['ref_cat']).astype(str) if 'ref_cat' in selected.colnames else np.full(len(selected),'unknown')
    for key in LEGACY_REQUIRED+LEGACY_OPTIONAL:
        if key in selected.colnames and key not in out.colnames:out[key]=selected[key]
    isolated=_isolated(ra[good],dec[good],quality.get('isolation_arcsec',1.5));out=out[isolated]
    out.meta.update(raw_count=len(raw),quality_rejected=int((~good).sum()),isolation_rejected=int((~isolated).sum()),quality_cuts=quality,
        coverage_source_count=int(np.count_nonzero(np.any(coverage,axis=0))) if coverage else None,
        observation_epoch=None,epoch_caveat='Tractor centroid epochs unavailable here; DR9-North is not newly reprocessed DR10; allow systematic/epoch uncertainty')
    return out


_ACTIVE_CACHE_LOCKS=set()

def _cached(catalog,parameters,loader,cache_dir,refresh,offline,allow_offline_load=False):
    if cache_dir is None or str(Path(cache_dir).resolve()) in _ACTIVE_CACHE_LOCKS:
        return _cached_unlocked(catalog,parameters,loader,cache_dir,refresh,offline,allow_offline_load)
    from ...common.target_runner import target_lock
    # A target's shared cache is small; hold this lock through the request to
    # avoid duplicate downloads and simultaneous reads of partially written files.
    with target_lock(Path(cache_dir)):
        token=str(Path(cache_dir).resolve());_ACTIVE_CACHE_LOCKS.add(token)
        try:return _cached_unlocked(catalog,parameters,loader,cache_dir,refresh,offline,allow_offline_load)
        finally:_ACTIVE_CACHE_LOCKS.remove(token)

def _cached_unlocked(catalog,parameters,loader,cache_dir,refresh,offline,allow_offline_load=False):
    query_id=hashlib.sha256(json.dumps(parameters,sort_keys=True,default=str).encode()).hexdigest()[:20]
    if cache_dir is None:
        if offline:raise CatalogOfflineError(f'{catalog}: no cache directory supplied')
        result,metadata=loader();result.meta.update(metadata);return result
    directory=Path(cache_dir);table=directory/f'{catalog}_{query_id}.ecsv';meta=directory/f'{catalog}_{query_id}.json'
    cache=StageCache(directory,prefix='catalog_cache_');key=cache.key([],parameters)
    if not refresh and cache.valid(query_id,key):
        result=Table.read(table,format='ascii.ecsv');result.meta['cache_file']=str(table);return result
    if offline and not allow_offline_load:raise CatalogOfflineError(f'{catalog}: no valid cached query {query_id}')
    result,metadata=loader();directory.mkdir(parents=True,exist_ok=True)
    metadata.update(catalog=catalog,parameters=parameters,source_count=len(result),versions=versions(),normalized_utc=datetime.now(timezone.utc).isoformat())
    metadata.setdefault('retrieved_utc',datetime.now(timezone.utc).isoformat())
    result.meta.update(metadata);temporary=table.with_suffix('.tmp');result.write(temporary,format='ascii.ecsv',overwrite=True);temporary.replace(table)
    save(meta,metadata);cache.record(query_id,key,[table,meta]);result.meta['cache_file']=str(table)
    return result


def query_gaia_dr3(ra,dec,radius,observation_epoch=None,*,quality=None,cache_dir=None,refresh=False,offline=False,timeout=60,retries=2,tap_urls=None):
    ra,dec,radius=_numbers(ra,dec,radius);query=gaia_adql(ra,dec,radius)
    endpoints=tuple(tap_urls or GAIA_TAP_URLS)
    epoch=_epoch(observation_epoch).isot if observation_epoch is not None else None
    parameters=dict(normalizer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),provider='astroquery.gaia',endpoints=endpoints,release='DR3',query=query,center=[ra,dec],radius_deg=radius,observation_epoch=epoch,quality=quality or {})
    def load():
        raw_parameters=dict(provider='astroquery.gaia',query=query,release='DR3',endpoints=endpoints)
        def raw_load():
            failures=[]
            for endpoint in endpoints:
                try:
                    table=request('gaia_dr3',query,tap_url=endpoint,timeout=None if timeout is None else timeout/len(endpoints),retries=retries)
                    _verify_cone(table,ra,dec,radius)
                    return table,dict(query=query,endpoint=endpoint,attempted_endpoints=endpoints,previous_endpoint_failures=failures)
                except CatalogQueryError as exc:failures.append(str(exc))
            raise CatalogQueryError('; '.join(failures))
        raw=_cached('gaia_dr3_raw',raw_parameters,raw_load,cache_dir,refresh,offline)
        _verify_cone(raw,ra,dec,radius)
        return standardize_gaia(raw,observation_epoch,quality),dict(query=query,endpoint=raw.meta.get('endpoint',endpoints[0]),release='DR3',retrieved_utc=raw.meta.get('retrieved_utc',datetime.now(timezone.utc).isoformat()))
    return _cached('gaia_dr3',parameters,load,cache_dir,refresh,offline,allow_offline_load=True)


def query_legacy_dr10(ra,dec,radius,*,quality=None,cache_dir=None,refresh=False,offline=False,timeout=60,retries=2,tap_url=TAP_URL):
    ra,dec,radius=_numbers(ra,dec,radius)
    parameters=dict(normalizer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),provider='pyvo',release='DR10-South/DR9-North',endpoint=tap_url,table='ls_dr10.tractor',center=[ra,dec],radius_deg=radius,quality=quality or {})
    def load():
        schema_query="SELECT column_name, datatype FROM TAP_SCHEMA.columns WHERE table_name='ls_dr10.tractor'"
        schema=_cached('legacy_dr10_schema',dict(query=schema_query,endpoint=tap_url),lambda:(request('legacy_schema',schema_query,tap_url=tap_url,timeout=timeout,retries=retries),dict(query=schema_query,endpoint=tap_url)),cache_dir,refresh,offline)
        if 'column_name' not in schema.colnames:raise CatalogSchemaError('TAP_SCHEMA missing column_name')
        available={str(c).lower() for c in schema['column_name']}
        missing=set(LEGACY_REQUIRED)-available
        if missing:raise CatalogSchemaError(f'Tractor schema missing {sorted(missing)}')
        columns=LEGACY_REQUIRED+[c for c in LEGACY_OPTIONAL if c in available]
        query=legacy_adql(ra,dec,radius,columns)
        def raw_load():
            table=request('legacy_dr10',query,tap_url=tap_url,timeout=timeout,retries=retries);_verify_cone(table,ra,dec,radius)
            return table,dict(query=query,endpoint=tap_url)
        raw=_cached('legacy_dr10_raw',dict(query=query,endpoint=tap_url),raw_load,cache_dir,refresh,offline)
        _verify_cone(raw,ra,dec,radius)
        return standardize_legacy(raw,quality),dict(query=query,schema_query=schema_query,available_columns=sorted(available),endpoint=tap_url,release='DR10-South/DR9-North',retrieved_utc=raw.meta.get('retrieved_utc',datetime.now(timezone.utc).isoformat()))
    return _cached('legacy_dr10',parameters,load,cache_dir,refresh,offline,allow_offline_load=True)


def query_footprint(image,margin_arcsec=5):
    """Spherical bounding cone from actual covered WCS edges, including RA wrap."""
    with fits.open(image) as hd:
        h=hd['SCI'];w=WCS(h.header,hd,key=h.header.get('HPORIGK',' '));ny,nx=h.shape
        if 'WHT' in hd:
            covered=hd['WHT'].data>0
            xx=np.flatnonzero(np.any(covered,axis=0));yy=np.flatnonzero(np.any(covered,axis=1))
            if not len(xx):raise ValueError('No covered footprint')
            x0,x1=int(xx[0]),int(xx[-1]);y0,y1=int(yy[0]),int(yy[-1])
        else:x0,x1,y0,y1=0,nx-1,0,ny-1
        fraction=np.linspace(0,1,17);x=x0+(x1-x0)*fraction;y=y0+(y1-y0)*fraction
        pixels=np.vstack([np.column_stack([x,np.full(17,y0)]),np.column_stack([x,np.full(17,y1)]),
                          np.column_stack([np.full(17,x0),y]),np.column_stack([np.full(17,x1),y])])
        world=w.all_pix2world(pixels,0);sky=SkyCoord(world[:,0]*u.deg,world[:,1]*u.deg)
        cart=sky.cartesian.xyz.mean(axis=1);center=SkyCoord(CartesianRepresentation(cart)).spherical
        coord=SkyCoord(center.lon,center.lat);radius=max(coord.separation(sky)).to_value(u.deg)+margin_arcsec/3600
        if not np.isfinite(radius) or radius<=0:raise ValueError('Invalid footprint WCS')
        return dict(ra=float(coord.ra.deg),dec=float(coord.dec.deg),radius=float(radius),margin_arcsec=float(margin_arcsec))


def gaia_at_epoch(reference,epoch):
    """Evaluate cached native Gaia astrometry at another actual HST epoch."""
    raw=reference.copy()
    raw['ra']=reference['ra_reference'];raw['dec']=reference['dec_reference']
    raw['ra_error']=reference['ra_error_reference'];raw['dec_error']=reference['dec_error_reference']
    raw['ref_epoch']=reference['reference_epoch']
    return standardize_gaia(raw,epoch,reference.meta.get('quality_cuts',{}))


def legacy_at_epoch(reference,epoch):
    """Propagate only Tractor stars with an explicit epoch and measured motion.

    Unanchored coadd centroids retain unknown effective epochs; their uncertainty
    is not disguised as Gaia DR3 precision.
    """
    out=reference.copy();native=np.ma.asarray(out['reference_epoch'])
    years=np.asarray(native.filled(np.nan),float)
    pmra=_column(out,'pmra');pmdec=_column(out,'pmdec');ivra=_column(out,'pmra_ivar');ivdec=_column(out,'pmdec_ivar')
    anchored=np.char.startswith(np.asarray(out['ref_cat']).astype(str),'G')
    good=anchored&np.isfinite(years)&np.isfinite(pmra)&np.isfinite(pmdec)&(ivra>0)&(ivdec>0)
    out['pm_applied']=np.zeros(len(out),bool)
    if good.any():
        sky=SkyCoord(ra=np.asarray(out['ra_reference'])[good]*u.deg,dec=np.asarray(out['dec_reference'])[good]*u.deg,
            pm_ra_cosdec=pmra[good]*u.mas/u.yr,pm_dec=pmdec[good]*u.mas/u.yr,obstime=Time(years[good],format='jyear',scale='tcb'))
        moved=sky.apply_space_motion(new_obstime=Time(epoch,format='jyear',scale='tcb'))
        out['ra'][good]=moved.ra.deg;out['dec'][good]=moved.dec.deg
        out['position_epoch'][good]=epoch;out['pm_applied'][good]=True
        for axis,ivar in [('ra',ivra),('dec',ivdec)]:
            error=np.asarray(out[axis+'_error_mas'])[good]
            out[axis+'_error_mas'][good]=np.sqrt(error**2+(epoch-years[good])**2/ivar[good])
    out.meta['epoch_caveat']='Measured Tractor motion applied only for Gaia-associated stars with explicit reference epoch and finite inverse variance; other coadd centroids retain unknown effective epoch'
    return out


def cached_footprint(catalog,ra,dec,radius,*,cache_dir,quality=None,observation_epoch=None):
    """Explicit same-provider raw-cache recovery, including partial coverage.

    Never expands a cached cone or treats missing objects outside it as absent.
    Receipts verify the original raw table. This does not silently switch survey.
    """
    import re
    ra,dec,radius=_numbers(ra,dec,radius);directory=Path(cache_dir);choices=[]
    pattern=(r"Q3C_RADIAL_QUERY\(\s*ra\s*,\s*dec\s*,\s*([0-9eE+.-]+)\s*,\s*([0-9eE+.-]+)\s*,\s*([0-9eE+.-]+)" if catalog=='legacy_dr10' else r"CIRCLE\('ICRS',\s*([0-9eE+.-]+)\s*,\s*([0-9eE+.-]+)\s*,\s*([0-9eE+.-]+)")
    for path in directory.glob(catalog+'_raw_*.ecsv'):
        try:
            meta=json.loads(path.with_suffix('.json').read_text());query=meta.get('query',meta.get('parameters',{}).get('query',''))
            parsed=re.search(pattern,query,re.IGNORECASE)
            if not parsed:continue
            x,y,r=map(float,parsed.groups());sep=SkyCoord(ra*u.deg,dec*u.deg).separation(SkyCoord(x*u.deg,y*u.deg)).deg
            if sep>r+radius:continue
            query_id=path.stem.rsplit('_',1)[1];receipt=json.loads(StageCache(directory,prefix='catalog_cache_').receipt(query_id).read_text())
            if not StageCache(directory,prefix='catalog_cache_').valid(query_id,receipt['key']):continue
            raw=Table.read(path,format='ascii.ecsv')
            table=standardize_legacy(raw,quality) if catalog=='legacy_dr10' else standardize_gaia(raw,observation_epoch,quality)
            if len(table):
                sky=SkyCoord(np.asarray(table['ra'])*u.deg,np.asarray(table['dec'])*u.deg)
                table=table[SkyCoord(ra*u.deg,dec*u.deg).separation(sky).deg<=radius]
            complete=bool(sep+radius<=r+1e-9)
            choices.append((complete,len(table),r,path,table,meta,query,[x,y,r]))
        except Exception:continue
    if not choices:raise CatalogOfflineError(f'{catalog}: no verified overlapping raw cache')
    complete,_,_,path,table,metadata,query,cone=max(choices,key=lambda c:(c[0],c[1],c[2]))
    table.meta.update(cache_fallback=True,complete_query_coverage=complete,coverage_caveat='Counts are lower bounds when cached coverage is partial; unavailable objects are not nondetections',
        original_cache_file=str(path),original_query=query,original_cone=cone,requested_cone=[ra,dec,radius],retrieved_utc=metadata.get('retrieved_utc'),endpoint=metadata.get('endpoint'),observation_epoch=observation_epoch)
    ident=hashlib.sha256(json.dumps(dict(source=str(path),query=query,requested=[ra,dec,radius],quality=quality,epoch=observation_epoch,code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),sort_keys=True).encode()).hexdigest()[:20]
    out=directory/f'{catalog}_identification_subset_{ident}.ecsv'
    import tempfile,os
    fd,name=tempfile.mkstemp(dir=directory,suffix='.ecsv');os.close(fd)
    temporary=Path(name)
    try:table.write(temporary,format='ascii.ecsv',overwrite=True);temporary.replace(out)
    finally:temporary.unlink(missing_ok=True)
    table.meta['cache_file']=str(out)
    return table
