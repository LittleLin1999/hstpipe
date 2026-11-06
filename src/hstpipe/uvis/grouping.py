"""Combine overlapping UVIS footprints without visit/proposal partitions."""
from collections import defaultdict
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

def overlaps(a,b):
    # Convex footprint intersection in a local sky-offset frame. Distortion
    # is included in the corner coordinates; no image-size/aperture assumption.
    for x in a.metadata['geometry']:
        for y in b.metadata['geometry']:
            xx=np.array(x['footprint']);yy=np.array(y['footprint'])
            ca=SkyCoord(xx[:,0]*u.deg,xx[:,1]*u.deg);cb=SkyCoord(yy[:,0]*u.deg,yy[:,1]*u.deg)
            if ca[0].separation(cb[0])>1*u.deg:continue
            frame=ca[0].skyoffset_frame()
            aa=ca.transform_to(frame);bb=cb.transform_to(frame)
            pa=np.column_stack([aa.lon.deg,aa.lat.deg]);pb=np.column_stack([bb.lon.deg,bb.lat.deg])
            separated=False
            for poly in (pa,pb):
                for edge in np.roll(poly,-1,axis=0)-poly:
                    axis=np.array([-edge[1],edge[0]])
                    ap=pa@axis;bp=pb@axis
                    if ap.max()<bp.min() or bp.max()<ap.min():separated=True;break
                if separated:break
            if not separated:return True
    return False

def group_exposures(exposures,field=None):
    base=defaultdict(list);result={}
    for e in exposures:
        if e.product!='flc':continue
        name=field or e.metadata.get('TARGNAME')
        if not name:raise ValueError('Missing TARGNAME: set --field or field in the UVIS configuration')
        from ..common.paths import label
        try:label(name)
        except ValueError as exc:raise ValueError('Set a filesystem-safe --field for these observations') from exc
        base[(name,e.metadata['FILTER'])].append(e)
    for (target,filt),rows in base.items():
        components=[]
        for e in rows:
            touching=[c for c in components if any(overlaps(e,q) for q in c)]
            merged=[e]
            for c in touching:merged+=c;components.remove(c)
            components.append(merged)
        for i,c in enumerate(components):
            import hashlib,json
            fingerprint=json.dumps(sorted(e.metadata['ROOTNAME'] for e in c)).encode()
            setup='field-'+hashlib.sha256(fingerprint).hexdigest()[:8] if len(components)>1 else None
            for e in c:e.metadata['grouping_reason']=('explicit field/filter' if field else 'same TARGNAME/filter')+'; connected WCS footprints; no visit/proposal split'+('; disjoint field separated' if setup else '')
            result[(target,filt,setup)]=c
    return result
