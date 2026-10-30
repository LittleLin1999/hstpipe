"""Read-only portable STIS product/checksum regression check."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
p.add_argument('--baseline',type=Path);args=p.parse_args()
manifest=json.loads((args.run/'manifest.json').read_text());report={}
for record in manifest['products']:
    path=args.run/'products'/Path(record['path']).name
    assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'],path
    with fits.open(path,checksum=True) as h:
        h.verify('exception');assert all(a.verify_checksum()==1 and a.verify_datasum()==1 for a in h)
        comparison={}
        if args.baseline:
            with fits.open(args.baseline/'products'/path.name) as before:
                for name in ('SCI','ERR','DQ','NEXP','EXPMAP','BKG','BKGERR'):
                    if name in h and name in before:
                        equal=bool(np.array_equal(h[name].data,before[name].data,equal_nan=True))
                        comparison[name]=equal;assert equal,(path,name)
                for name in ('FLT','DRZ','X2D'):
                    if name not in h:continue
                    for field in ('FLUX_BOX','FLUX_OPT','ERR_BOX','ERR_OPT'):
                        a=h[name].data[field];b=before[name].data[field]
                        assert np.array_equal(np.isfinite(a),np.isfinite(b))
                        good=np.isfinite(a)&np.isfinite(b)
                        error=float(np.max(abs(a[good]-b[good]))/max(np.max(abs(b[good])),1e-30))
                        comparison[name+'/'+field]=error
                        assert error<3e-6,(path,name,field,error)
    report[path.name]=dict(checksums_valid=True,comparison=comparison)
print(json.dumps(report,indent=2))
