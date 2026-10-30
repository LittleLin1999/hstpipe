"""Read-only standalone environment check; does not reduce observations."""
import argparse,importlib,json,sys
from pathlib import Path
from importlib.metadata import version
p=argparse.ArgumentParser();p.add_argument('--hstcal',type=Path,default=Path(sys.prefix)/'bin');args=p.parse_args()
result={'python':sys.version,'prefix':sys.prefix,'dependencies':{}}
for package,module in [('numpy','numpy'),('astropy','astropy'),('scipy','scipy'),('matplotlib','matplotlib'),('PyYAML','yaml'),('drizzle','drizzle.resample'),('stistools','stistools.x1d'),('drizzlepac','drizzlepac.astrodrizzle'),('stsci.tools','stsci.tools')]:
 importlib.import_module(module);result['dependencies'][package]=version(package)
for name in ('cs6.e','cs7.e'):
 executable=args.hstcal/name
 if not executable.is_file():raise FileNotFoundError(executable)
 result[name]=str(executable.resolve())
print(json.dumps(result,indent=2))
