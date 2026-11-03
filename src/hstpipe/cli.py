"""One thin argument parser for both instrument APIs."""
import argparse,json
from .common.config import ReductionConfig
from .common.inventory import discover
from .common.grouping import groups

def main(argv=None):
    parser=argparse.ArgumentParser(prog='hstpipe');instruments=parser.add_subparsers(dest='instrument',required=True)
    for instrument in ('stis','uvis'):
        commands=instruments.add_parser(instrument).add_subparsers(dest='command',required=True)
        for command in ('inspect','reduce'):
            p=commands.add_parser(command)
            p.add_argument('--input');p.add_argument('--output');p.add_argument('--config')
            p.add_argument('--target');p.add_argument('--grating' if instrument=='stis' else '--filter',dest='element');p.add_argument('--setup')
            p.add_argument('--dry-run',action='store_true');p.add_argument('--force',action='store_true');p.add_argument('--resume',action='store_true')
        if instrument=='uvis':
            p=commands.add_parser('migrate')
            p.add_argument('--output',required=True);p.add_argument('--target');p.add_argument('--dry-run',action='store_true')
            p=commands.add_parser('cleanup',help='Preview or remove regenerable UVIS comparison images')
            p.add_argument('--output',required=True);p.add_argument('--target')
            mode=p.add_mutually_exclusive_group()
            mode.add_argument('--apply',action='store_true',help='Remove the eligible diagnostic files')
            mode.add_argument('--dry-run',action='store_true',help='Preview only (default)')
    args=parser.parse_args(argv)
    if args.command=='cleanup':
        from .common.uvis_cleanup import cleanup_uvis
        print(json.dumps(cleanup_uvis(args.output,target=args.target,apply=args.apply),indent=2));return
    if args.command=='migrate':
        from .uvis.migration import migrate_uvis
        print(json.dumps(migrate_uvis(args.output,target=args.target,dry_run=args.dry_run),indent=2));return
    if args.command=='inspect':
        if not args.input:parser.error('inspect requires --input')
        rows=discover(args.input,args.instrument)
        print(json.dumps({'exposures':[dict(path=str(e.path),product=e.product,metadata=e.metadata) for e in rows],
                         'groups':[dict(target=k[0],element=k[1],count=len(v)) for k,v in groups(rows,args.instrument).items()]},indent=2));return
    if args.config:
        c=ReductionConfig.read(args.config)
        if c.instrument!=args.instrument:parser.error('Config instrument differs from command')
        for key in ('input','output','target','element','setup'):
            value=getattr(args,key)
            if value is not None:setattr(c,key,value)
        c.force=c.force or args.force;c.resume=c.resume or args.resume;c.__post_init__()
    else:
        if not args.input or not args.output:parser.error('reduce requires --config or --input and --output')
        c=ReductionConfig(args.input,args.output,args.instrument,args.target,args.element,args.setup,resume=args.resume,force=args.force)
    if args.instrument=='stis':
        from .stis import reduce_stis as reduce
    else:
        from .uvis import reduce_uvis as reduce
    result=reduce(c,dry_run=args.dry_run);print(json.dumps(result,default=str,indent=2))
