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
            if instrument=='uvis':p.add_argument('--field',help='Name one field/mosaic from the input exposures, independent of TARGNAME')
            p.add_argument('--target');p.add_argument('--grating' if instrument=='stis' else '--filter',dest='element');p.add_argument('--setup')
            p.add_argument('--dry-run',action='store_true');p.add_argument('--force',action='store_true');p.add_argument('--resume',action='store_true')
        if instrument=='stis':
            p=commands.add_parser('psf',help='Build empirical calibration-star PSF products')
            p.add_argument('--config',required=True)
            p.add_argument('--stage',choices=['all','discover','preview','download','recalibrate','screen','reduce','extract','library','validate','report','plot'],default='all')
        if instrument=='stis':
            p=commands.add_parser('psf-collapse',help='Collapse explicitly selected empirical PSFs onto a requested spatial grid')
            p.add_argument('--library',required=True)
            p.add_argument('--request',required=True)
            p.add_argument('--output',required=True)
            p=commands.add_parser('acq-profile',help='Collapse optical acquisition target images across the nominal slit width')
            p.add_argument('--metadata',required=True)
            p.add_argument('--output',required=True)
            p.add_argument('--spatial-limit',type=float,default=2.)
            p=commands.add_parser('acq',help='Plot acquisition images with nominal science-slit footprints')
            p.add_argument('--acquisition',nargs='+',required=True)
            p.add_argument('--science',nargs='+',required=True)
            p.add_argument('--output',required=True)
            p.add_argument('--cutout-size',type=float,default=5.,help='Target-centered square field width in arcseconds')
        if instrument=='uvis':
            p=commands.add_parser('report',help='Refresh index.html from saved reduction records without reducing images')
            p.add_argument('--directory',required=True,help='Target UVIS directory containing manifest.json')
            p=commands.add_parser('migrate')
            p.add_argument('--output',required=True);p.add_argument('--target');p.add_argument('--dry-run',action='store_true')
            p=commands.add_parser('cleanup',help='Preview or remove regenerable UVIS comparison images')
            p.add_argument('--output',required=True);p.add_argument('--target')
            mode=p.add_mutually_exclusive_group()
            mode.add_argument('--apply',action='store_true',help='Remove the eligible diagnostic files')
            mode.add_argument('--dry-run',action='store_true',help='Preview only (default)')
            p=commands.add_parser('review',help='Inspect or approve recorded UVIS field-source QA')
            p.add_argument('--directory',required=True);p.add_argument('--filter',required=True)
            p.add_argument('--approve',action='store_true');p.add_argument('--reviewer');p.add_argument('--reason')
    args=parser.parse_args(argv)
    if args.command=='psf-collapse':
        from .stis.psf.collapse import collapse_library_profile
        collapse_library_profile(args.library,args.request,args.output);return
    if args.command=='acq-profile':
        from .stis.acquisition_profiles import collapse_acquisitions
        collapse_acquisitions(args.metadata,args.output,spatial_limit=args.spatial_limit);return
    if args.command=='acq':
        from .stis.acquisition import plot_acquisitions
        plot_acquisitions(args.acquisition,args.science,args.output,cutout_size_arcsec=args.cutout_size);return
    if args.command=='psf':
        from .stis.psf.cli import run
        run(args.config,args.stage);return
    if args.command=='report':
        from pathlib import Path
        from .common.target_runner import target_lock
        from .uvis.reporting import generate_report
        directory=Path(args.directory).expanduser().resolve()
        if not (directory/'manifest.json').is_file():parser.error('report requires a target UVIS directory containing manifest.json')
        state=json.loads((directory/'manifest.json').read_text())
        from .common.paths import ReductionPaths
        try:
            ReductionPaths.from_directory(directory,'report',target=state['target'])
            if state.get('instrument','uvis')!='uvis':raise ValueError('Not a UVIS reduction')
        except (KeyError,ValueError) as exc:parser.error(str(exc))
        with target_lock(directory):result=generate_report(directory)
        print(result)
        return
    if args.command=='review':
        from .common.uvis_review import review_uvis
        result=review_uvis(args.directory,args.filter,approve=args.approve,reviewer=args.reviewer,reason=args.reason)
        print(f"{result['target']} / {result['filter']} — {result['status']}")
        print('Report:',result['report']);print('Review image:',result['image'])
        for source in result['sources']:
            print('Source',source['source_id'],':',', '.join(source['failures']))
        if result['approved']:print('Published:',result['product'])
        else:print('Preview only. Inspect the current report; use --approve --reviewer NAME --reason TEXT to accept all listed field-source flux flags. Other scientific gates still apply.')
        return
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
                         'groups':[dict(target=k[0],element=k[1],count=len(v)) for k,v in groups(rows,args.instrument,field=getattr(args,'field',None)).items()]},indent=2));return
    if args.config:
        c=ReductionConfig.read(args.config)
        if c.instrument!=args.instrument:parser.error('Config instrument differs from command')
        for key in ('input','output','target','element','setup'):
            value=getattr(args,key)
            if value is not None:setattr(c,key,value)
        if getattr(args,'field',None) is not None:c.field=args.field
        c.force=c.force or args.force;c.resume=c.resume or args.resume;c.__post_init__()
    else:
        if not args.input or not args.output:parser.error('reduce requires --config or --input and --output')
        c=ReductionConfig(args.input,args.output,args.instrument,args.target,args.element,args.setup,resume=args.resume,force=args.force,field=getattr(args,'field',None))
    if args.instrument=='stis':
        from .stis import reduce_stis as reduce
    else:
        from .uvis import reduce_uvis as reduce
    result=reduce(c,dry_run=args.dry_run);print(json.dumps(result,default=str,indent=2))
