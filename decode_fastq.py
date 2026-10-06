"""Decode in a fresh process: only FASTQ and optional public settings are read."""
import argparse
import json
from pathlib import Path
import read_pipeline as rp


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r1',required=True,type=Path)
    parser.add_argument('--r2',type=Path)
    parser.add_argument('--settings',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    report=args.output.with_suffix(args.output.suffix+'.report.json')
    if args.output.exists() or report.exists():parser.error('Output exists. Choose a new path.')
    c=rp.load_settings(args.settings.read_bytes()) if args.settings else rp.Settings()
    result=rp.run(args.r1.read_bytes(),c,args.r2.read_bytes() if args.r2 else b'')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if result.get('data') is not None:args.output.write_bytes(result['data'])
    report.write_text(json.dumps({k:v for k,v in result.items() if k not in ('data','restored_data')},indent=2),encoding='utf-8')
    print(json.dumps({'status':result['status'],'stored_sha256':result.get('sha256'),'reads':result['read_stats']}))


if __name__=='__main__':main()
