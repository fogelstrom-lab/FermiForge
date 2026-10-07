#!/usr/bin/env python3
"""Prepare (default) or execute bounded, fixed-state MPI strong-scaling tests.

No compilation, scheduler submission, plotting, or production relaxation.
Uses Python standard library only. Run inside a granted compute allocation.
"""
import argparse
import csv
import datetime as dt
import hashlib
import itertools
import json
import math
import os
import platform
from pathlib import Path
import re
import shlex
import signal
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def set_value(text, key, value):
    pattern = rf"(?im)^([ \t]*{re.escape(key)}[ \t]*=[ \t]*)(?:'(?:[^'\n]|'')*'|[^,!\n/]+)"
    text, count = re.subn(pattern, lambda m: m[1]+value, text, count=1)
    if count:
        return text
    text, count = re.subn(r'(?m)^[ \t]*/[ \t]*$', lambda m: f'  {key} = {value}\n/', text, count=1)
    if count != 1:
        raise ValueError('expected a namelist ending on its own / line')
    return text.rstrip()+'\n'

def quote(value):
    return "'"+str(value).replace("'", "''")+"'"

def stage_inputs(text, folder, project):
    """Archive external inputs once; relative solver paths are project-relative."""
    records = []
    for key in ('order_parameter_file','current_field_file','gauss_file','ozaki_file','restart_field_file'):
        if key == 'ozaki_file' and re.search(r"""(?im)^\s*ozaki_mode\s*=\s*['"]generate['"]""", text):
            continue
        match = re.search(rf"(?im)^\s*{key}\s*=\s*'((?:[^']|'')*)'", text)
        if not match or not match[1]:
            if key != 'restart_field_file':
                raise ValueError(f'explicit quoted {key} required for reproducibility')
            continue
        path = Path(match[1].replace("''", "'"))
        path = path if path.is_absolute() else project/path
        target = folder/key
        target.write_bytes(path.read_bytes())
        records.append(dict(key=key, source=str(path.resolve()), sha256=digest(target)))
        text = set_value(text,key,quote(target))
    return text, records

def compare_fields(reference, candidate, atol, rtol):
    def rows(path):
        with Path(path).open() as f:
            for line in f:
                if line.strip() and not line.lstrip().startswith('#'):
                    yield [float(s.replace('D','E')) for s in line.split()]
    largest = 0.0
    count = 0
    for a,b in itertools.zip_longest(rows(reference),rows(candidate)):
        count += 1
        if a is None or b is None or len(a)!=len(b):
            raise ValueError('mapped-field shapes differ between ranks/repetitions')
        for x,y in zip(a,b):
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError('nonfinite mapped field')
            largest=max(largest,abs(x-y))
            if abs(x-y)>atol+rtol*abs(x):
                raise ValueError(f'mapped field disagreement: {x} versus {y}')
    if count == 0:
        raise ValueError('empty mapped field')
    return largest

def launch(command, folder, environment, timeout, cwd):
    start=time.monotonic()
    with (folder/'run.log').open('w') as out, (folder/'mpi.log').open('w') as err:
        process=subprocess.Popen(command,cwd=cwd,env=environment,stdout=out,stderr=err,start_new_session=True)
        try:
            code=process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            # Only the process group created by this invocation is signalled.
            os.killpg(process.pid,signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL)
                process.wait()
            raise RuntimeError(f'launch timed out; see {folder}')
    if code:
        raise RuntimeError(f'launcher exited {code}; see {folder}/mpi.log')
    return time.monotonic()-start

def summarize(rows):
    grouped={}
    for row in rows:
        if not row['warmup']:
            grouped.setdefault(row['ranks'],[]).append(row)
    result=[]
    for p,values in sorted(grouped.items()):
        result.append(dict(ranks=p, repeats=len(values),
                           map_seconds_median=statistics.median(v['map_seconds'] for v in values),
                           map_seconds_min=min(v['map_seconds'] for v in values),
                           map_seconds_max=max(v['map_seconds'] for v in values),
                           elapsed_seconds_median=statistics.median(v['elapsed_seconds'] for v in values)))
    if result:
        base=result[0]
        for r in result:
            r['map_speedup_vs_baseline']=base['map_seconds_median']/r['map_seconds_median']
            r['relative_map_efficiency']=r['map_speedup_vs_baseline']*base['ranks']/r['ranks']
            r['elapsed_speedup_vs_baseline']=base['elapsed_seconds_median']/r['elapsed_seconds_median']
            r['relative_elapsed_efficiency']=r['elapsed_speedup_vs_baseline']*base['ranks']/r['ranks']
            r['rank_hours_per_launch']=r['ranks']*r['elapsed_seconds_median']/3600
    return result

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'examples/hpc_scaling_smoke.nml')
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--project-dir',type=Path,default=ROOT)
    parser.add_argument('--output-root',type=Path,default=ROOT/'work/hpc-scaling')
    parser.add_argument('--ranks',type=int,nargs='+',default=[1,2,4])
    parser.add_argument('--warmups',type=int,default=1)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--launcher',default='mpiexec -n {ranks}',
                        help='argv template, e.g. "srun --exclusive --ntasks={ranks} --cpus-per-task=1"')
    parser.add_argument('--timeout',type=float,default=1800,help='seconds per launch')
    parser.add_argument('--budget-seconds',type=float,default=7200,help='total elapsed campaign limit')
    parser.add_argument('--atol',type=float,default=1e-12)
    parser.add_argument('--rtol',type=float,default=1e-10)
    parser.add_argument('--execute',action='store_true',help='explicitly run; otherwise prepare only')
    args=parser.parse_args(argv)
    if min(args.ranks)<1 or len(set(args.ranks))!=len(args.ranks) or args.warmups<0 or args.repeats<1:
        parser.error('positive unique ranks, nonnegative warmups and positive repeats required')
    if not all(math.isfinite(x) and x>0 for x in [args.timeout,args.budget_seconds]):
        parser.error('finite positive time limits required')
    if not all(math.isfinite(x) and x>=0 for x in [args.atol,args.rtol]):
        parser.error('finite nonnegative comparison tolerances required')
    if '{ranks}' not in args.launcher:
        parser.error('launcher must contain {ranks}')
    executable=args.executable.resolve()
    if not executable.is_file():
        parser.error(f'executable not found: {executable}')
    project=args.project_dir.resolve()
    stamp=dt.datetime.now().strftime('%y%m%d-%H%M%S-%f')
    root=args.output_root.resolve()/stamp
    root.mkdir(parents=True,exist_ok=False)
    archive=root/'inputs'; archive.mkdir()
    cache=executable.parent/'CMakeCache.txt'
    if cache.is_file():
        (archive/'CMakeCache.txt').write_bytes(cache.read_bytes())
    text,dependencies=stage_inputs(args.input.read_text(),archive,project)
    (archive/'source.nml').write_text(args.input.read_text())
    for key,value in dict(maximum_iterations='0',checkpoint_interval='0',
                          prepare_only='.false.',require_convergence='.false.').items():
        text=set_value(text,key,value)
    input_file=archive/'benchmark.nml'; input_file.write_text(text)
    environment=os.environ.copy()
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
        environment[key]='1'
    # Do not impose macOS shared-memory transports or PMIx options on a cluster.
    jobs=[]
    for p in sorted(args.ranks):
        for repeat in range(args.warmups+args.repeats):
            folder=root/f'p{p:05d}-r{repeat:03d}'; folder.mkdir()
            outputs=[folder/name for name in ['input_fields.dat','mapped_fields.dat','metrics.txt',
                                             'final_fields.dat','history.dat','checkpoint.dat']]
            command=[token.format(ranks=p) for token in shlex.split(args.launcher)]
            command += [str(executable),str(input_file),*map(str,outputs)]
            jobs.append(dict(ranks=p,repeat=repeat,warmup=repeat<args.warmups,folder=str(folder),command=command))
    def git(*arguments):
        r=subprocess.run(['git','-C',str(project),*arguments],capture_output=True,text=True)
        return r.stdout.strip() if r.returncode==0 else 'unavailable'
    manifest=dict(created=stamp,status='prepared',input_sha256=digest(input_file),
                  platform=platform.platform(),host=platform.node(),
                  build_cache_sha256=digest(cache) if cache.is_file() else None,
                  executable=str(executable),executable_sha256=digest(executable),dependencies=dependencies,
                  git_commit=git('rev-parse','HEAD'),git_status=git('status','--short'),
                  tracked_diff_sha256=hashlib.sha256(git('diff','HEAD').encode()).hexdigest(),
                  launcher=args.launcher,timeout=args.timeout,budget_seconds=args.budget_seconds,
                  atol=args.atol,rtol=args.rtol,
                  environment={k:v for k,v in environment.items() if k.startswith(('OMPI_MCA_','PMIX_MCA_','SLURM_'))
                               or k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','LOADEDMODULES')},jobs=jobs)
    def save():
        (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    print(f'Prepared {len(jobs)} single-map launches: {root}',flush=True)
    if not args.execute:
        print('Nothing executed. Add --execute inside a compute allocation to run a fresh campaign.')
        return 0
    start=time.monotonic(); rows=[]; reference=None
    try:
        for job in jobs:
            remaining=args.budget_seconds-(time.monotonic()-start)
            if remaining<=0:
                raise RuntimeError('campaign time budget exhausted')
            folder=Path(job['folder'])
            elapsed=launch(job['command'],folder,environment,min(args.timeout,remaining),project)
            metrics=dict(line.split('=',1) for line in (folder/'metrics.txt').read_text().splitlines() if '=' in line)
            if metrics.get('passed','').strip()!='T' or metrics.get('terminal_status','').strip()!='single_map':
                raise ValueError(f'invalid map status in {folder}')
            mapped=folder/'mapped_fields.dat'
            difference=compare_fields(reference or mapped,mapped,args.atol,args.rtol)
            reference=reference or mapped
            seconds=float(metrics['maximum_rank_elapsed_seconds'])
            if not math.isfinite(seconds) or seconds<=0:
                raise ValueError('nonpositive/nonfinite map timing')
            active=int(metrics['active_points'])
            rows.append(dict(ranks=job['ranks'],repeat=job['repeat'],warmup=job['warmup'],
                             active_points=active,map_seconds=seconds,elapsed_seconds=elapsed,
                             max_field_difference=difference))
            (root/'measurements.json').write_text(json.dumps(rows,indent=2)+'\n')
            print(f"{job['ranks']} ranks, repeat {job['repeat']}: map {seconds:.3f}s, elapsed {elapsed:.3f}s",flush=True)
        summary=summarize(rows)
        with (root/'summary.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
        manifest['status']='completed'
    except (ValueError,RuntimeError,OSError) as error:
        manifest.update(status='failed',error=str(error));save()
        print(f'FAILED: {error}',flush=True)
        return 1
    save()
    print(f'Results: {root}/summary.csv')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
