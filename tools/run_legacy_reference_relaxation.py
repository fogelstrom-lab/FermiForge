#!/usr/bin/env python3
"""Full-domain legacy map, short relaxation, and independent final-state map."""
import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iterations', type=int, default=3)
    parser.add_argument('--ranks', type=int, default=10)
    parser.add_argument('--batch-directory', type=Path, required=True)
    args = parser.parse_args()
    if args.iterations < 1 or args.ranks < 1:
        parser.error('iterations and ranks must be positive')
    batch = args.batch_directory.resolve()
    batch.mkdir(parents=True, exist_ok=True)
    status_path = batch / 'batch_status.json'
    if status_path.exists():
        parser.error('batch directory already has a status file; choose a new directory')
    status = {'started': dt.datetime.now().isoformat(), 'status': 'preparing',
              'stages': [], 'iterations': args.iterations, 'ranks': args.ranks}

    def save():
        status_path.write_text(json.dumps(status, indent=2) + '\n')

    save()
    try:
        prepared = subprocess.run([
            sys.executable, str(ROOT / 'tools/run_converged_reference.py'),
            '--prepare-only', '--fine-spacing', '0.2', '--ranks', str(args.ranks),
        ], cwd=ROOT, check=True, text=True, capture_output=True)
        print(prepared.stdout, flush=True)
        input_path = Path(prepared.stdout.split('Prepared input: ', 1)[1].strip())
        status['input'] = str(input_path)
        restart = None
        for name, iterations in [('01-initial-map', 0),
                                 ('02-relaxation', args.iterations),
                                 ('03-final-map', 0)]:
            stage_root = batch / name
            stage_root.mkdir()
            status['status'] = 'running'
            status['current_stage'] = name
            save()
            command = [sys.executable, '-u', str(ROOT / 'tools/run_double_core_from_scratch.py'),
                       '--input', str(input_path), '--iterations', str(iterations),
                       '--ranks', str(args.ranks), '--run-root', str(stage_root),
                       '--run-name', name]
            if name != '01-initial-map':
                command.append('--no-build')
            if name == '03-final-map':
                command.extend(['--sparse-restart', str(restart)])
            print('Starting ' + name, flush=True)
            subprocess.run(command, cwd=ROOT, check=True)
            reports = list(stage_root.glob('*/scratch_run_report.json'))
            if len(reports) != 1:
                raise RuntimeError('Expected exactly one stage report')
            report = json.loads(reports[0].read_text())
            if not report['passed']:
                raise RuntimeError('Stage report failed: ' + str(reports[0]))
            status['stages'].append({'name': name, 'report': str(reports[0]),
                                     'result': report['result']})
            if name == '02-relaxation':
                restart = reports[0].parent / 'sparse_final_state.dat'
                if not restart.is_file():
                    raise RuntimeError('Missing final restart state')
            save()
        initial = status['stages'][0]['result']
        final = status['stages'][-1]['result']
        status['comparison'] = {
            'initial_relative_l2': initial['sample_relative_l2'],
            'final_relative_l2': final['sample_relative_l2'],
            'initial_separation': initial['measured_half_core_separation'],
            'final_separation': final['measured_half_core_separation'],
            'note': 'Execution success is not a physical convergence criterion.',
        }
        status['status'] = 'complete'
        print(json.dumps(status['comparison'], indent=2), flush=True)
    except Exception as error:
        status['status'] = 'failed'
        status['error'] = str(error)
        raise
    finally:
        save()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
