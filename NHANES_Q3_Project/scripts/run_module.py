"""One-command entry, with a real-data fit gate and compatible-stage reuse."""
from runtime import limit_threads
limit_threads()
import argparse
import subprocess
import sys
from runtime import ROOT, OUT, read_json, event, write_json, now

def command(script):
    cmd=[sys.executable,str(ROOT/'scripts'/script)]
    event('command_start',command=cmd)
    subprocess.run(cmd,check=True)
    event('command_complete',script=script)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--report-only','--outputs-only',dest='report_only',action='store_true',
                        help='Render saved CSV/JSON results; never enter fitting, model loading or cache checks.')
    parser.add_argument('--report-dir',help='Destination for report-only artifacts (default: output/Patterns_Integrated).')
    args=parser.parse_args()
    if args.report_only:
        if args.prepare_only:
            parser.error('--prepare-only and --report-only are mutually exclusive')
        from saved_reports import generate
        generate(args.report_dir)
        return
    if args.report_dir:
        parser.error('--report-dir requires --report-only')
    if not read_json(ROOT/'protocol/PROTOCOL_RECORD.json').get('cadence_amendment_utc'):
        command('initialize_execution.py')
    if not (OUT/'checks/preflight.json').exists():
        command('preflight.py')
    if not (OUT/'checks/artificial_modeling.json').exists():
        command('check_modeling.py')
    if args.prepare_only:
        if read_json(ROOT/'protocol/PROTOCOL_RECORD.json').get('first_study_fit_utc'):
            print('Study fitting already occurred; retain actual logs and completed drafts. Preparation does not reset history.')
            return
        if (ROOT/'scripts/write_manuscript.py').exists():
            from write_manuscript import prepare
            prepare()
        write_json(OUT/'logs/preparation_complete.json',dict(utc=now(),status='prepared_pending_remote_repository',study_fits=0))
        print('Preparation complete. No study model has been fitted. Supply and verify the designated remote protocol record before continuing.')
        return
    command('run_benchmark.py')
    # Always evaluate the predictions used by this invocation. A saved summary
    # or a previous rendering marker does not establish correspondence.
    command('evaluate_fixed_predictions.py')
    command('make_outputs.py')
    print('Analysis completed. Inspect the generated tables and Figure 6 in outputs.')

if __name__=='__main__':
    main()
