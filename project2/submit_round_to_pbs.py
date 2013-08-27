import os
import sys
import subprocess
from jinja2 import Environment
from argparse import ArgumentParser


template = """#!/bin/bash
#MSUB -N msma-{{name}}
#MSUB -l nodes={{n_nodes}}:ppn=12:gpus={{n_gpus}}
#MSUB -l walltime={{time}}
#MSUB -q {{queue}}
#MSUB -d /home/harrigan/projects/msmaccelerator2.villin/project/

cd /home/harrigan/projects/msmaccelerator2.villin/project2/
mkdir -p logs
nvidia-smi &>logs/msmaccelerator-{{name}}.out

for i in `seq 1 {{repeats}}`
do
  {% for j in range(n_gpus*n_nodes) %}
  accelerator OpenMM --device_index={{ j }} &>logs/msmaccelerator-{{ name }}-$i-{{ j }}.out &
  sleep 10
  {% endfor %}

  wait
  accelerator model >>logs/msmaccelerator-{{ name }}-$i.out 2>&1
done
"""

queues = {
          'shortq': '13:00:00',
          'standard': '23:00:00',
          'rlongq': '440:00:00',
          'longq': '240:00:00'
          }

def submit_round_to_pbs(args):
    env = Environment()
    script = env.from_string(template).render(args.__dict__)

    fn = '{dir}/msmaccelerator-{name}.pbs'.format(name=args.name, dir=args.dir)
    with open(fn, 'w') as f:
        print script
        f.write(script)

    cmd = ['msub', os.path.abspath(fn)]
    print ' '.join(cmd)

    if not args.dry_run:
        print 'SUBMITTING JOB'
        code = subprocess.check_output(cmd)
        print 'JOB_ID: %s' % code

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--name', help='name of the job', type=str, required=True)
    parser.add_argument('--n_gpus', help='number of gpus', type=int, default=7)
    parser.add_argument('--n_nodes', help='number of nodes', type=int, default=7)
    parser.add_argument('--repeats', help='number of  repeats', type=int, default=10)
    parser.add_argument('--dry_run', help='dont actually submit', type=bool, default=False)
    parser.add_argument('--dir', help='Directory to store the job file', default='jobs')
    parser.add_argument('--queue', help='The queue to submit to', default='longq')

    args = parser.parse_args()

    if args.n_gpus <= 0:
        parser.exit('n_gpus must be positive')
    if args.repeats <= 0:
        parser.exit('repeats must be positive')
    if args.queue not in queues.keys():
        parser.exit('please supply a valid queue')

    args.time = queues[args.queue]

    submit_round_to_pbs(args)
