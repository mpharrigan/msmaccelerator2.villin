import os
import sys
import subprocess
from jinja2 import Environment
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--name',    help='name of the job', type=int, required=True)
parser.add_argument('--n_gpus',  help='number of gpus', type=int, default=6)
parser.add_argument('--repeats', help='number of  repeats', type=int, default=10)
parser.add_argument('--dry_run', help='dont actually submit', type=bool, default=False)
parser.add_argument('--time',    help='time', type=str, default='99:59:59')

args = parser.parse_args()

if args.n_gpus <= 0:
    parser.exit('n_gpus must be positive')
if args.repeats <= 0:
    parser.exit('n_gpus must be positive')


template = """#!/bin/bash
#MSUB -N msmaccelerator-{{name}}
#MSUB -l nodes=1:ppn={{n_gpus}}:gpus={{n_gpus}}
#MSUB -l walltime={{time}}
#MSUB -q longq
#MSUB -d /home/rmcgibbo/projects/msmaccelerator2.villin/project/

cd /home/rmcgibbo/projects/msmaccelerator2.villin/project/
mkdir -p logs
nvidia-smi &>logs/msmaccelerator-{{name}}.out

for i in `seq 1 {{repeats}}`
do
  {% for j in range(n_gpus) %}
  accelerator simulate --device_index={{ j }} &>logs/msmaccelerator-{{ name }}-$i-{{ j }}.out &
  sleep 10
  {% endfor %}
  
  wait
  accelerator model >>logs/msmaccelerator-{{ name }}-{{ i }}.out 2>&1
done
"""
env = Environment()
script = env.from_string(template).render(args.__dict__)

fn = 'msmaccelerator-{name}.pbs'.format(name=args.name)
with open(fn, 'w') as f:
    print script
    f.write(script)

cmd = ['/usr/local/bin/msub', os.path.abspath(fn)]
print ' '.join(cmd)

if not args.dry_run:
    print 'SUBMITTING JOB'
    code = subprocess.check_output(cmd)
    print 'JOB_ID: %s' % code
