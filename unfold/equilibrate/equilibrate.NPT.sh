#PBS -n equilibration_NPT
#PBS -l nodes=1:gpus=6
#PBS -l walltime=5:00:00

source /home/rmcgibbo/.bashrc
cd /home/rmcgibbo/projects/msmaccelerator2.villin/unfold/equilibrate
python -u equilibrate.py > log.equil 2> err.equil
