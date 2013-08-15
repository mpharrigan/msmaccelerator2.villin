#PBS -n equilibration_NVT
#PBS -l nodes=1:gpus=6
#PBS -l walltime=22:59:00
#PBS -q default

source /home/rmcgibbo/.bashrc
cd /home/rmcgibbo/projects/msmaccelerator2.villin/unfold/equilibrate
python -u equilibrate.NVT.py > NVT.log
