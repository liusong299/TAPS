#!/bin/bash
#SBATCH --job-name=TAPs
#SBATCH --partition=workq
#SBATCH --time=24:00:00
#SBATCH --account=k1068
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --array=1

#HOME=$PWD
#GMXDIR="/lustre/project/k1068/hjiang/software/gromacs_5.0.4/bin"
#export LD_LIBRARY_PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4/lib/${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
#GMXDIR="/lustre/project/k1068/hjiang/software/gromacs_4.5.4/bin"

HOME=$PWD
#GMXDIR=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/bin
export PATH=/lustre/project/k1068/hjiang/software/anaconda2/bin:$PATH
export PATH=/opt/intel/impi/2017.3.196/bin64:$PATH
#export PATH=/lustre/project/k1068/hjiang/software/plumed_2.1.3_sl/bin:$PATH
export GROMACS_USE=_sp_plumed
module load gromacs/5.0.5
module swap PrgEnv-cray PrgEnv-intel
module swap PrgEnv-gnu PrgEnv-intel

module load plumed/2.3.1
export PYTHON_EGGS_CACHE=/lustre/project/k1068/hjiang/software/anaconda2_eggs
#export PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/bin:$PATH
#export LD_LIBRARY_PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/lib64:$LD_LIBRARY_PATH
#for concorde, Stephen added 20170927
export PATH=/lustre/project/k1068/hjiang/software/concorde/TSP:$PATH

#srun -n 2 --mpi=pmi2 -t 600 /lustre/project/k1068/hjiang/software/anaconda2/bin/python runTAPS.py

#mpirun -n $SLURM_NTASKS python runTAPS.py
/opt/intel/impi/2017.3.196/bin64/mpirun -np 32 -ppn 32 python runTAPS.py
