#!/bin/bash
#SBATCH --job-name=TAPs_test_plumed
#SBATCH --partition=workq
#SBATCH --time=24:00:00
#SBATCH --account=k1068
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --array=1

#HOME=$PWD
#GMXDIR="/lustre/project/k1068/hjiang/software/gromacs_5.0.4/bin"
#export LD_LIBRARY_PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4/lib/${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
#GMXDIR="/lustre/project/k1068/hjiang/software/gromacs_4.5.4/bin"

HOME=$PWD
#GMXDIR=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/bin
export PATH=/lustre/project/k1068/hjiang/software/anaconda2/bin:$PATH
#export PATH=/lustre/project/k1068/hjiang/software/plumed_2.1.3_sl/bin:$PATH
#export GROMACS_USE=_dp_plumed
#module load gromacs/5.0.5
#module swap PrgEnv-cray PrgEnv-intel

module load plumed/2.3.1
export PYTHON_EGGS_CACHE=/lustre/project/k1068/hjiang/software/anaconda2_eggs
#export PATH=/lustre/sw/xc40/plumed/2.1.2/cnl5.2_intel15.0.2/install/bin:$PATH
#export PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/bin:$PATH
#export LD_LIBRARY_PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/lib64:$LD_LIBRARY_PATH
#for concorde, Stephen added 20170927
export PATH=/lustre/project/k1068/hjiang/software/concorde/TSP:$PATH

#srun --ntasks=1 --hint=nomultithread --ntasks-per-node=1 --ntasks-per-socket=1 --ntasks-per-core=1 --mem_bind=v,local  ${GROMACS_HOME}/bin/gmx_mpi grompp -c /lustre/project/k1068/118/v1.5/node000/node.pdb -f /lustre/project/k1068/118/v1.5/node000/run.mdp -pp /lustre/project/k1068/118/v1.5/node000/processed.top -o /lustre/project/k1068/118/v1.5/node000/run.tpr -p /lustre/project/k1068/118/v1.5/pars/topol.top -maxwarn 50

srun --ntasks=1 --hint=nomultithread --ntasks-per-node=1 --ntasks-per-socket=1 --ntasks-per-core=1 --mem_bind=v,local plumed driver --mf_xtc /lustre/project/k1068/118/v1.5/node000/tmp.xtc --plumed /lustre/project/k1068/118/v1.5/node000/iter000.plu
