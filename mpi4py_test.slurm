#!/bin/bash
#SBATCH --job-name=mpi4py_test
#SBATCH --partition=workq
#SBATCH --time=24:00:00
#SBATCH --account=k1068
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --array=1

HOME=$PWD
export PATH=/lustre/project/k1068/hjiang/software/anaconda2/bin:$PATH
module swap  PrgEnv-cray PrgEnv-gnu
#module load openmpi/1.10.6

#srun -n $SLURM_NTASKS --mpi=pmi2 python mpi4py_test.py
mpirun -n $SLURM_NTASKS python mpi4py_test.py
#srun --ntasks=1 --hint=nomultithread --ntasks-per-node=1 --ntasks-per-socket=1 --ntasks-per-core=1 --mem_bind=v,local  ${GROMACS_HOME}/bin/gmx_mpi grompp -c /lustre/project/k1068/118/v1.5/node000/node.pdb -f /lustre/project/k1068/118/v1.5/node000/run.mdp -pp /lustre/project/k1068/118/v1.5/node000/processed.top -o /lustre/project/k1068/118/v1.5/node000/run.tpr -p /lustre/project/k1068/118/v1.5/pars/topol.top -maxwarn 50
