#!/bin/bash
#SBATCH --partition=workq
#SBATCH --account=k1068
#SBATCH --job-name="gromacs"
#SBATCH --nodes=1
#SBATCH --extra-node-info=2:16:1
#SBATCH --time=4:00:00
#SBATCH --exclusive
#----------------------------------------------------------#
export GROMACS_USE=_sp_plumed
#export GROMACS_USE=_dp_plumed
module load gromacs/5.0.4
module swap PrgEnv-cray PrgEnv-intel
export OMP_NUM_THREADS=1
#----------------------------------------------------------#
echo "The job "${SLURM_JOB_ID}" is running on "${SLURM_JOB_NODELIST}
#----------------------------------------------------------#
srun --ntasks=1 --hint=nomultithread --ntasks-per-node=64 --ntasks-per-socket=16 --ntasks-per-core=1 --mem_bind=v,local ${GROMACS_HOME}/bin/gmx_mpi grompp -c /lustre/project/k1068/118/v1.5/node000/node.pdb -f /lustre/project/k1068/118/v1.5/node000/run.mdp -pp /lustre/project/k1068/118/v1.5/node000/processed.top -o /lustre/project/k1068/118/v1.5/node000/run.tpr -p /lustre/project/k1068/118/v1.5/pars/topol.top -maxwarn 50
