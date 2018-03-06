#!/bin/sh
#
# Run a gromacs MD job on the local machine -- THIS IS JUST AN EXAMPLE
# This script is very simplistic and assumes that you can simply run
# gromacs with mpiexec (see below)

NCORES=4
#for Gromacs 5.0.4, Stephen added 20170927
#GMXDIR=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/bin
#export PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/bin:$PATH
#export LD_LIBRARY_PATH=/lustre/project/k1068/hjiang/software/gromacs_5.0.4_craympi_plumed/lib64:$LD_LIBRARY_PATH
export GROMACS_USE=_sp_plumed
module load gromacs/5.0.5
module swap PrgEnv-cray PrgEnv-gnu
export PATH=/lustre/project/k1068/hjiang/software/anaconda2/bin:$PATH
export PYTHON_EGGS_CACHE=/lustre/project/k1068/hjiang/software/anaconda2_eggs
GMXDIR=$GROMACS_HOME
# deffnm line is possibly modified by gromacs.setup
# (leave it as it is in the template)
DEFFNM=run

TPR=${DEFFNM}.tpr
OUTPUT=${DEFFNM}.out
PDB=${DEFFNM}.pdb

MDRUN_OPTS=""

$GMXDIR/gmx_mpi mdrun -ntomp $NCORES -nice 19 -v -deffnm ${DEFFNM} -c ${PDB} -cpi -append \
    $MDRUN_OPTS >$OUTPUT 2>&1

