#! /bin/bash

for dire in `ls -d tol020*`
do
  cd $dire/paths
    cp ../../pars/contour.dat .
    cp ../../pars/tors.plu .
    for i in {000..002}
    do
      plumed driver --mf_xtc iter$i.xtc --plumed tors.plu
      sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > cvs_iter$i.cv
      rm -rf COLVAR_TORS
    done
  cd ../../
done
