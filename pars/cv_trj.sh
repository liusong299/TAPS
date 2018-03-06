#! /bin/bash
cp ../../../pars/contour.dat .
cp ../../../pars/tors.plu .
for dire in `ls -d node*`
do
  cd $dire
  plumed driver --mf_xtc run.xtc --plumed ../tors.plu
  sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > tors.cv
  rm -rf COLVAR_TORS
  cd ..
done
