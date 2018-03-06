#! /bin/bash
cp ../../../pars/contour.dat .
cp ../../../pars/tors.plu .

plumed driver --mf_xtc mz.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp_tr.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp_tr.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp_in.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp_in.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp_in_rc.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp_in_rc.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp_in_st.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp_in_st.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp_in_st_tmd.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp_in_st_tmd.cv
rm -rf COLVAR_TORS

plumed driver --mf_xtc mz_tsp_in_st_tmd_rc.xtc --plumed ./tors.plu
sed /\#/d COLVAR_TORS | awk '{print $2,$3}' > mz_tsp_in_st_tmd_rc.cv
rm -rf COLVAR_TORS

