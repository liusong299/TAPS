#! /bin/bash
make_ndx -f $1.pdb -o index<<eof
q
eof

trjconv -s $1.pdb -pbc whole -o whole.xtc -f $1.xtc -n index <<eof
1
1
eof

trjconv -s $1.pdb -fit rot+trans -f whole.xtc -o protein.xtc -n index <<eof
1
1
eof

trjconv -s $1.pdb -o protein.gro -f protein.xtc -dump 0 <<eof
1
eof

#trjconv -s protein.gro -f protein.xtc -skip 10 -o short.xtc <<eof
#0
#eof

rm *#
