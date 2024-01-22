TAPS released in 2024 (Travelling-salesman based Automatic Path-Searching)
===
Requirements
==
  * GROMACS >=2019.4
  * Plumed >=2.7
  * Cuda >=10.2
  * Openmpi >=4.0.2
  * python >=3.7.4 
  * numpy >=1.17.2
  * MDTraj >=1.9.3
  * GromacsWrapper >=0.8.0
  * mpi4py >=3.0.3
  * [concorde](http://www.math.uwaterloo.ca/tsp/concorde.html) :
  * These requirements were tested and proved to work. Feel free to test

Install Python and Python Packages
==
We highly recommend that you download the Python 3.x version of Anaconda, which is a completely free enterprise-ready Python distribution for large-scale data processing, predictive analytics, and scientific computing.

Design and Usage
==
Travelling-salesman based Automatic Path-Searching (TAPS)
===
  * No static coordinate space (CVs): ordered high dimensional conformations
  * Perpendicular relaxation: Quickly find MFEP segments
  * Automatic re-order of path nodes by Travelling-salesman
  * Enhanced sampling along path on PCV-z by MetaD
  * Validated for three protein systems (76-303 residues and total 30000-80000 atoms, [TAPStest](https://pubs.acs.org/journal/jctcce))      
            
Tutorial
==
  1. The parameters used for TAPS can be modified in pars/taps.par file
  2. Serial Running:
```
     1. Change runMode to "serial" in taps.par
     2. > python runTAPS.py                                  
```
  3. Parallel Running:
```
     1. Change runMode to "Parallel" in taps.par
     2. > mpirun -np 8 python runTAPS.py          
```

TODO
=

[  ] Update Tutorial

[  ] Test more
 
Authors:
=

Lizhe Zhu: *zhulizhe@cuhk.edu.cn*

Song Liu: *sliubu@connect.ust.hk*

Kun Xi: *xikun@cuhk.edu.cn*

Contributors:
=
 Maybe you !
