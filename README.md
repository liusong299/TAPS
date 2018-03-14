# Travelling-salesman based Automatic Path-Searching (TAPS)
## Requirements
- python 2.7.8
- numpy 1.9.2
- scipy 0.15.1
- MDTraj 1.5.0
- scikit-learn 0.16.1
- matplotlib 1.4.3
- GromacsWrapper 0.5.2
- mpi4py 2.0.0
- [concorde](http://www.math.uwaterloo.ca/tsp/concorde.html) (Travelling-salesman solver, implemented in C++)

**OBS:**
- These requirements were tested and proved to work. Feel free to test

## Install Python and Python Packages
We highly recommend that you download the Python 2.7 version of Anaconda, which is a completely free enterprise-ready Python distribution for large-scale data processing, predictive analytics, and scientific computing.

## Design and Usage
### Travelling-salesman based Automatic Path-Searching (TAPS)
- No static coordinate space (CVs): ordered high dimensional conformations
- Perpendicular relaxation: Quickly find MFEP segments
- Automatic re-order of path nodes by Travelling-salesman
- Enhanced sampling along path by MetaD
- 5-10x speedup against string method (swarms-of-trajectories)

## Tutorial
1. The parameters of TAPs could be modified in pars/taps.ini file
2. Serial Running:

		1.Change runMode to "Serial" in taps.ini
		2. > python runTAPs.py

3. Parallel Running:

		1. Change runMode to "Parallel" in taps.ini
		2. > mpirun -np 4 python runTAPs.py
		
## TODO
- [ ] Update Tutorial
- [ ] Test more

## labels:
- [DON] - Done
- [IMP] - To be improoved
- [BUG] - Buggy and experimental

## Authors:
Lizhe Zhu: *lizhezhu11@gmail.com*

[Stephen Liu](/https://github.com/liusong299/): *liusong299@gmail.com*

## Contributors:
Maybe you !