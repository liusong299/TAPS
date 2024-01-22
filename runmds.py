from __future__ import print_function, division
import mdtraj as md
import numpy as np
import os
import re
import shutil
import time
import errno
import random as rdm
from numpy import *
from Confs import Confs
from copy import deepcopy

#Gromacs
import gromacs.setup
import gromacs.run
import gromacs.tools
import gromacs

def digits(s1):
	s2 = "%.3d" % s1
	return s2

#MDS Analysis
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
 
from sklearn import manifold
from sklearn.metrics import euclidean_distances


n_iter = 72	#number of paths

dire = './pars'
pdire = './PathInd0/paths'
topfile = dire + '/' + 'protein.pdb'

alignName = dire + '/' + 'align.ndx'
f = open(alignName, 'r')
fr = f.readlines()
align = []
for i in range( len(fr) ):
	align.append( int(fr[i]) )
f.close()

print('align 0')

rmsName = dire + '/' + 'rms.ndx'
f = open(rmsName, 'r')
fr = f.readlines()
rms = []
for i in range( len(fr) ):
	rms.append( int(fr[i]) )
f.close()

print('rms 0')


trj_list = []
n_conftrj = np.ones(n_iter+1, dtype= int)
totConf = 0
for i in range(n_iter+1):
	pname = pdire + '/' + 'iter' + digits(i) + '.xtc'
	pnode = Confs.traj2conf(md.load(pname, top=topfile))
	trj_list.append( pnode )
	n_conftrj[i] = len(pnode)
	totConf += len(pnode)

print('totConf is %d' % totConf)

trj2one = Confs.merge(trj_list)
print(len(trj_list))

rmsdf = np.ones( (totConf, totConf), dtype=float)
for i in range (totConf):
	trj2one.superpose(trj2one, i, align)
	rmscal = md.rmsd( trj2one, trj2one, i, rms)
	for j in range(totConf):
		rmsdf[j][i] = rmscal[j]
		rmsdf[i][j] = rmscal[j]

print("Next: perform MDS Cal\n")
 
mds=manifold.MDS(n_components=2, max_iter=300, eps=0.00001, random_state=3, dissimilarity="precomputed", n_jobs=1)
rmsdf_mds=mds.fit(rmsdf).embedding_

print("finished MDS Cal\n")

num_re = 0
mdsName = './mdsdata'
os.mkdir(mdsName)
for i in range(n_iter+1):
	fwName = mdsName + '/' + 'iter' + digits(i) + '.dat'
	fw = open(fwName, 'w+')
	for j in range(num_re, num_re+n_conftrj[i]):
		line = " Coord1OfConf %f  Coord2OfConf %f " % (rmsdf_mds[j][0], rmsdf_mds[j][1])
		print(line, file= fw)
	fw.close()
	num_re += n_conftrj[i]

print("Analysis Finished")
