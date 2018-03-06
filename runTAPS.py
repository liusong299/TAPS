from __future__ import print_function, division

import ConfigParser
import os

cf=ConfigParser.ConfigParser()
config_path='./pars/taps.ini'
if os.path.exists(config_path):
    cf.read(config_path)
else:
    raise ValueError("Noparameter file %s" % (parFile))
print (cf.sections())
for item in cf.items('taps'):
    print (item[0], '=', item[1])

n_start = cf.getint('taps','n_start')
n_taps = cf.getint('taps','n_taps')
iter_start = cf.getint('taps','iter_start')
n_iter = cf.getint('taps','n_iter')
dirPars = cf.get('taps','dirPars')          # dirPar
parFile = cf.get('taps','parFile')        # parameters filename
topFile = cf.get('taps','topFile')      # topology filename
p0File = cf.get('taps','p0File')    # initial path file
alignFile = cf.get('taps','alignFile')       # atom index file for alignment
rmsFile = cf.get('taps','rmsFile')          # atom index file for rmsd computation
tapsName = cf.get('taps','tapsName')

print (n_start, n_taps, iter_start, n_iter)
print (dirPars, parFile, topFile, p0File, alignFile, rmsFile, tapsName)

# =======================================================================================================================
#                           import mpi4py for parallel computing
# =======================================================================================================================
import multiprocessing
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# ======================================================================================================================
#                                       digits formater for iterations: "3 --> 03"
# ======================================================================================================================
def digits(s1):
    s2 = "%.3d" % s1
    return s2

from TAPS import *
from Confs import Confs
import time
import errno
import copy
import mdtraj as md
import numpy as np
import shutil
from copy import deepcopy


# =========================================================================================================
#                                         multiple independent taps
# =========================================================================================================
#n_start = 0
#n_taps = 1

# =========================================================================================================
#                                       number of iterations per taps
# =========================================================================================================
#n_iter = 2
#iter_start = 1

# =========================================================================================================
#                                                input files
# =========================================================================================================
#dirPars = 'pars_2'            # dirPar
#parFile = 'taps.par'        # parameters filename
#topFile = 'system.pdb'     # topology filename
#p0File = 'tmdAB_k1e5_rc050.xtc'   # initial path file
#alignFile = 'align.ndx'      # atom index file for alignment
#rmsFile = 'rms.ndx'         # atom index file for rmsd computation

#rint ("###DEBUG### Size:", size, "Rank:", rank, "begining.")
for i in range(n_start,n_taps+n_start):
    #tapsName = 'metEn_tol070_500ps_' + str(i)
    tapsName = tapsName + str(i).zfill(2)

    #print("###DEBUG### Barrier before iteration")
    #comm.Barrier()

    #if rank == 0:
    #    if not os.path.exists(tapsName):
    #        os.makedirs(tapsName)
    
    if rank == 0 and not os.path.exists(tapsName):
        try:
            os.makedirs(tapsName)
        except OSError as error:
            if error.errno != errno.EEXIST:
                raise
    #else:
    #    time.sleep(5)
    comm.Barrier()
    #f_log = open(tapsName + '/' + tapsName + '_' + str(rank) + '.log', 'w+')
    print(tapsName, ":")
    #print(tapsName, ":", file=f_log)

    print("  data initialization")
    #print("  data initialization", file=f_log)
    
    if rank == 0:
        t0 = time.time()
        #print("###DEBUG### Size:", size, "Rank:", rank, "running TAPS")
        #print("###DEBUG###", dirPars, parFile, topFile, p0File, alignFile, rmsFile)
        #print("###DEBUG###", dirPars, parFile, topFile, p0File, alignFile, rmsFile, file=f_log)
        taps = TAPS(dirPars, parFile, topFile, p0File, alignFile, rmsFile)
        te = time.time()
        print("    time-cost: ", te - t0, ' sec')
        #print("    time-cost: ", te - t0, ' sec', file=f_log)
        pathList = []
        refPath = copy.copy(taps.refPath)
        pathList.append(refPath)
        dirEvol = tapsName + '/paths'
        if not os.path.exists(dirEvol):
            os.makedirs(dirEvol)
        refPath.pathName = 'iter' + digits(0)
        refPath.exportPath(dirEvol)
    else:
        taps = None
        refPath = None
    taps = comm.bcast(taps, root=0)
    refPath = comm.bcast(refPath, root=0)
    comm.Barrier()

    for j in range(iter_start, iter_start + n_iter):
        # ==================================================================================================
        #                                          iteration index
        # ==================================================================================================
        iter = 'iter' + digits(i)
        # ==================================================================================================
        #                                        one taps iteration
        # ==================================================================================================
        dirMeta = tapsName + '/sampling/' + iter
        print("  ", iter, ": Preparing MetaD")
        #print("  ", iter, ": Preparing MetaD", file=f_log)
        #dirRUNs = taps.meta_dirs(refPath, dirMeta)
        comm.Barrier()
        if rank == 0:
            dirRUNs = taps.meta_dirs(refPath, dirMeta)
            t0 = time.time()
            #print("###DEBUG### Size:", size, "Rank:", rank, "running taps.meta_setup")
            taps.meta_setup(refPath, dirMeta, dirRUNs)
            t1 = time.time()
            print('   timecost: ', t1 - t0, ' sec')
            #print('   timecost: ', t1 - t0, ' sec', file=f_log)
            print("  ", iter, ": Sampling MetaD")
            #print("  ", iter, ": Sampling MetaD", file=f_log)
        else:
            #time.sleep(240)
	    dirRUNs = None
        dirRUNs = comm.bcast(dirRUNs, root=0)
        #print("###DEBUG### Barrier before meta_sample")
        comm.Barrier()
        if rank != 0:
            time.sleep(10)
            #print(" I am rank", rank, ", I slept 10 secs before meta_sample")
        #print("###DEBUG### Size:", size, "Rank:", rank, "running taps.meta_sample")
        taps.meta_sample(dirMeta, dirRUNs)
        comm.Barrier()
        #print("###DEBUG### Barrier after meta_sample")
        #
        t0 = time.time()
        if rank == 0:
            print('   timecost: ',t0 - t1, ' sec')
            #print('   timecost: ',t0 - t1, ' sec', file=f_log)
            print("  ", iter, ": Finding median(z) conformations, update path")
            #print("  ", iter, ": Finding median(z) conformations, update path", file=f_log)
        #print("###DEBUG### Barrier before meta_analyze")
        comm.Barrier()
        p_meta = taps.meta_analyze(dirMeta, dirRUNs)
        #print("###DEBUG### Barrier after meta_analyze")
        comm.Barrier()
        t1 = time.time()

        if rank == 0:
            print('   timecost: ', t1 - t0, ' sec')
            #print('   timecost: ', t1 - t0, ' sec', file=f_log)
            p_meta.pathName = iter
            p_meta.exportPath(dirEvol)
            print(' ')
            #print(' ', file=f_log)
            refPath = deepcopy(p_meta)
        else:
            p_meta.pathName = None
            p_meta.exportPath = None
            refPath = None
        p_meta = comm.bcast(p_meta, root=0)
        refPath = comm.bcast(refPath, root=0)
        comm.Barrier()
'''