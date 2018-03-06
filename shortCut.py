from __future__ import print_function, division
from TAPS import *
from Confs import Confs
import time
import mdtraj as md
import numpy as np
import shutil

dirPars = 'pars'            # dirPar
parFile = 'taps.par'        # parameters filename
topFile = 'system.pdb'     # topology filename
p0File = 'p0_bb_rc020.xtc'   # initial path file
alignFile = 'align.ndx'      # atom index file for alignment
rmsFile =  'rms.ndx'         # atom index file for rmsd computation

taps = TAPS(dirPars, parFile, topFile, p0File, alignFile, rmsFile)
p0 = taps.refPath
p1 = p0.rmClose(0.028)
p1.nodes.save('p2_rc028.xtc')
p1.pcv()

