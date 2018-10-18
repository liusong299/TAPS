# =======================================================================================================================
#                        global import: future_print, numpy, mdtraj, os, re, argparse
# =======================================================================================================================
from __future__ import print_function, division
import numpy as np
import mdtraj as md
import os
import re
import shutil
import time
import errno
# =======================================================================================================================
#                           import python wrappers for MD engines (GROMACS, NAMD, AMBER)
# =======================================================================================================================
import gromacs.setup
import gromacs.run
import gromacs.tools
import gromacs

# For treating gromacs warnings at mdrun as exception in python
#   such that mdrun can be terminated if
import warnings

warnings.simplefilter('error', gromacs.AutoCorrectionWarning)
warnings.simplefilter('error', gromacs.BadParameterWarning)
warnings.simplefilter('error', gromacs.GromacsValueWarning)
warnings.simplefilter('error', gromacs.GromacsFailureWarning)

from Confs import Confs
from PcvInd import PcvInd
from Path import Path


# =======================================================================================================================
#                           import mpi4py for parallel computing
# =======================================================================================================================
import multiprocessing
from mpi4py import MPI
comm = MPI.COMM_WORLD
#print("hello world")
size = comm.Get_size()
rank = comm.Get_rank()
threads = 1
# ======================================================================================================================
#                                       digits formater for iterations: "3 --> 03"
# ======================================================================================================================
def digits(s1):
    s2 = "%.3d" % s1
    return s2

# ==================================================================================================================
#                                       wrapper for a single run of metaD
#                     engine-specific implementation of sampling is realized in this function
# ==================================================================================================================
def runMeta(dire, engine, runName, pluName, cvOut, trjName):
    #lock is for parallel process
    if engine == 'GROMACS':
        meta = gromacs.run.MDrunner(dire, ntomp=threads, deffnm=runName, plumed=pluName)
        meta.run()
        # make sure xtc is complete, under restrains, high energy confs may be generated
        # use trjconv to disregard "unphysical (CV values is nan)" frames
        # meanwhile, there are two other cases to consider
        #  1. when there is no output or only one line in CV file (sampling crashed in the first step)
        #       this can be dealt with by put the endtime at 0
        #  2. when the last line is incomplete
        #       just remove the last line
        keepLineIndex = []
        input = open(dire + '/' + cvOut, 'r+')
        for i, line in enumerate(input):
            if line[0] != "#":
                if not re.match('.*nan.*', line):
                    keepLineIndex.append(i)
        lineCount = len(keepLineIndex)
        output = open(dire + '/colvar_filter', 'w+')
        if lineCount == 1:  # only one line: keep this line
            endTime = 0
            lines = input.readlines()
            line = lines[keepLineIndex[0]]
            output.write(line)
        else:  # many lines in CV file, only remove the last line (incomplete when sampling crashes)
            input = open(dire + '/' + cvOut, 'r+')
            lines = input.readlines()
            output = open(dire + '/colvar_filter', 'w+')
            for k in range(lineCount - 1):  # remove the last line (incomplete when sampling crashes)
                line = lines[keepLineIndex[k]]
                endTime = line.split()[0]
                output.write(line)
        input.close()
        output.close()
        shutil.move(dire + '/' + trjName, dire + '/bak_' + trjName)
        trjconv = gromacs.tools.Trjconv(s=dire + '/' + runName + '.tpr', f=dire + '/bak_' + trjName,
                                        o=dire + '/' + trjName, e=endTime, \
                                        ur='compact', center=True, pbc='mol', input=('Protein', 'System'))
        trjconv.run()
        os.remove(dire + '/bak_' + trjName)
    else:
        raise ValueError("MD engines other than GROMACS are not support yet")


# ==================================================================================================================
#                                       wrapper for a single run of metaD
#                     engine-specific implementation of sampling is realized in this function
# ==================================================================================================================
def runTMD(dire, engine, runName, pluName, trjName):
    if engine == 'GROMACS':
        tmd = gromacs.run.MDrunner(dire, ntomp=threads, deffnm=runName, plumed=pluName)
        tmd.run()
        shutil.move(dire + '/' + trjName, dire + '/bak_' + trjName)
        trjconv = gromacs.tools.Trjconv(s=dire + '/' + runName + '.tpr', f=dire + '/bak_' + trjName, \
                                            o=dire + '/' + trjName, ur='compact', center=True, pbc='mol', \
                                        input=('Protein', 'System'))
        trjconv.run()
        os.remove(dire + '/bak_' + trjName)
    else:
        raise ValueError("MD engines other than GROMACS are not support yet")
# ======================================================================================================================
#                               class TAPS: encoding methods for each iteration of TAPS
# ======================================================================================================================
class TAPS(object):
    # default structure file, these names are important during sampling and plumed computation
    nodeName = 'node.pdb'
    runName = 'run'
    trjName = 'run.xtc'
    trjFilter = 'filtered.xtc'
    pluName = 'plumed.dat'
    cvOut = 'COLVAR'

    # ==================================================================================================================
    #    constructor: read in taps parameters and relevant files (system topology, initial path, PCV definition)
    # ==================================================================================================================
    def __init__(self, dire='pars', parFile='taps.par', topFile='protein.pdb', p0='path0.xtc', alignFile='align.ndx', \
                 rmsFile='rms.ndx'):

        # check if inputs exists
        if not os.path.isdir(dire):
            raise ValueError("Directory %s for initial path & parameters does not exist" % dire)
        if not os.path.exists(dire + '/' + parFile):
            raise ValueError("Parameters file %s is not found in directory %s" % (parFile, dire))
        if not os.path.exists(dire + '/' + topFile):
            raise ValueError("Structure file %s is not found in directory %s" % (topFile, dire))
        if not os.path.exists(dire + '/' + p0):
            raise ValueError("Trajectory of initial path (%s) is not found in directory '%s'" % (p0, dire))
        if not os.path.exists(dire + '/' + alignFile):
            raise ValueError("Atom index file for alignment (%s) is not found in directory %s" % (alignFile, dire))
        if not os.path.exists(dire + '/' + rmsFile):
            raise ValueError("Atom index file for rms computation (%s) is not found in directory %s" % (rmsFile, dire))

        # record root directory
        self.dirRoot = os.getcwd()

        # record directory for initial path and parameters
        self.dirPar = self.dirRoot + '/' + dire

        # record topology file name and position
        self.topNAME = topFile
        self.topFile = self.dirPar + '/' + topFile

        # record alignment index file position
        self.alignFile = self.dirPar + '/' + alignFile

        # record rms index file position
        self.rmsFile = self.dirPar + '/' + rmsFile

        # load atom indices for PCV definition (alignment & rmsd calculation)
        align = np.loadtxt(self.dirPar + '/' + alignFile, dtype=np.int32)
        rms = np.loadtxt(self.dirPar + '/' + rmsFile, dtype=np.int32)
        self.pcvInd = PcvInd(align, rms)

        # load initial refPath (compute initial s, included)
        self.refPath = Path('iter' + digits(0), self.pcvInd)
        self.refPath.loadFromTRJ(self.dirPar + '/' + p0, self.dirPar + '/' + topFile)

        # initialize initial node (extracting from initial path)
        self.initNode = self.refPath.nodes.slice(0)

        # initialize final node (extracting from initial path)
        self.finalNode = self.refPath.nodes.slice(self.refPath.n_nodes - 1)

        # read in parameters for MD and metaD
        fr = open(self.dirPar + '/' + parFile, 'r+')
        pars = fr.read()
        fr.close()

        # MD parameters
        # engine specific input check
        match = re.search("engine=.*\n", pars)
        if match is not None:
            self.engine = re.split('=', match.group(0).rstrip('\n'))[1]
        else:
            raise ValueError("MD engine not given in parameter file %s" % (parFile))
        if self.engine == 'GROMACS':
            match = re.search("groTOP=.*\n", pars)
            if match is not None:
                self.groTOP = re.split('=', match.group(0).rstrip('\n'))[1]
                if not os.path.exists(self.dirPar + '/' + self.groTOP):
                    raise ValueError("GROMACS topology file %s is not found in directory %s" % (self.groTOP, \
                                                                                                self.dirPar))
            else:
                raise ValueError("GROMACS topology file not given in %s" % (parFile))
            match = re.search("groMDP=.*\n", pars)
            if match is not None:
                self.groMDP = re.split('=', match.group(0).rstrip('\n'))[1]
                if not os.path.exists(self.dirPar + '/' + self.groMDP):
                    raise ValueError("GROMACS template mdp file %s is not found in directory %s" % (self.groMDP, \
                                                                                                    self.dirPar))
            else:
                raise ValueError("gromacs mdp file %s not given in %s" % (parFile))
        elif self.engine == 'NAMD':
            raise ValueError('NAMD is not supported yet')
        elif self.engine == 'AMBER':
            raise ValueError('AMBER is not supported yet')
        else:
            raise ValueError("unknown MD engine %s" % self.engine)

        # mode = {serial, parallel, qjob}
        match = re.search("runMode=.*\n", pars)
        if match is not None:
            self.mode = re.split('=', match.group(0).rstrip('\n'))[1]
        else:
            raise ValueError("Mode of running (runMode) not given in parameter file %f" % (parFile))

        # time step
        match = re.search("timeStep=.*\n", pars)
        if match is not None:
            self.timeStep = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("MD timestep (timestep, unit: ps) not given in parameter file %f" % (parFile))

        match = re.search("lenSample=.*\n", pars)
        if match is not None:
            self.lenSample = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Amount of sampling per taps iteration ('lenSample', unit: ps) not given in \
            parameter file %f" % (parFile))

        # gaussian height for MetaD
        match = re.search("gauHeight=.*\n", pars)
        if match is not None:
            self.gauHeight = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Height of gaussian for MetaDynamics (gh) not give in parameter file %s" % (parFile))

        # gaussian width for MetaD
        match = re.search("gauWidth=.*\n", pars)
        if match is not None:
            self.sigma = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Width of gaussian for MetaDynamics (sigma) not given in parameter file %s" % (parFile))

        # deposition interval for MetaD
        match = re.search("tauMetaD=.*\n", pars)
        if match is not None:
            self.tauMetaD = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Period for adding gaussians (tauMetaD) not given in parameter file %s" % (parFile))

        # biasFactor for well-tempered MetaD
        match = re.search("biasFactor=.*\n", pars)
        if match is not None:
            self.bf = int(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("BiasFactor for wt-MetaDynamics (biasFactor, 2-10) not given  in parameter file %s" % (parFile))

        # system temperature for well-tempered MetaD
        match = re.search("temp=.*\n", pars)
        if match is not None:
            self.temp = int(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Temperature for wt-MetaDynamics (temp) not given in parameter file %s" % (parFile))

        # output frequency of trajectories
        match = re.search("freqTRJ=.*\n", pars)
        if match is not None:
            self.freqTRJ = int(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Output frequency of sampling trajectories (freqTRJ) not given in parameter file %f" % (parFile))

        fr = open(self.dirPar + '/' + self.groMDP, 'r+')
        linesMDP = fr.readlines()
        fr.close()
        mdpFile = 'md.mdp'
        fw = open(self.dirPar + '/' + mdpFile, 'w+')
        fw.writelines(linesMDP)
        print('nstxout-compressed= %d' % self.freqTRJ, file=fw)
        fw.close()
        self.groMDP = mdpFile

        # output frequency of trajectories
        match = re.search("kappa=.*\n", pars)
        if match is not None:
            self.kappa = int(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Wall strength on PCV-s (kappa, 10-50) not given in parameter file %f" % (parFile))

        # tolerable restraining potential to ensure "physically irrelevant" conformations are selected
        # selecting frames with small restrain potential is a more direct approach than ds-s[0]<sTol
        # because it makes the choice independent from the kappa of the restraining potential
        match = re.search("tolRS=.*\n", pars)
        if match is not None:
            self.rsTol = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Tolerable restraining potential (rsTol) not found in parameter file %s \n  This parameter\
             is crucial for selecting frames from MetaD trajectories" % (parFile))

	    # parameters for path-reparameterization
	    # tolerable distance between neighbor nodes, used for reparameterization
        match = re.search("tolDist=.*\n", pars)
        if match is not None:
            self.tolDist = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            raise ValueError("Tolerable maximum distance (tolDist) between neighbor nodes not given in parameter\
             file %s\n  This parameter is crucial for path reparameterzation" % (parFile))

        # tolerable asymmetry factor, determines how much deviation from the used for path reparameterization
        match = re.search("devMID=.*\n", pars)
        if match is not None:
            self.devMID = float(re.split('=', match.group(0).rstrip('\n'))[1])
            if self.devMID > 1 or self.devMID <= 0:
                raise ValueError("Parameter devMID out of range ( 0<devMID<=1 required )")
        else:
            raise ValueError(
                "Tolerable deviation from vertical line between two distant nodes (devMID) is not given in parameter \
                file %s\n  This parameter is crucial for path reparameterzation" % (parFile))

        # tolerable cosTheta, used for reparameterization
        match = re.search("tolCos=.*\n", pars)
        if match is not None:
            self.tolCos = float(re.split('=', match.group(0).rstrip('\n'))[1])
            if self.tolCos > 0.5:
                self.tolCos = 0.5
                print("Tolerable cos(theta) in parameter file %s must be <=0.5\n  setting to 0.5" % (parFile))
        else:
            raise ValueError(
                "Tolerable cos(theta) to select \"middle\" conformations between neighbor nodes is not given in \
                parameter file %s" % (parFile))

        # straightening factor
        sub_i = self.initNode.atom_slice(self.pcvInd.atomSlice)
        sub_f = self.finalNode.atom_slice(self.pcvInd.atomSlice)
        sub_f.superpose(sub_i, 0, self.pcvInd.align)
        dist_term = md.rmsd(sub_f, sub_i, 0, self.pcvInd.rms)
        match = re.search("stf=.*\n", pars)
        if match is not None:
            self.stf = float(re.split('=', match.group(0).rstrip('\n'))[1])
            if ((self.stf < 1) or (self.stf > (dist_term/self.tolDist/2.5))):
                 print("Straightening factor (stf) is out of range (must be 1 <= stf <= d[0,end]/tolDist )")
                 self.stf = dist_term / self.tolDist / 3
                 print("Setting stf as d[0,end]/tolDist/3: stf=", self.stf)
        else:
            print("Straightening Factor for path reparameterization (stf) not given in \
                    parameter file %s" % (parFile))
            self.stf = dist_term / self.tolDist / 3
            print("Setting stf as d[0,end]/tolDist/3: stf=", self.stf)

        # wall position of PCV-Z for MetaD
        match = re.search("zw=.*\n", pars)
        if match is not None:
            self.zw = float(re.split('=', match.group(0).rstrip('\n'))[1])
            self.zw = (self.zw * self.tolDist) ** 2
        else:
            raise ValueError("Wall position of PCV-Z for MetaDynamics (zw, unit: nm^2) not given in parameter file %s" % (parFile))

        # wall strength of PCV-Z for MetaD
        match = re.search("zwK=.*\n", pars)
        if match is not None:
            self.zwK = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            self.zwK = self.rsTol / (self.tolDist / 20) ** 2
            # raise ValueError("Kappa for wall on PCV-Z is not given for MetaD in parameter file %s" % (parFile))

        # kappa for targeted MD
        match = re.search("kTMD=.*\n", pars)
        if match is not None:
            self.kTMD = int(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            print("Kappa of targeted MD (kTMD) for path reparameterization is not given in \
                parameter file %s" % (parFile))

        # length of targeted MD 
        # default length of targeted MD
        self.lenTMD = 10
        match = re.search("lenTMD=.*\n", pars)
        if match is not None:
            self.lenTMD = float(re.split('=', match.group(0).rstrip('\n'))[1])
        else:
            print("Length of targeted MD (lenTMD) for path reparameterization is not given in \
                parameter file %s" % (parFile))

        # for storing the directories of meta trajecories of last iteration
        self.lastMeta = None
        self.lastSamples = None


    # ==================================================================================================================
    #                                    Prepare directories & files for MetaD
    #                   1. make directories
    #                   2. store node.pdb for sampling under each directory
    #                   3. specify the MetaD length by self.lenSample / path.n_nodes
    # ==================================================================================================================
    def meta_dirs(self, path, dirMeta):
        # input dirMeta is the directory under which, the MetaD sampling and analysis will be performed
        # make sure the path is not empty for MetaD sampling
        if path is None:
            raise ValueError("Path '%s' is empty, can not be sampled" % path.pathName)
        # list to record directories for running
        dirRUNs = []
        for n in range(path.n_nodes):
            dirNode = 'node' + digits(n)
            longDirNode = self.dirRoot + '/' + dirMeta + '/' + dirNode
            if not os.path.exists(longDirNode):
                try:
                    os.makedirs(longDirNode)
                except OSError as error:
                    if error.errno != errno.EEXIST:
                        raise
            nd = path.nodes.slice(n)
            dirRUNs.append(dirNode)
            nodeFile = longDirNode + '/' + self.nodeName
            nd.save(nodeFile)
        # deciding length of each meta using the total amount of sampling
        self.lenMetaD = self.lenSample / path.n_nodes
        return dirRUNs

    # ==================================================================================================================
    #                                   Prepare plumed files for metaD sampling
    #                       1. plumed input file
    #                       2. path pdb file for PCV definition in plumed2 format
    #                       NOTE: engine-specific running files is implemented in prepSampling()
    # ==================================================================================================================
    def meta_setup(self, p_bak, dirMeta, dirRUNs):
        if not os.path.exists(dirMeta):
            os.makedirs(dirMeta)
        for i in range(len(dirRUNs)):
            # prepare MD files
            runDir = self.dirRoot + '/' + dirMeta + '/' + dirRUNs[i]
            #print("###DEBUG### prepSampling")
            self.prepSampling(runDir + '/' + self.nodeName, runDir)
            # prepare plumed path files
            #print("###DEBUG### p.exportPCV")
            #print(dirMeta + '/' + dirRUNs[i])
            p = self.refPath
            p.exportPCV(dirMeta + '/' + dirRUNs[i])
            #print("###EEBUG### p.pcv")
            p.pcv(dirMeta + '/' + dirRUNs[i]) # compute lamda for this path
            # compute self PCV for restraining position on PCV-s
            #   Here we must use the node.pdb generated by meta_dirs(), it is the exact starting conformation
            #   if we use the xtc file, there might be problems
            #print("###DEBUG### md.load")
            node = md.load(runDir + '/' + self.nodeName, top=self.topFile)
            s0,z0 = p.pcv(dirMeta + '/' + dirRUNs[i],node)                        
            # write plumed parameters
            # prepare plumed input file for distance calculation
            pluInput = dirMeta + '/' + dirRUNs[i] + '/' + self.pluName
            # sigma = (p.sw[i+1]-p.sw[i])/2.0 # sigma
            # print('sigma=',sigma)
            f = open(pluInput, 'w+')
            atoms = ''
            for j in range(len(self.pcvInd.atomSlice) - 1):
                atoms = atoms + str(self.pcvInd.atomSlice[j] + 1) + ','
            atoms = atoms + str(self.pcvInd.atomSlice[len(self.pcvInd.atomSlice) - 1] + 1)
            print("WHOLEMOLECULES STRIDE=1 ENTITY0=%s" % atoms, file=f)
            print("p1: PATHMSD REFERENCE=%s LAMBDA=%f NEIGH_STRIDE=4 NEIGH_SIZE=8" \
                  % (p.pathName + '_plu.pdb', p.lamda), file=f)
            print("METAD ARG=p1.sss SIGMA=%f HEIGHT=%f PACE=%d TEMP=%f BIASFACTOR=%d LABEL=metaU" \
                  % (self.sigma, self.gauHeight, self.tauMetaD, self.temp, self.bf), file=f)
            print("UPPER_WALLS ARG=p1.zzz AT=%f KAPPA=%f EXP=2 EPS=1 OFFSET=0 LABEL=zwall" \
                  % (self.zw, self.zwK), file=f)
            print("RESTRAINT ARG=p1.sss KAPPA=%f AT=%f LABEL=res" % (self.kappa, s0), file=f)
            print("PRINT ARG=p1.sss,p1.zzz,metaU.bias,res.bias,zwall.bias STRIDE=" \
                  + str(self.freqTRJ) + " FILE=" + self.cvOut + " FMT=%8.16f", file=f)
            f.close()

    def prepSampling(self, node, dire):
        if self.engine == 'GROMACS':
            gromacs.setup.MD(dire, mdp=self.dirPar + '/' + self.groMDP, mainselection=None, struct=node, \
                             top=self.dirPar + '/' + self.groTOP, deffnm=self.runName, runtime=self.lenMetaD, \
                             dt=self.timeStep, maxwarn=50)
        else:
            raise ValueError("MD engines other than GROMACS are not support yet")

    def prepTMD(self, node, dire):
        if self.engine == 'GROMACS':
            gromacs.setup.MD(dire, mdp=self.dirPar + '/' + self.groMDP, mainselection=None, struct=node, \
                             top=self.dirPar + '/' + self.groTOP, deffnm=self.runName, runtime=self.lenTMD, \
                             dt=self.timeStep, maxwarn=50)
        else:
            raise ValueError("MD engines other than GROMACS are not support yet")


    # ==================================================================================================================
    #                                       perform the actual MetaD sampling
    # ==================================================================================================================
    def meta_sample(self, dirMeta, dirRUNs):
        totTRJ = len(dirRUNs)  # the total number of trajectories to run
        # print("%d trajectories to sample in total." % totTRJ)
        if self.mode == 'serial':
            for i in range(totTRJ):
                runDir = self.dirRoot + '/' + dirMeta + '/' + dirRUNs[i]
                #self.runMeta(runDir)
                runMeta(dire=runDir, engine=self.engine, runName=self.runName, pluName=self.pluName, cvOut=self.cvOut, trjName=self.trjName)
        elif self.mode == 'parallel':
            # record = []
            # lock = multiprocessing.Lock()
            # runDir_list = []
            # cpus = multiprocessing.cpu_count()
            # pool = multiprocessing.Pool(processes=cpus)
            # for i in range(totTRJ):
            #     runDir = self.dirRoot + '/' + dirMeta + '/' + dirRUNs[i]
            #     pool.apply_async(runMeta, args=(runDir, self.engine, self.runName, self.pluName, self.cvOut, self.trjName))
            # pool.close()
            # pool.join()
            N_jobs = totTRJ
            for iter in xrange(0, int(N_jobs / size) + 1):
                tid = iter * size + rank
                if tid < N_jobs:
                    runDir = self.dirRoot + '/' + dirMeta + '/' + dirRUNs[tid]
                    print("+++DEBUG+++ runMeta", tid, "of", N_jobs, ", size", size, ", rank", rank, "iter,", iter)
                    runMeta(dire=runDir, engine=self.engine, runName=self.runName, pluName=self.pluName,
                            cvOut=self.cvOut, trjName=self.trjName)
        elif self.mode == 'qjobs':
            raise ValueError("qjobs version to be implemented")

    # ==================================================================================================================
    #                                    Prepare directories & files for tMD (for reparameterization)
    #                   after inserting nodes, generate an list of nodes to insert
    #                   1. make directories
    #                   2. store node.pdb for sampling under each directory
    #                   3. specify the MetaD length by self.lenSample / path.n_nodes
    # ==================================================================================================================
    def tmd_dirs(self, list_pairs, dirMeta):
        # input dirMeta is the directory under which, the targeted MD sampling will be performed
        if list_pairs is None:
            raise ValueError("list_pairs is empty, no tMD will be performed for path-reparameterization")
        # list to record directories for running
        dirRUNs = []
        for i in range(len(list_pairs)):
            dirPair = 'pair' + digits(i)
            longDirPair = self.dirRoot + '/' + dirMeta + '/tmd4repar/' + dirPair
            if not os.path.exists(longDirPair):
                os.makedirs(longDirPair)
            # store 1st node of the pair as starting conformation for targeted MD
            nd = list_pairs[i].slice(0)
            dirRUNs.append('tmd4repar/' + dirPair)
            nodeFile = longDirPair + '/' + self.nodeName
            nd.save(nodeFile)
            # store 2nd node of the pair as target conformation for targeted MD
            targetFile = longDirPair + '/target.pdb'
            list_pairs[i].slice(1).save_plu2(targetFile,self.pcvInd)
        return dirRUNs

    # ==================================================================================================================
    #                                   Prepare plumed files for targeted MD sampling
    #                       1. plumed input file
    #                       NOTE: engine-specific running files is implemented in prepTMD()
    # ==================================================================================================================
    def tmd_setup(self, dirMeta, dirRUNs):
        if not os.path.exists(dirMeta):
            os.makedirs(dirMeta)
        for i in range(len(dirRUNs)):
            # prepare MD files
            runDir = self.dirRoot + '/' + dirMeta + '/' + dirRUNs[i]
            self.prepTMD(runDir + '/' + self.nodeName, runDir)
            # prepare plumed input file for distance calculation
            pluInput = runDir + '/' + self.pluName
            f = open(pluInput, 'w+')
            atoms = ''
            for j in range(len(self.pcvInd.atomSlice) - 1):
                atoms = atoms + str(self.pcvInd.atomSlice[j] + 1) + ','
            atoms = atoms + str(self.pcvInd.atomSlice[len(self.pcvInd.atomSlice) - 1] + 1)
            print("WHOLEMOLECULES STRIDE=1 ENTITY0=%s" % atoms, file=f)
            print("rmsd: RMSD REFERENCE=target.pdb TYPE=OPTIMAL", file=f)
            print("restraint: ...", file=f)
            print("MOVINGRESTRAINT", file=f)
            print("  ARG=rmsd", file=f)
            print("  AT0=0 STEP0=0 KAPPA0=0", file=f)
            # length of targeted MD
            numSteps = int(self.lenTMD / self.timeStep / 2)
            print("  AT1=0 STEP1=%d KAPPA1=%d" % (numSteps,self.kTMD), file=f)
            print("  AT2=0 STEP2=%d KAPPA2=%d" % (numSteps*2,self.kTMD), file=f)
            print("...", file=f)
            print("PRINT ARG=rmsd STRIDE=" + str(self.freqTRJ) + " FILE=" + self.cvOut + " FMT=%8.16f", file=f)
            f.close()


    # ==================================================================================================================
    #                                       perform the actual MetaD sampling
    # ==================================================================================================================
    def tmd_sample(self, dirMeta, dirTMD):
        totTRJ = len(dirTMD)  # the total number of trajectories to run
        # print("%d trajectories to sample in total." % totTRJ)
        if self.mode == 'serial':
            for i in range(totTRJ):
                runDir = self.dirRoot + '/' + dirMeta + '/' + dirTMD[i]
                runTMD(dire=runDir, engine=self.engine, runName=self.runName, pluName=self.pluName, trjName=self.trjName)
        elif self.mode == 'parallel':
            # record = []
            # lock = multiprocessing.Lock()
            # pool = multiprocessing.Pool(processes=8)
            # runDir_list = []
            # for i in range(totTRJ):
            #     runDir = self.dirRoot + '/' + dirMeta + '/' + dirTMD[i]
            #     pool.apply_async(runTMD, args=(runDir, self.engine, self.runName, self.pluName, self.trjName))
            # pool.close()
            # pool.join()
            N_jobs = totTRJ
            for iter in xrange(0, int(N_jobs / size) + 1):
                tid = iter * size + rank
                if tid < N_jobs:
                    runDir = self.dirRoot + '/' + dirMeta + '/' + dirTMD[tid]
                    runTMD(dire=runDir, engine=self.engine, runName=self.runName, pluName=self.pluName,
                           trjName=self.trjName)
                    print("+++DEBUG+++ runTMD", tid, "of", N_jobs, ", size", size, ", rank", rank, "iter,", iter)
        elif self.mode == 'qjobs':
            raise ValueError("qjobs version to be implemented")
        ## cat all tmd conformations together into one trajectory for reparameterization
        #listTMD = []
        #for i in range(totTRJ):
        #    trjDir = self.dirRoot + '/' + dirMeta + '/' + dirTMD[i]
        #    trjTMD = Confs.traj2conf(md.load(trjDir + '/' + self.trjName , top=self.topFile))
        #    listTMD.append(trjTMD)
        #samplesTMD = Confs.merge(listTMD)
        #trjTMD = self.dirRoot + '/' + dirMeta + '/tmd.xtc'
        #samplesTMD.save(trjTMD)
        ## return the trajectories location of the sampled TMD conformations
        #return trjTMD

    # ==================================================================================================================
    #                                     Analyze metaD trajectories to update path
    #   1. filter metaD trajectories
    #       select only frames whose restraining potential on PCV-s are wall potenial on PCV-z are within tolerance
    #       This helps remove frames with high restraining or wall potential that may be ""unphysical".
    #   2. select median(z) from filtered data, find geometric centroid of the median(z) conformations in each metaD
    #       NOTE: pre-process median(z) conformations for clustering (TODO: to be removed in future versions)
    #   3. reorder median(z) nodes via concorde (Travelling-salesman solver, implemented in C++)
    #   4. path reparameterization (truncate at terminal nodes and insert conformations between distant neighbor nodes)
    # ==================================================================================================================
    def meta_analyze(self, dirMeta, dirSamples):

        # moving to metaD directory
        os.chdir(self.dirRoot + '/' + dirMeta)

        listFilter = []
        list_medz = []
        if rank == 0:
            print("      Filtering metaD trajectories: restraining potential must not exceed ", self.rsTol)
            print("      Finding median(z) conformations")
        print("+++DEBUG+++ meta_analyze", "size", size, ", rank", rank)

        comm.Barrier()
        if rank == 0:
            for i in range(len(dirSamples)):
                trjDir = self.dirRoot + '/' + dirMeta + '/' + dirSamples[i]
                print("        in " + dirSamples[i])
                meta = md.load(trjDir + '/' + self.trjName, top=self.topFile)
                # colvar_filter is generated after sampling
                cvs = np.loadtxt(trjDir + '/colvar_filter', dtype=float)
                z = cvs[:, 2]  # third column is pcv-z
                rsPol = cvs[:, 4]  # fifth column is the value of the restraining potential
                zwPol = cvs[:, 5]  # sixth column is the value of the wall potential on pcv-z
                with np.errstate(invalid='ignore'):
                    cuts = np.where((rsPol >= self.rsTol) | (zwPol >= self.rsTol))[0]
                if len(cuts) > 0:
                    cut = cuts[0]
                    if cut < meta.n_frames:
                        z_filter = z[0:cut]
                    else:
                        z_filter = z[0:meta.n_frames]
                else:
                    cut = meta.n_frames
                    z_filter = z[0:meta.n_frames]
                ranz = np.absolute(np.max(z_filter) - np.min(z_filter)) / 10.0

                # ==========================================================================================================
                #                       find median(z) conformations from filtered trajectories
                #   NOTE: median(z) values given by 'numpy.median(z_filter)' is not an element of the array 'z_filter'
                #   NOTE: This causes troubles when z is unevenly distributed in the array 'z_filter'
                #   NOTE: Therefore, we use an straightforward implementation for median(z) as the following:
                # ==========================================================================================================
                medz = np.sort(z_filter)[len(z_filter) // 2]
                ind_medz = np.where(np.absolute(z_filter - medz) < ranz)[0]

                # ==========================================================================================================
                # extract median(z) confs
                # ==========================================================================================================
                conf_medz = meta.slice(ind_medz)
                list_medz.append(Confs.traj2conf(conf_medz).geoCentroid(self.pcvInd))

                # ==========================================================================================================
                # extract "physical" conformations from  as input for path re-parameterization
                # ==========================================================================================================
                ind_filter = np.array(range(cut))
                filtTRJ = meta.slice(ind_filter)
                # print("      Storing filtered conformations")
                filtTRJ.save(trjDir + '/' + self.trjFilter)

            # ==============================================================================================================
            # store medz centroids per metaD trajectory
            # ==============================================================================================================
            print("      Storing median(z) nodes of this iteration")
            mz = Confs.merge(list_medz)
            pmz = Path('mz_tsp', self.pcvInd, mz)
            pmz.nodes.save(self.dirRoot + '/' + dirMeta + '/mz.xtc')
            pmz.reOrder(truncate=False, dire=self.dirRoot + '/' + dirMeta)
            pmz.exportPath(self.dirRoot + '/' + dirMeta)

            # ==============================================================================================================
            # path re-parameterization, step 1. truncation
            # ==============================================================================================================
            print("      Truncating path: remove segments beyond the two fixed ends")
            p_trunc = pmz.truncate(self.initNode, self.finalNode)
            p_trunc.exportPath(self.dirRoot + '/' + dirMeta)

            # ==============================================================================================================
            # path re-parameterization, step 2. insert conformation between distant neighbor nodes
            # ==============================================================================================================
            # For inserting comformations between distant nodes
            #   filtered samples of both the current and last round should be used as candidates
            #     this is to ensure that there are always sufficient input conformations for path reparameterization
            #       and avoids the path to be broken ( which allows a larger zwall and quicker convergence)
            #         as long as the first rounds gives connected path, this strategy should work fine,
            #           because although this round has broken path, once conformations are inserted from previous round
            #             sampling of next round will definitely include connecting conformations
            # It is also possible to store all sampled data for path-reparameterization, but this is too memory-consuming
            # ==============================================================================================================
            print("      Inserting MetaD conformations into median(z) nodes from the latest two rounds")
            p_in = p_trunc
            for i in range(len(dirSamples)):
                trjDir = self.dirRoot + '/' + dirMeta + '/' + dirSamples[i]
                for data in md.iterload(trjDir + '/' + self.trjFilter, chunk=1000, top=self.topFile):
                    p_in = p_in.insert(data, self.tolDist, self.devMID, self.tolCos)
            if self.lastMeta is not None:
                for i in range(len(self.lastSamples)):
                    trjDir = self.dirRoot + '/' + self.lastMeta + '/' + self.lastSamples[i]
                    for data in md.iterload(trjDir + '/' + self.trjFilter, chunk=1000, top=self.topFile):
                        p_in = p_in.insert(data, self.tolDist, self.devMID, self.tolCos)
            p_in.pathName="mz_tsp_in"
            p_in.exportPath(self.dirRoot + '/' + dirMeta)

            # ==============================================================================================================
            # path re-parameterization, step 3. increase tolDist*=2 shortcut the path
            # ==============================================================================================================
            print("      Straightening path:")
            print("        [a] Short-cutting path by ", self.tolDist, " x ", self.stf)
            p_rc = p_in.rmClose(self.tolDist * self.stf)
            p_rc.pathName = "mz_tsp_in_rc"
            p_rc.exportPath(self.dirRoot + '/' + dirMeta)

            # ==============================================================================================================
            # path re-parameterization, step 4. insert conformation between distant neighbor nodes
            # ==============================================================================================================
            print("        [b] Re-inserting MetaD conformation into straightened path")
            p_in = p_rc
            for i in range(len(dirSamples)):
                trjDir = self.dirRoot + '/' + dirMeta + '/' + dirSamples[i]
                for data in md.iterload(trjDir + '/' + self.trjFilter, chunk=1000, top=self.topFile):
                    p_in = p_in.insert(data, self.tolDist, self.devMID, self.tolCos)
            if self.lastMeta is not None:
                for i in range(len(self.lastSamples)):
                    trjDir = self.dirRoot + '/' + self.lastMeta + '/' + self.lastSamples[i]
                    for data in md.iterload(trjDir + '/' + self.trjFilter, chunk=1000, top=self.topFile):
                        p_in = p_in.insert(data, self.tolDist, self.devMID, self.tolCos)
            p_in.pathName = "mz_tsp_in_st"
            p_in.exportPath(self.dirRoot + '/' + dirMeta)

        # ==============================================================================================================
        # path re-parameterization, step 5. if there are still distant neighbor nodes, extra targeted MD is performed
        # ==============================================================================================================
            p_tmd = p_in
            listFar = p_in.distantNeighbors(self.tolDist)
        else:
            listFar = None
        listFar = comm.bcast(listFar, root=0)
        #comm.Barrier()

        if len(listFar) > 0:
            if rank == 0:
                print("      Distant nodes are still present\n        Perform targeted MD for reparameterization.")
                dirTMDs = self.tmd_dirs(listFar, dirMeta)
                self.tmd_setup(dirMeta,dirTMDs)
            else:
                dirTMDs = None
            dirTMDs = comm.bcast(dirTMDs, root=0)
            comm.Barrier()
            print("###DEBUG### Barrier before tmd_sample")
            dataTMD = self.tmd_sample(dirMeta,dirTMDs)
            comm.Barrier()
            print("###DEBUG### Barrier after tmd_sample")
            # use the tMD samples to insert
            if rank == 0:
                for i in range(len(dirTMDs)):
                    trjDir = self.dirRoot + '/' + dirMeta + '/' + dirTMDs[i]
                    for data in md.iterload(trjDir + '/' + self.trjName, chunk=1000, top=self.topFile):
                        p_tmd = p_tmd.insert(data, self.tolDist, self.devMID, self.tolCos)
        if rank == 0:
            p_tmd.pathName = "mz_tsp_in_st_tmd"
            p_tmd.exportPath(self.dirRoot + '/' + dirMeta)

            # ==============================================================================================================
            # path re-parameterization, step 6. use tolDist to shortcut the path
            # ==============================================================================================================
            print("      Short-cutting path by ", self.tolDist)
            p_meta = p_tmd.rmClose(self.tolDist)
            p_meta.exportPath(self.dirRoot + '/' + dirMeta)

            # ==============================================================================================================
            # store the directory of this round for next iterations
            # ==============================================================================================================
            self.lastMeta = dirMeta
            self.lastSamples = dirSamples

            # return to root directory
            os.chdir(self.dirRoot)
        else:
            p_meta = None
        p_meta = comm.bcast(p_meta, root=0)
        comm.Barrier()
        return p_meta

