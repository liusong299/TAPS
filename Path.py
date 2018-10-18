#=======================================================================================================================
#                                                   global imports
#=======================================================================================================================
from __future__ import print_function, division
import numpy as np
import mdtraj as md
import os
import re
from Confs import Confs
from copy import deepcopy
import glob
import subprocess

#=======================================================================================================================
#                               global names for plumed command: 'plumed' and 'driver'
#=======================================================================================================================
#plu="srun --ntasks=1 --hint=nomultithread --ntasks-per-node=1 --ntasks-per-socket=1 --ntasks-per-core=1 --mem_bind=v,local plumed"
#  TODO: put binary name of 'plumed' as input
plu="plumed"
dri="driver" # TODO: put binary name of 'driver' as input

class Path(object):
    # ==================================================================================================================
    # constructor: one input at least - the name of the path
    # ==================================================================================================================
    def __init__(self, pName="path_node", pInd = None, nodes = None):
        self.pathName = pName
        self.pcvInd = deepcopy(pInd)
        self.nodes = deepcopy(nodes)
        if (pInd is not None) and (nodes is not None):
            self.n_nodes=self.nodes.n_frames
            # self.pcv()
        else:
            self.n_nodes = None
            self.lamda = None
            self.pcvs = None
            # self.sw = None

    # # ==================================================================================================================
    # # compute the location of wall potential on pcv-s for each node
    # # ==================================================================================================================
    # def sWall(self):
    #     self.sw = (self.pcvs[0:(self.n_nodes-1)]+self.pcvs[1:(self.n_nodes)])/2.0
    #     self.sw = np.append(np.insert(self.sw, 0, 1.0), self.n_nodes)
    # ==================================================================================================================
    # compute lamda, PCV using plumed (METAPROGRAMING), assuming plumed is correctly installed
    # ==================================================================================================================
    def pcv(self, dire = None, trj = None):

        # temporary file names
        
        tmpCVfile = dire + '/COLVAR'
        # nodeFile = "node.pdb"
        pluPdbFile = dire + "/tmpFrames.pdb"
        self.nodes.save_plu2(pluPdbFile,self.pcvInd)

        # if not trajectory is given, compute lamda and sw, etc for self
        if trj is None:
            # compute lamda
            pluAtoms = self.nodes.atom_slice(self.pcvInd.atomSlice)
            # temporary array to store all neighbor RMSD
            nb = np.zeros(((self.n_nodes-1), 1))
            for frame in range(0, (self.n_nodes-1)):
                # align frames using alignIndex
                pluAtoms.superpose(pluAtoms,frame, self.pcvInd.align)  # here align and rmsd are both for sliced trajectory
                # compute RMSD
                rms = md.rmsd(pluAtoms, pluAtoms, frame, self.pcvInd.rms)
                nb[frame]=rms[(frame+1)]
            avNB2=np.divide(np.sum(np.square(nb)),frame)
            self.lamda=2.3/avNB2

            # prepare plumed input file for distance calculation
            pluInput =  dire + "/" + self.pathName + ".plu"
            f = open(pluInput, 'w+')
            atoms = ''
            for i in range(len(self.pcvInd.atomSlice) - 1):
                atoms = atoms + str(self.pcvInd.atomSlice[i] + 1) + ','
            atoms = atoms + str(self.pcvInd.atomSlice[len(self.pcvInd.atomSlice) - 1] + 1)
            print("WHOLEMOLECULES STRIDE=1 ENTITY0=%s" % atoms, file=f)
            print('p1: PATHMSD REFERENCE=%s LAMBDA=%f NEIGH_STRIDE=4 NEIGH_SIZE=8' % \
                  (pluPdbFile, self.lamda), file=f)
            print('PRINT ARG=p1.sss,p1.zzz STRIDE=1 FILE=%s FMT=%s' % (tmpCVfile, '%8.8f'), file=f)
            f.close()

            # write nodes to disk
            trjFile = dire + "/tmp.xtc"
            self.nodes.save(trjFile)

            # launch plumed for pcv calculation
            tmpOutFile = dire + "/pluOut.tmp"
            cmd = plu + " " + dri + " --mf_xtc " + trjFile + " --plumed " + pluInput \
                  + " 1>" + tmpOutFile + " 2>" + tmpOutFile
            os.system(cmd)

            # read plumed output
            fr = open(tmpCVfile, 'r+')
            lines = fr.read()
            fr.close()
            fw = open(tmpCVfile, 'w')
            # remove any line that includes '#'
            pdata = re.sub("#.*\n", "", lines)
            print(pdata, file=fw)
            fw.close()
            # store in self.pcvs
            cvs = np.loadtxt(tmpCVfile, dtype=float)
            self.pcvs = cvs[:, 1]
            self.pcvz = cvs[:, 2]
            os.remove(tmpCVfile)
            os.remove(tmpOutFile)
            os.remove(pluInput)
            os.remove(trjFile)
            os.remove(pluPdbFile)

        # if a trajectory is given, compute the PCV-s and PCV-z of this trajectory
        else:
            if isinstance(trj,md.Trajectory):
                # write md.trajectory to disk
                trjFile = "tmp.xtc"
                trj.save(trjFile)

                # prepare plumed input file for distance calculation
                pluInput = self.pathName + ".plu"
                f = open(pluInput, 'w+')
                atoms = ''
                for i in range(len(self.pcvInd.atomSlice) - 1):
                    atoms = atoms + str(self.pcvInd.atomSlice[i] + 1) + ','
                atoms = atoms + str(self.pcvInd.atomSlice[len(self.pcvInd.atomSlice) - 1] + 1)
                print("WHOLEMOLECULES STRIDE=1 ENTITY0=%s" % atoms, file=f)
                print('p1: PATHMSD REFERENCE=%s LAMBDA=%f NEIGH_STRIDE=4 NEIGH_SIZE=8' % \
                      (pluPdbFile, self.lamda), file=f)
                print('PRINT ARG=p1.sss,p1.zzz STRIDE=1 FILE=%s FMT=%s' % (tmpCVfile, '%8.8f'), file=f)
                f.close()

                # launch plumed for pcv calculation
                tmpOutFile = 'pluOut.tmp'
                cmd = plu + " " + dri + " --mf_xtc " + trjFile + " --plumed " + pluInput \
                      + " 1>" + tmpOutFile + " 2>" + tmpOutFile
                os.system(cmd);

                # read plumed output
                fr = open(tmpCVfile, 'r+')
                lines = fr.read()
                fr.close()
                fw = open(tmpCVfile, 'w')
                # remove any line that includes '#'
                pdata = re.sub("#.*\n", "", lines)
                print(pdata, file=fw)
                fw.close()
                # store in self.pcvs
                cvs = np.loadtxt(tmpCVfile, dtype=float)
                if trj.n_frames == 1:
                    s = cvs[1]
                    z = cvs[2]
                else:
                    s = cvs[:, 1]
                    z = cvs[:, 2]
                os.remove(tmpCVfile)
                os.remove(tmpOutFile)
                os.remove(pluInput)
                os.remove(trjFile)
                os.remove(pluPdbFile)
                return s, z
            else:
                raise ValueError("nodes must be an instance of mdtraj.Trajectory or TAPS.Confs")

    # ==================================================================================================================
    # load path from a trajectory, computes PCV of this path automatically after loading
    # ==================================================================================================================
    def loadFromTRJ(self, trajName, topName):
        self.nodes = Confs.traj2conf(md.load(trajName, top=topName))
        self.n_nodes = self.nodes.n_frames
        # self.pcv()

    # ==================================================================================================================
    # export each node as pdb for MD/MetaD, must include all atoms of the system
    # ==================================================================================================================
    def exportFrames(self):
        for i in range(self.n_nodes):
            nd=self.nodes.slice(i)
            nd.save(self.pathName+'_node'+str(i)+'.pdb')

    # ==================================================================================================================
    # export all nodes into a trajectory file under a specified directory
    # ==================================================================================================================
    def exportPath(self, dire):
        self.nodes.save(dire + '/' + self.pathName+'.xtc')

    # ==================================================================================================================
    # export all nodes as a Plumed pdb file for PCV computation
    # ==================================================================================================================
    def exportPCV(self, dire):
        pluPdbFile = self.pathName + "_plu.pdb"
        if type(self.nodes) is Confs:
            self.nodes.save_plu2(dire+'/'+pluPdbFile,self.pcvInd)
        else:
            raise ValueError("nodes must be an instance of the class TAPS/Confs")

    # ==================================================================================================================
    # use concorde to generate new order of the nodes
    # ==================================================================================================================
    def reOrder(self, truncate=False, source=None, target=None, doPBC=False, dire=None):
        # generate RMSD matrix for concorde
        dist = np.zeros((self.n_nodes+1, self.n_nodes+1))
        if doPBC:
            self.nodes.image_molecules()
        pluAtoms = self.nodes.atom_slice(self.pcvInd.atomSlice)
        pluAtoms.superpose(pluAtoms, 0, self.pcvInd.align)
        for i in range(self.n_nodes):
            dist[i][0:self.n_nodes] = md.rmsd(pluAtoms, pluAtoms, i, self.pcvInd.rms)
        np.savetxt('rmsd.conc', dist*1000, fmt='%d')
        fc = open('head.conc', 'w')
        print("NAME: RMSD",file=fc)
        print("TYPE: TSP", file=fc)
        print("DIMENSION: %d" % (self.n_nodes+1), file=fc)
        print("EDGE_WEIGHT_TYPE: EXPLICIT", file=fc)
        print("EDGE_WEIGHT_FORMAT: FULL_MATRIX", file=fc)
        print("EDGE_WEIGHT_SECTION:", file=fc)
        fc.close()

	# ==============================================================================================================
	# combine two files as input for concorde
	# ==============================================================================================================
	rmsdINT = self.pathName + '.conc'
        files = ['head.conc', 'rmsd.conc']
        with open(rmsdINT, 'w') as combine:
            for file_ in files:
                for line in open(file_, 'r'):
                    combine.write(line)
        cmd = "concorde %s" % rmsdINT
        os.system(cmd) # run concorde

        result = self.pathName + '.sol'
        fr = open(result,'r')
        lines = fr.readlines()
        fr.close()
        lines.pop(0) # remove the first line
        tmp=''.join(lines).replace('\n', '').split(' ')
        tmp.remove('')
        nds=np.array(map(int,tmp),dtype=int)
        brk=np.argmax(nds) # maxvalue is the virtual node
        order = np.append(nds[(brk + 1):len(nds)], nds[0:brk])
        if truncate:
            if source is None:
                raise ValueError("To truncate path, source node must be provided")
            if target is None:
                raise ValueError("To truncate path, target node must be provided")
            si = np.where(order == source)[0]
            ti = np.where(order == target)[0]
            if si > ti:
                order = order[ti:(si+1)]
            else:
                order = order[si:(ti+1)]
        newList = []
        for i in range(len(order)):
            newList.append(self.nodes.slice(order[i]))
        self.nodes = Confs.merge(newList)
        self.n_nodes = len(self.nodes)
        self.pcv(dire=dire)
        # clean up concorde files
        for f in glob.glob("*.sav"):
            os.remove(f)
        for f in glob.glob("*.pul"):
            os.remove(f)
        for f in glob.glob("*.mas"):
            os.remove(f)
        for f in glob.glob("*.conc"):
            os.remove(f)
        for f in glob.glob("*.sol"):
            os.remove(f)
        return order

    # ==================================================================================================================
    # Truncate the nodes of path beyond the two fixed terminals
    # ==================================================================================================================
    def truncate(self, iNode, fNode):
        if (iNode is None) or (fNode is None):
            raise ValueError("the two terminal nodes must be provided for truncation")
        else:
            initNode = deepcopy(iNode)
            finNode = deepcopy(fNode)
            sub_nodes=self.nodes.atom_slice(self.pcvInd.atomSlice)
            sub_i = initNode.atom_slice(self.pcvInd.atomSlice)
            sub_f = finNode.atom_slice(self.pcvInd.atomSlice)
            sub_nodes.superpose(sub_i,0,self.pcvInd.align)
            r1 = md.rmsd(sub_nodes,sub_i,0,self.pcvInd.rms)
            si = np.argmin(r1)
            sub_nodes.superpose(sub_f,0, self.pcvInd.align)
            r2 = md.rmsd(sub_nodes,sub_f,0,self.pcvInd.rms)
            ti = np.argmin(r2)
            if si > ti:
                # add self.initNode and self.finNode at terminals
                listNodes = []
                listNodes.append(finNode)
                listNodes.append(self.nodes.slice(np.arange(ti, si + 1)))
                listNodes.append(initNode)
                newNodes = Confs.merge(listNodes)
            else:
                # add self.initNode and self.finNode at terminals
                listNodes = []
                listNodes.append(initNode)
                listNodes.append(self.nodes.slice(np.arange(si, ti + 1)))
                listNodes.append(finNode)
                newNodes = Confs.merge(listNodes)
            newPath = Path(self.pathName+'_tr',self.pcvInd,newNodes)
        return newPath

    # ==================================================================================================================
    # Shorten the path by skiping nodes i-j if d[i-1,j+1] is already shorter than tolerable maximum distance
    #       finding short-cuts (if any) avoid unnecessary curvature and loop;
    #       This will make the path as straight as the sampled data allows
    # ==================================================================================================================
    def rmClose(self,tolDist,doPBC=False):
        # prepare for rmsd computation
        pathConfs = self.nodes
        if doPBC:
            pathConfs.image_molecules()
        sub_confs = pathConfs.atom_slice(self.pcvInd.atomSlice)
        sub_confs.superpose(sub_confs, 0, self.pcvInd.align)
        i = 0  # starting from the first node
        while (i < pathConfs.n_frames):
            distances = md.rmsd(sub_confs, sub_confs, i, self.pcvInd.rms)
            # make initial self-distance ultra-large
            distances[i]=1000000000
            shortcuts = np.where(distances <= tolDist)[0]
            # print('number of nodes in path: ', pathConfs.n_frames)
            # print('short-cuts = ',shortcuts)
            if (len(shortcuts) <= 0):  # no short-cuts found, move to the next node
                i += 1
            else:
                scID = np.max(shortcuts)
                if (scID>i):
                    # print('node', i, ' short-cut to node', scID)
                    # scID is not yet the last frame, cut the frames inbetween out
                    if (scID < (pathConfs.n_frames - 1)):
                        cuts = np.append(np.array(range(i + 1)), np.array(range(scID, pathConfs.n_frames)))
                        # print('cuts=',cuts)
                        sub_confs = sub_confs.slice(cuts)
                        pathConfs = pathConfs.slice(cuts)
                        i += 1
                    # already the last frame in path, jump out of loop
                    else:
                        cuts = np.append(np.array(range(i + 1)), np.array(range(scID, pathConfs.n_frames)))
                        sub_confs = sub_confs.slice(cuts)
                        pathConfs = pathConfs.slice(cuts)
                        # print('cuts=', cuts)
                        break
                else:
                    i+=1 # closest node appears before the current node; move to the next node
        newPath = Path(self.pathName + '_rc', self.pcvInd, pathConfs)
        return newPath

    # ==================================================================================================================
    # Insert conformations between neighbor nodes i,i+1, if d[i,i+1] > tolerable maximum distance
    #   1. ensure conformation f is closer to node i and i+1 than other nodes on path
    #   2. ensure conformation f is on the vertical lines between node i,i+1 |d[f,i]-d[f,i+1]| < (d[i,i+1]*devMID)
    #   3. ensure  cosTheta = (d[i, f]*2 + d[i+1,f]*2 - d[i,i+1]) < tolCos
    #   4. finding short-cuts (if any) avoid unnecessary curvature and loop, using function rmClose
    # ==================================================================================================================
    def insert(self, traj, tolDist=0.015, devMID=0.1, tolCos=0, doPBC=False):
        # traj is a mdTraj.Trajectory object or the Confs object

        # copy current nodes into a path, prepare for rmsd computation
        pathConfs = deepcopy(self.nodes)
        if doPBC:
            pathConfs.image_molecules()
        sub_opath = pathConfs.atom_slice(self.pcvInd.atomSlice)
        sub_opath.superpose(sub_opath, 0, self.pcvInd.align)

        # prepare data for rmsd computation
        # align sampled data to the first node of path
        trj = traj
        if doPBC:
            trj.image_molecules()
        sub_data = trj.atom_slice(self.pcvInd.atomSlice)
        sub_data.superpose(sub_opath, 0, self.pcvInd.align)

        # scan the initial path, label the neighbors whose distance is larger than tolDist
        tmp = []
        for i in range((sub_opath.n_frames - 1)):
            dist = md.rmsd(sub_opath.slice(i), sub_opath.slice(i + 1), 0, self.pcvInd.rms)
            if dist > tolDist:
                tmp.append(i)
        toInsert = np.array(tmp)
        # print('toInsert=',toInsert)

        # remove the path node from dataset
        for i in range(sub_opath.n_frames):
            dist = md.rmsd(sub_data, sub_opath.slice(i), 0, self.pcvInd.rms)
            ind_zero = np.where(dist > 0)[0]
            sub_data = sub_data.slice(ind_zero)
            trj = trj.slice(ind_zero)

        while len(toInsert) > 0:
            i1 = toInsert[0]
            i2 = i1 + 1
            # print('Searching for conformations to insert between node %d & %d' % (i1,i2))
            # Criterion 1: conf to insert must be closer to the i1 or i2 than any other nodes
            # loop over all path nodes
            dists = np.empty([sub_opath.n_frames, sub_data.n_frames])
            for i in range(sub_opath.n_frames):
                dists[i] = md.rmsd(sub_data, sub_opath, i, self.pcvInd.rms)
            # logic_close = np.logical_or((np.argmin(dists, axis=0) == i1), (np.argmin(dists, axis=0) == i2))
            # print(' logic_close=: ', np.where(logic_close)[0].size)
            r1 = dists[i1]
            r2 = dists[i2]
            ind_data = np.arange(sub_data.n_frames)
            r12 = md.rmsd(sub_opath.slice(i1), sub_opath.slice(i2), 0, self.pcvInd.rms)
            # Criterion 2: |r1-r2|< r12/scale_r12, i.e. conf to insert must be in the middle of i1,i2
            scale_r12 = 0.01
            while (scale_r12 < devMID):
                logic_mid = (np.absolute(r1 - r2) < (r12 * scale_r12))  # r1-r2->0
                # print('  # logic_mid at scale_r12=', scale_r12, ': ',np.where(logic_mid)[0].size)
                cosTheta = (r1 ** 2 + r2 ** 2 - r12 ** 2) / (2 * r1 * r2)
                # print('cosTheta=',cosTheta)
                # Criterion 3: cos(Theta) must be smaller than preset-value
                logic_cos = (cosTheta < tolCos)  # cosTheta < 0, theta > 90 degrees
                # logic_candidate = ((logic_close & logic_mid) & logic_cos)
                logic_candidate = (logic_mid & logic_cos)
                # print('  # logic_both at scale_r12=', scale_r12, ': ',np.where(logic_candidate)[0].size)
                n_candidate = np.where(logic_candidate)[0].size
                if n_candidate > 0:
                    # if jump out of loop and proceed
                    break
                else:
                    scale_r12 *= 2
            # print('number of candidates: %d' % n_candidate)
            # if there is data meeting both criteria, select the min(cos(Theta))
            if n_candidate > 0:
                ind_candidate = ind_data[logic_candidate]
                sel = ind_candidate[np.argmin(cosTheta[logic_candidate])]  # find the candidate index in original data
                # print('Found conformations between %d-%d, in data[%d]' % (i1, i2, sel))
                tmp = []
                tmp.append(Confs.traj2conf(pathConfs.slice(range(i2))))
                tmp.append(Confs.traj2conf(trj.slice(sel)))
                tmp.append(Confs.traj2conf(pathConfs.slice(range(i2, len(pathConfs)))))
                pathConfs = Confs.merge(tmp)

                # PathConfs.insert(i2, trj.slice(sel))
                tmp = []
                tmp.append(Confs.traj2conf(sub_opath.slice(range(i2))))
                tmp.append(Confs.traj2conf(sub_data.slice(sel)))
                tmp.append(Confs.traj2conf(sub_opath.slice(range(i2, len(sub_opath)))))
                sub_opath = Confs.merge(tmp)
                # print('inserted pathConfs=', pathConfs)
                # remove the conformation from dataset
                ind_data = np.delete(ind_data, sel)
                # print('Removing selected conformation from dataset, %d conformations left' % len(ind_data))
                trj = trj.slice(ind_data)
                sub_data = sub_data.slice(ind_data)
                # update the list toInsert
                # if r1,r2 > tolDist, extra insertion is still needed in both pairs
                if ((r1[sel] > tolDist) and (r2[sel] > tolDist)):
                    # print(r1[sel],'>',tolDist,'\t',r2[sel],'>',tolDist)
                    toInsert[1:len(toInsert)] = toInsert[1:len(toInsert)] + 1
                    toInsert = np.insert(toInsert, 1, i2)
                # if r1,r2 <= tolDist, extra insertion is not necessary
                elif ((r1[sel] <= tolDist) and (r2[sel] <= tolDist)):
                    # print(r1[sel], '<=', tolDist, '\t', r2[sel], '<=', tolDist)
                    toInsert = toInsert + 1
                    toInsert = np.delete(toInsert, 0) 

                # if r1<=tolDist, r2>tolDist, only insert between newly inserted conf and i2
                elif ((r1[sel] <= tolDist) and (r2[sel] > tolDist)):
                    # print(r1[sel], '<=', tolDist, '\t', r2[sel], '>', tolDist)
                    toInsert = toInsert + 1
                # if r1>tolDist, r2<=tolDist, only insert between i1 and newly inserted
                else:
                    # print(r1[sel], '>', tolDist, '\t', r2[sel], '<', tolDist)
                    toInsert[1:len(toInsert)] = toInsert[1:len(toInsert)] + 1
                    # print(toInsert)
            else:
                # print('No conformation found between %d-%d' % (i1,i2))
                toInsert = np.delete(toInsert, 0)
                # print('after trial toInsert=',toInsert)
                # print('')

        return Path(self.pathName + '_in', self.pcvInd, pathConfs)

    # ==================================================================================================================
    # Find neighbours that are distant (rmsd > cutoff)
    #   return a list of mdtraj.trajectory object, each containing only the two distant neighbor nodes
    # =================================================================================================================
    def distantNeighbors(self, tolDist=0.015, doPBC=False):
        listDistant = []
        # copy path nodes
        nodes = deepcopy(self.nodes)
        n = nodes.n_frames
        # compute rmsd between neighbor nodes
        if doPBC:
            nodes.image_molecules()
        sub_nodes = nodes.atom_slice(self.pcvInd.atomSlice)
        sub_nodes.superpose(sub_nodes, 0, self.pcvInd.align)
        distances = np.zeros((n, n), dtype=np.float)
        for i in range(n-1):
            dist = md.rmsd(sub_nodes.slice(i), sub_nodes.slice(i+1), 0, self.pcvInd.rms)
            if dist > tolDist:
                listDistant.append(nodes.slice([i,i+1]))
        return listDistant
