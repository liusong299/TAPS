from __future__ import print_function, division
import mdtraj as md
import numpy as np
import PluPDB as plu
from mdtraj.utils import in_units_of
from copy import deepcopy

class Confs(md.Trajectory):  # inherit class Trajectory in mdtraj
    """Container object for a molecular dynamics conformations, built upon Class Trajectory from MDTraj

    # ==============
    Child Attributes
    # ==============
    gc : int [the index for the geometrically central confomration within this trajectory]

    # ==============
    Child Method
    # ==============
    atom_slice() overrides parent class mdtraj.Trajectory
    save_plu2()
    geoCentroid()
    """
    @property
    def gc(self):
        """
        The index of geometric central conformation in the this Trajectory
        :return:
        gc : int
            The index of geometric central conformation in the this Trajectory
        """
        return self._gc

    @gc.setter
    def gc(self, value):
        """
        The index of geometric central conformation in the this Trajectory
        :return:
        gc : int
            The index of geometric central conformation in the this Trajectory
        """
        self._gc = value

    def __init__(self, xyz, topology, gc=None, time=None, unitcell_lengths=None, unitcell_angles=None):

        # call constructor of parent class md.Trajectory
        """
        :rtype: Confs
        """
        super(Confs,self).__init__(xyz, topology, time, unitcell_lengths, unitcell_angles)
        #md.Trajectory.__init__(self, xyz, topology, time, unitcell_lengths, unitcell_angles) # equal to last line

        # geometric center has no default
        self._gc = gc

    # ==================================================================================================================
    # convert a mdtraj.Trajectory object into Confs to enable TAPS-relevant functions
    # ==================================================================================================================
    @classmethod
    def traj2conf(self,trj):
        tmp=Confs(trj.xyz,trj.topology,None)
        tmp.__dict__=trj.__dict__
        tmp.gc=None
        return tmp

    # ==================================================================================================================
    # merges a few confs/mdtraj.trajectories into one, regardless of how many frames are in each Conf
    # ==================================================================================================================
    @classmethod
    def merge(self, listConfs):
        if len(listConfs)>0:
            tot = 0
            merged = deepcopy(listConfs[0])
            for i in range(len(listConfs)):
                n = listConfs[i].n_frames
                listConfs[i].time = range(tot, tot + n)
                if i > 0:
                    merged=merged.join(listConfs[i])
                tot += n
            return merged
        else:
            return None

    # ==================================================================================================================
    # find the geometric centroid of a Confs object containing multiple conformations, return the gc as one-frame Confs
    # ==================================================================================================================
    def geoCentroid(self, pcvInd=None, doPBC=False):
        """
        find the geometric central conformation and extract it into a new Confs instance
        and update the index into

        pcvInd

        :return:
        gcc : Confs
        """
        if pcvInd is None:
            raise ValueError("Atoms for defining PathCV not given")

        else:
            if self.n_frames > 1:
                if doPBC:
                    self.image_molecules(True)
                pluAtoms = self.atom_slice(pcvInd.atomSlice)
                # numpy array to store rowsum
                rowSum = np.zeros((self.n_frames, 1))
                for frame in range(0, self.n_frames):
                    # align frames using alignIndex
                    pluAtoms.superpose(pluAtoms,frame,pcvInd.align)  # here align and rmsd are both for sliced trajectory
                    # compute RMSD
                    val = md.rmsd(pluAtoms, pluAtoms, frame, pcvInd.rms)
                    rowSum[frame] = np.sum(val)
                # store index of frame in geometric Center in self._gc
                self._gc = np.argmin(rowSum)
                # extract the central frame
                tmp = self.slice(self._gc)
                # return a Confs object containing the central frame
                gcc = deepcopy(tmp)
            else:
                gcc = deepcopy(self)
            return gcc

    # ==================================================================================================================
    # Export the atoms defining PCV as a PDB file in plumed format
    # ==================================================================================================================
    def save_plu2(self, filename, pcvInd=None, force_overwrite=True):
        """Save trajectory to plumed PDB format
        Parameters
        ----------
        filename : str
            filesystem path in which to save the trajectory
        force_overwrite : bool, default=True
            Overwrite anything that exists at filename, if its already there
        bfactors : array_like, default=None, shape=(n_frames, n_atoms) or (n_atoms,)
            Save bfactors with pdb file. If the array is two dimensional it should
            contain a bfactor for each atom in each frame of the trajectory.
            Otherwise, the same bfactor will be saved in each frame.
        """
        self._check_valid_unitcell()

        if pcvInd is None:
            raise ValueError("Atoms for defining PCV not given")

        #substract the plumed atoms from original trajectory
        pluAtoms=self.atom_slice(pcvInd.atomSlice)

        if len(pcvInd.atomInd) != pluAtoms.n_atoms:
            raise ValueError("number of atom index %s should equal n_atoms %s" % str(len(pcvInd.atomInd)), str(pluAtoms.n_atoms))
        if len(pcvInd.alignPLU) != pluAtoms.n_atoms:
            raise ValueError("number of atoms to align %s should equal n_atoms %s" % str(len(pcvInd.alignPLU)), str(pluAtoms.n_atoms))
        if len(pcvInd.rmsPLU) != pluAtoms.n_atoms:
            raise ValueError("number of atoms for rmsd %s should equal n_atoms %s" % str(len(pcvInd.rmsPLU)), str(pluAtoms.n_atoms))

        with plu.PluPDBfile(filename, 'w', force_overwrite=force_overwrite) as f:
            for i in xrange(pluAtoms.n_frames):
                f.write(in_units_of(pluAtoms._xyz[i], Confs._distance_unit, f.distance_unit),
                            pluAtoms.topology,
                            frame_ind=(i+1),
                            pcv_ind=pcvInd)
