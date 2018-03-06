from __future__ import print_function, division
import numpy as np
from mdtraj.utils import ilen
import mdtraj.formats.pdb.pdbfile as pdb

def _format_83(f):
    """Format a single float into a string of width 8, with ideally 3 decimal
    places of precision. If the number is a little too large, we can
    gracefully degrade the precision by lopping off some of the decimal
    places. If it's much too large, we throw a ValueError"""
    if -999.999 < f < 9999.999:
        return '%8.3f' % f
    if -9999999 < f < 99999999:
        return ('%8.3f' % f)[:8]
    raise ValueError('coordinate "%s" could not be represnted '
                     'in a width-8 field' % f)

class PluPDBfile(pdb.PDBTrajectoryFile):
    def __init__(self, filename, mode='r', force_overwrite=True, standard_names=True):
        pdb.PDBTrajectoryFile.__init__(self, filename, mode, force_overwrite)  #, standard_names)
        # pdb.PDBTrajectoryFile.__init__(self, filename, mode, force_overwrite, standard_names)
        # do not write footer; footer will add an extra END at the last line and causes crash of plumed
        self._footer_written = True

    def write(self, positions, topology, frame_ind=None, pcv_ind=None, unitcell_lengths=None,
              unitcell_angles=None):
        """Write a PDB file to disk using plumed2 PCV format

        Parameters
        ----------
        positions : array_like
            The list of atomic positions to write.
        topology : mdtraj.Topology
            The Topology defining the model to write.
        frame_ind : {int, None}
            If not None, the index of frames will be surrounded by REMARK X=? and END
        unitcell_lengths : {tuple, None}
            Lengths of the three unit cell vectors, or None for a non-periodic system
        unitcell_angles : {tuple, None}
            Angles between the three unit cell vectors, or None for a non-periodic system
        """
        if not self._mode == 'w':
            raise ValueError('file not opened for writing')
        if not self._header_written:
            self._write_header(unitcell_lengths, unitcell_angles)
            self._header_written = True

        if ilen(topology.atoms) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if np.any(np.isnan(positions)):
            raise ValueError('Particle position is NaN')
        if np.any(np.isinf(positions)):
            raise ValueError('Particle position is infinite')

        self._last_topology = topology  # Hack to save the topology of the last frame written, allows us to output CONECT entries in write_footer()

        posIndex = 0
        if frame_ind is not None:
            print("REMARK X=%d" % frame_ind, file=self._file)
            for (chainIndex, chain) in enumerate(topology.chains):
                chainName = self._chain_names[chainIndex % len(self._chain_names)]
                residues = list(chain.residues)
                for (resIndex, res) in enumerate(residues):
                    if len(res.name) > 3:
                        resName = res.name[:3]
                    else:
                        resName = res.name
                    for atom in res.atoms:
                        if len(atom.name) < 4 and atom.name[:1].isalpha() and (
                                        atom.element is None or len(atom.element.symbol) < 2):
                            atomName = ' ' + atom.name
                        elif len(atom.name) > 4:
                            atomName = atom.name[:4]
                        else:
                            atomName = atom.name
                        coords = positions[posIndex]
                        if atom.element is not None:
                            symbol = atom.element.symbol
                        else:
                            symbol = ' '
                        line = "ATOM  %5d %-4s %3s %1s%4d    %s%s%s %5.2f %5.2f      %-4s%-2s  " % (
                            pcv_ind.atomInd[posIndex] % 100000, atomName, resName, chainName,
                            (res.resSeq) % 10000, _format_83(coords[0]),
                            _format_83(coords[1]), _format_83(coords[2]),
                            pcv_ind.alignPLU[posIndex], pcv_ind.rmsPLU[posIndex], atom.segment_id[:4], symbol[-2:])
                        assert len(line) == 80, 'Fixed width overflow detected'
                        print(line, file=self._file)
                        posIndex += 1
                    if resIndex == len(residues) - 1:
                        print("END", file=self._file)


