from rotamer.topology.residue_sidechains import amino_acids
from rotamer.dihedral.find_dihedrals import Dihedral, map_dihedrals
from rotamer.io.gmin import read_lowest
import rotamer.topology.read_amber_prmtop as ra
import networkx as nx

class BaseRotamerState(object):
    """
    Base class for rotamer states for the GMIN rotamer library.
    """
    pass

class FullRotamerState(BaseRotamerState):
    """
    Defines a rotamer state for use with the GMIN rotamer library.

    A rotamer state corresponds to a single minimum from the rotamer library simulations. Each state has the values of
    all the dihedrals in the corresponding tripeptide and the energy (relative to the global minimum for that
    tripeptide).
    """
    pass


class RotamerStateFactory(object):
    """
    Creates new rotamer states, either from raw GMIN data or from a HDF5 database.
    """
    def from_lowest(self, prmtop_filename, lowest_filename, state_type):
        """
        Creates a new RotamerState of appropriate type from a GMIN lowest file.

        :param prmtop_filename: File name of the AMBER topology file
        :param lowest_filename: File name of the GMIN lowest file
        :param state_type: Type of rotamer state (i.e. what sort of dependence on neighbours, backbone dihedrals etc.)
        :return:
        """
        molecule = ra.parse_topology_file(prmtop_filename)
        print molecule.identify_residues()[1]
        # print map_dihedrals(*(molecule.identify_residues()))
        configs = read_lowest(lowest_filename)


if __name__ == "__main__":
    test = RotamerStateFactory()
    test.from_lowest("/home/khs26/coords.prmtop", "/home/khs26/lowest", 1)