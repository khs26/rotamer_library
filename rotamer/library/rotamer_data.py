from rotamer.topology.residue_sidechains import amino_acids
import rotamer.dihedral.find_dihedrals as find_dihedrals
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
        # molecule.identify_residues()[1]
        find_dihedrals.map_dihedrals(*(molecule.identify_residues()))
        configs = read_lowest(lowest_filename)
        find_dihedrals.phi_psi_dihedrals(molecule)
        find_dihedrals.sidechain_dihedrals(molecule)
        for res in sorted(molecule.residues.nodes()):
            print res.name, res.dihedrals


if __name__ == "__main__":
    test = RotamerStateFactory()
    test.from_lowest("./coords.prmtop", "./lowest", 1)