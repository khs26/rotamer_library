import networkx as nx

import rotamer.topology.read_amber_prmtop as amber
import rotamer.topology.chirality as chirality
import transforms
from measure_dihedral import dihedral_with_symmetry, symmetric_atoms
from rotamer.io.amber import read_amber_restart


class Dihedral(object):
    def __init__(self, dihedral_atoms):
        self.atoms = dihedral_atoms
        self.residue = dihedral_atoms[0].residue
        self.atom_graph = dihedral_atoms[0].molecule.atoms
        # self.atom_map = {v: k for k, v in self.residue.atom_map.items()}
        self.atom_map = self.residue.atom_map
        self.get_moving_atoms()

    def measure_dihedral(self, coords):
        return dihedral_with_symmetry(coords, self)

    def set_dihedral(self, coords, value):
        """ Returns coordinates which set the dihedral to a particular angle, value.

        :param coords: Coordinates of the system
        :param value: Target angle for the dihedral (-pi < x < pi)
        :return new_coords: New coordinates (with chosen dihedral angle)
        """
        current_angle = self.measure_dihedral(coords)
        rotation_angle = value - current_angle
        # Define one end of the axis and the axis itself
        coords_1 = coords[self.atoms[1].index]
        bond_axis = coords[self.atoms[2].index] - coords_1
        # Create the relevant AffineTransform object
        trans_to_origin = transforms.translation(-1.0 * coords_1)
        trans_back = transforms.translation(+1.0 * coords_1)
        rotation = transforms.proper_rotation(bond_axis,
                                              rotation_angle,
                                              affine=True)
        overall = trans_back * rotation * trans_to_origin
        return overall(coords)

    def get_moving_atoms(self):
        """ Returns a list of moving atoms for a given dihedral.

        :return moving_atoms: list of moving atoms.
        """
        self.atom_graph.remove_edge(self.atoms[1], self.atoms[2])
        moving_atoms = [g for g in nx.connected_components(self.atom_graph) if self.atoms[2] in g][0]
        self.moving_atoms = [atom for atom in moving_atoms if not isinstance(atom, chirality.GhostAtom)]
        self.atom_graph.add_edge(self.atoms[1], self.atoms[2])
        return self.moving_atoms


def map_dihedrals(residue_identities, residue_atom_map):
    """ Generates tuples of atoms in each dihedral for every residue.

    :param residue_identities: Dictionary of atom identities.
    :param residue_atom_map: Dictionary of atom mappings from the molecule.atoms graph to corresponding
    residue_sidechains atom names.
    :return dihedral_dict: Dictionary of dihedrals for each residue in molecule.residues. Entries are lists of lists,
    where each element of the sublist is an atom in molecule.residue.atoms.
    """
    from rotamer.topology.residue_sidechains import dihedral_chain

    dihedral_dict = {}
    for residue, res_id in residue_identities.items():
        # If the residue doesn't have an entry in dihedral_chain, skip to the next residue.
        if not dihedral_chain[res_id]:
            # print "Skipping:", residue, res_id
            continue
        unmapped_dihedrals = [dihedral_chain[res_id][i:i + 4] for i in range(0, len(dihedral_chain[res_id]) - 3)]
        mapped_dihedrals = []
        for dihedral in unmapped_dihedrals:
            # We need to convert "N-1" to the N neighbouring "C0". The others can be converted using the map.
            atom_map = residue_atom_map[residue]
            mapped_dihedral = [
                [nb for nb in nx.neighbors(residue.molecule.atoms, atom_map["C0"]) if nb.element == "N"][0]
                if atom_name == "N-1" else atom_map[atom_name] for atom_name in dihedral]
            mapped_dihedrals.append(mapped_dihedral)
        dihedral_dict[residue] = mapped_dihedrals
    return dihedral_dict


if __name__ == "__main__":
    import os.path
    import numpy as np

    # pr = cProfile.Profile()
    # pr.enable()
    topology_data = amber.read_topology(os.path.normpath("/home/khs26/flu.prmtop"))
    coords = np.array(read_amber_restart(os.path.normpath("/home/khs26/flu.inpcrd"))).reshape((-1, 3))
    molecule = amber.create_molecule(topology_data)
    ress, maps = molecule.identify_residues()
    dihedrals = []
    for v in map_dihedrals(ress, maps).values():
        for dihedral in v:
            dihedrals.append(Dihedral(dihedral))
    for dihe in sorted(dihedrals, key=lambda x: x.residue.index):
        if dihe.residue.identity in symmetric_atoms:
            print dihe.residue, dihe.atoms, dihe.residue.atoms
            print "Angle:", dihe.measure_dihedral(coords) * 180.0 / np.pi
            # pr.disable()
            # pr.print_stats('cumulative')