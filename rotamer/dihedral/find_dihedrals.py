import networkx as nx
import rotamer.topology.read_amber_prmtop as amber
import rotamer.topology.chirality as chirality
import transforms
from measure_dihedral import dihedral_with_symmetry, symmetric_atoms
from rotamer.io.amber import read_amber_restart
from rotamer.topology.identify_backbone import find_backbone_carbons


class Dihedral(object):
    def __init__(self, dihedral_atoms):
        self.atoms = dihedral_atoms
        # Use atom with index 1, because it should be inside the residue (e.g. phi/psi angles).
        self.residue = dihedral_atoms[1].residue
        self.atom_graph = dihedral_atoms[1].molecule.atoms
        # self.atom_map = {v: k for k, v in self.residue.atom_map.items()}
        # self.atom_map = self.residue.atom_map
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
        # Define the axis
        axis = np.vstack([coords[self.atoms[1].index], coords[self.atoms[2].index]])
        # Create the transform object
        dihedral_rotation = transforms.general_rotation(axis, rotation_angle, affine=True)
        # Get the moving coordinates and make sure to only apply the transform to those.
        moving_indices = [atom.index for atom in self.get_moving_atoms()]
        moved_coords = dihedral_rotation(coords[moving_indices]).reshape(-1, 3)
        # Since we have used advanced indexing, we have a copy of the values and need to copy them to the new_coords
        # array.
        new_coords = coords.copy()
        for i, moved_coord in enumerate(moved_coords):
            new_coords[moving_indices[i]] = moved_coord
        return new_coords

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

def phi_psi_dihedrals(molecule):
    backbone_dict = find_backbone_carbons(molecule)
    # Chain should be C'(-1) - N(0) - Ca(0) - C'(0) - N(1)
    phis, psis = [], []
    for res in sorted(backbone_dict.keys()):
        c_alpha = res.backbone_map["CA"]
        c_prime = res.backbone_map["C'"]
        nitrogen = res.backbone_map["N"]
        try:
            c_minus_1 = [nb for nb in nx.neighbors(molecule.atoms, nitrogen)
                         if any([nb2.element == 'O' for nb2 in nx.neighbors(molecule.atoms, nb)])][0]
            phi = Dihedral([c_minus_1, nitrogen, c_alpha, c_prime])
            res.dihedrals["phi"] = phi
            phis.append(phi)
        except IndexError:
            c_minus_1 = None
        try:
            n_plus_1 = [nb for nb in nx.neighbors(molecule.atoms, c_prime) if nb.element == 'N'][0]
            psi = Dihedral([nitrogen, c_alpha, c_prime, n_plus_1])
            res.dihedrals["psi"] = psi
            psis.append(psi)
        except IndexError:
            n_plus_1 = None
    return phis, psis

def sidechain_dihedrals(molecule):
    ress, maps = molecule.identify_residues()
    sc_dihedrals = []
    for k, v in map_dihedrals(ress, maps).items():
        for idx, dihedral in enumerate(v):
            k.dihedrals["chi" + str(idx+1)] = Dihedral(dihedral)
            sc_dihedrals.append(dihedral)
    return sc_dihedrals


if __name__ == "__main__":
    import os.path
    import numpy as np

    topology_data = amber.read_topology(os.path.normpath("../tests/data/ARG_LYS_ASN.prmtop"))
    coords = np.array(read_amber_restart(os.path.normpath("../tests/data/ARG_LYS_ASN.inpcrd"))).reshape((-1, 3))
    molecule = amber.create_molecule(topology_data)
    phi_psi_dihedrals(molecule)
    sidechain_dihedrals(molecule)
    for res in sorted(molecule.residues, key=lambda x: x.index)[1:4]:
        print "========================"
        print res
        for k, v in res.dihedrals.items():
            # print k, v.measure_dihedral(coords) * 180.0 / np.pi
            new_coords = v.set_dihedral(coords, -45.0 * np.pi / 180.0)
            print k, v.measure_dihedral(coords) * 180.0 / np.pi, v.measure_dihedral(new_coords) * 180.0 / np.pi