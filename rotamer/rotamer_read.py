import pele.amber.read_amber as ra
import playground.group_rotation.amino_acids as amino
import playground.group_rotation.chirality as chir
import playground.group_rotation.group_rotation as gr
import numpy as np

dirname = '/scratch/khs26/rotamer_lib_igb2/ARG/TRP/ALA/'

# First test the order of atoms in lowest and in the topology file
# with open('/scratch/khs26/rotamer_lib_igb2/ALA/ARG/TRP/lowest', 'r') as lowest:
#     num_atoms = int(lowest.readline())
#     energy = float(lowest.readline().split()[4])
#     coords = []
#     for line in lowest:
#         if len(coords) >= num_atoms:
#             break
#         coords.append(line.split())
#     print num_atoms, energy
#     print coords

# Open and parse the topology file
topology = ra.read_topology(''.join((dirname, 'coords.prmtop')))
mol_graph = ra.create_atoms_and_residues(topology)

# We need to know the residues of interest and their dihedrals
res1 = next(residue for residue in mol_graph.residues.nodes() if residue.index == 1)
res2 = next(residue for residue in mol_graph.residues.nodes() if residue.index == 2)
res3 = next(residue for residue in mol_graph.residues.nodes() if residue.index == 3)

# Get the list of possible bonds from the amino_acids module
dihedral_bonds = {}
for res in (res1, res2, res3):
    dihedral_bonds[res] = [bond[1] for bond in amino.def_parameters.keys() if bond[0] == res.name]
    # print chir.chiral_order(mol_graph.atoms, dihedral_bonds[res][0][0], depth=2)
print dihedral_bonds

# Define all the dihedrals
# Get the atom node in the networkx graph
def get_residue_atom_node(atom_name, residue, graph):
    return next(node for node in graph.nodes() if node.name == atom_name and node.residue == residue)

dihedrals = {}
for res in (res1, res2, res3):
    for atom_pair in dihedral_bonds[res]:
        atom0 = get_residue_atom_node(atom_pair[0], res, mol_graph.atoms)
        atom1 = get_residue_atom_node(atom_pair[1], res, mol_graph.atoms)
        atom_1 = [atom for atom in chir.chiral_order(mol_graph.atoms, atom0, depth=2) if atom != atom1][0]
        atom2 = [atom for atom in chir.chiral_order(mol_graph.atoms, atom1, depth=2) if atom != atom0][0]
        dihedrals[(res, atom0, atom1)] = (atom_1, atom0, atom1, atom2)

# print dihedrals

# Read coords and energies from the lowest file
class LowestFile(object):
    """
    Reads in a lowest file and outputs a dictionary of structure coords, keyed by energy.
    """
    def __init__(self, input_filename):
        self.structures = self.read(input_filename)
        self.structures = self.normalise_structures(self.structures)
    def read(self, input_filename):
        '''
        Read the file...
        '''
        structures = {}
        with open(input_filename, 'r') as input_file:
            num_atoms = 0
            energy = []
            coords = []
            for line in input_file:
            # Look for num_atoms
                try:
                    if coords:
                        structures[energy] = coords
                    num_atoms = int(line)
                    coords = []
                except ValueError:
                # Look for energy
                    try:
                        energy = line.split()[4]
                # Otherwise, read in the coords
                    except IndexError:
                        coords.append(line.split())
        return structures
    def normalise_structures(self, structures):
        '''
        Normalise the energy of a series of structures
        '''
        min_energy = min(map(float, structures.keys()))
        normalised = {}
        for energy in structures:
            norm_energy = float(energy) - min_energy
            normalised[norm_energy] = structures[energy]
        return normalised

myLowestFile = LowestFile(''.join((dirname, 'lowest')))

def to_degrees(rads):
    return 180.0 * rads / np.pi

# Measure the dihedrals
for energy in sorted(myLowestFile.structures, key=lambda x: float(x)):
    coords = np.array([map(float, coord[1:]) for coord in myLowestFile.structures[energy]]).flatten()
    print energy
    for dihedral in sorted(dihedrals, key=lambda x: map(lambda item: item.index, x)):
        angle = gr.measure_dihedral(coords, map(lambda x: x.index, dihedrals[dihedral]))
        if ['O', 'C', 'CA', 'N'] == [di.name for di in dihedrals[dihedral]] and dihedral[0].index == 2:
            print '   ', dihedral[0], ':', dihedrals[dihedral], ':', to_degrees(angle)