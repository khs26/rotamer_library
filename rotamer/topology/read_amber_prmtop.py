import operator as op
import exceptions as ex
import itertools
import re

import networkx as nx
import numpy as np

import elements as elem
from rotamer.io.amber import read_amber_restart


class Atom(object):
    """ Atom defined from the AMBER topology file. """

    def __init__(self, index, name, mass, amber_atom_type, charge, molecule, **kwargs):
        self.index = int(index)
        self.name = name.strip()
        self.mass = float(mass)
        self.amber_atom_type = amber_atom_type
        self.charge = float(charge)
        self.residue = None
        self.set_element(mass)
        self.molecule = molecule
        self.coords = []

    def set_element(self, mass):
        """ Defines the element, based on the mass. """
        self.element = elem.lookup_element_by_mass(mass)

    def set_residue(self, residue):
        """ Defines the Residue to which the Atom belongs. """
        self.residue = residue

    def __repr__(self):
        return str(self.index) + " " + self.element + " " + self.name

        # def __cmp__(self, other):
        # return cmp(self.mass, other.mass)


class Residue(object):
    """ Residue defined from the AMBER topology file. """

    def __init__(self, index, name, molecule):
        self.index = int(index)
        self.name = name.strip()
        self.molecule = molecule
        self.identity = None

    def add_atoms(self, atoms):
        """ 
        Adds Atoms in atoms to the Residue and sets the residue for those
        Atoms accordingly.
        """
        self.atoms = atoms
        for atom in atoms:
            atom.set_residue(self)

    def __repr__(self):
        return str(self.index) + " " + self.name

    def get_identity(self):
        if not self.identity:
            try:
                self.identity = self.molecule.res_identities[self]
                self.atom_map = self.molecule.atom_map[self]
            except AttributeError:
                ident, atom_map = self.molecule.identify_residues()
                self.identity, self.atom_map = ident[self], atom_map[self]
        return self.identity


class Molecule(object):
    """ Molecule defined from the AMBER topology file. """

    def __init__(self):
        self.residues = nx.Graph()
        self.atoms = nx.Graph()

    def read_coords(self, inpcrd_filename):
        """ Reads coordinates from a file into the Atoms of the Molecule. """
        coords = np.array(read_amber_restart(inpcrd_filename)).reshape((-1, 3))
        # Sanity check that there are 3 coordinates per atom
        if len(coords) != self.num_atoms():
            raise ex.RuntimeError("Incorrect number of coordinates read into molecule.")
        # Now run through the list of atoms (sorted by index) and coordinates
        for atom, xyz in itertools.izip(sorted(self.atoms.nodes(), key=op.attrgetter("index")), coords):
            atom.coords = xyz

    def num_atoms(self):
        return len(self.atoms.nodes())

    def identify_residues(self):
        import identify_residue as id_res

        sidechains = id_res.find_sidechains(self)
        self.res_identities, self.atom_maps = id_res.residue_from_sidechain(sidechains)
        for res in self.res_identities:
            res.identity = self.res_identities[res]
            res.atom_map = self.atom_maps[res]
        return self.res_identities, self.atom_maps


def split_len(seq, length):
    """ Returns a list containing the elements of seq, split into length-long chunks. """
    return [seq[i:i + length] for i in range(0, len(seq), length)]


def data_cast(datum, type):
    """
    Casts the datum as defined by type.

    a = string
    i = integer
    e = float (in exponent notation)
    """
    if type == "a" or type == "A":
        return str(datum)
    if type == "i" or type == "I":
        return int(datum)
    if type == "e" or type == "E":
        return float(datum)


def read_topology(filename):
    """ Reads a topology file, filename, and returns data blocks containing the topology data. """
    # Read the topology file to a series of lines. 
    # prmtop_input_lines = [[]]
    with open(filename, "r") as topology_file:
        prmtop_input_lines = []
        for line in topology_file:
            line = line.strip("\n")
            if line.startswith("%VERSION"):
                continue
            if line.startswith("%FLAG"):
                prmtop_input_lines.append([line])
            else:
                prmtop_input_lines[-1].append(line)
    # Parse the input into data blocks, a list containing dictionaries with the following entries.
    # "flag": flag describing the data block (e.g. atom number, residue type);
    # "format": describes the format of the data block (e.g. 8a10 is 8 data of type char, length 10);
    # "data": contains the data, split into blocks of appropriate length and formatted to the correct type.
    prmtop_data_blocks = {}
    for input_line in prmtop_input_lines:
        # Get rid of %FLAG from front and strip whitespace.
        flag = input_line[0][6:].rstrip(" ")
        # Find length and data type.
        data_length = int(re.split("[aEI\.\)]", input_line[1])[1])
        data_type = re.findall("([a-zA-Z])", input_line[1][8:])[0]
        # Split the string into chunks of the appropriate size and cast to appropriate data type.
        split_string = split_len("".join(input_line[2:]), data_length)
        data = map(data_cast, split_string, itertools.repeat(data_type, len(split_string)))
        prmtop_data_blocks[flag] = data
    return prmtop_data_blocks


def create_molecule(topology_data):
    """
    This function takes a topology data set (topology_data) and returns a Molecule object containing
    a graph of the atoms and residues and their connectivity.
    """
    molecule = Molecule()
    # Create a list of atoms from the topology data and add them to the molecule's graph. 
    # The first element of POINTERS is the atom count.
    atoms = list(itertools.imap(Atom,
                                range(0, (topology_data["POINTERS"][0])),
                                topology_data["ATOM_NAME"],
                                topology_data["MASS"],
                                topology_data["AMBER_ATOM_TYPE"],
                                topology_data["CHARGE"],
                                itertools.repeat(molecule)))
    atom_list = [(atom, atom.__dict__) for atom in atoms]
    molecule.atoms.add_nodes_from(atom_list)
    # Create a list of residues and add them to the molecule's list of residues.
    # The 11th element of POINTERS is the residue count.
    residues = list(itertools.imap(Residue,
                                   range(0, (topology_data["POINTERS"][11])),
                                   topology_data["RESIDUE_LABEL"],
                                   itertools.repeat(molecule)))
    molecule.residues.add_nodes_from(residues)
    # Go through the BONDS_INC_HYDROGEN and BONDS_WITHOUT_HYDROGEN lists to extract lists of bonded atoms.
    # Index is i / 3 + 1, because AMBER still uses coordinate indices for runtime speed.
    first_atoms = map(lambda x: (x / 3) + 1,
                      topology_data["BONDS_INC_HYDROGEN"][0::3] + topology_data["BONDS_WITHOUT_HYDROGEN"][0::3])
    second_atoms = map(lambda x: (x / 3) + 1,
                       topology_data["BONDS_INC_HYDROGEN"][1::3] + topology_data["BONDS_WITHOUT_HYDROGEN"][1::3])
    bond_list = zip(first_atoms, second_atoms)
    # Next put the appropriate atoms into the appropriate residues.  AMBER specifies the first atom index of
    # each residue (i.e. no end, hence the exception).
    residue_indices = topology_data["RESIDUE_POINTER"]
    for i, residue in enumerate(residues):
        try:
            start = residue_indices[i] - 1
            end = residue_indices[i + 1] - 1
        except IndexError, ex:
            end = None
        residue.add_atoms(atoms[start:end])
    # Now go through bond list and create bonds between the relevant atoms.
    for bond in bond_list:
        bonded_atoms = (atoms[(bond[0] - 1)], atoms[(bond[1] - 1)])
        molecule.atoms.add_edge(*bonded_atoms)
        if bonded_atoms[0].residue != bonded_atoms[1].residue:
            molecule.residues.add_edge(bonded_atoms[0].residue, bonded_atoms[1].residue)
    return molecule


def group_rotation_file(molecule, params, filename):
    """
    Defines group rotation for bonds between atoms defined in params.  This function DOESN'T check
    whether the atoms are actually bonded at any point (so make sure you're referencing appropriate
    atom pairs.

    atoms and residues are the atoms and residues lists returned from create_molecule()

    params has the format:

    params[(residue, (atom_x, atom_y))] = (rotation_scale, probability)
    e.g. params[("ARG", ("CA", "CB"))] = (0.5, 0.2) rotates about the ARG CA-CB bond at most 90 degs each 
         time, with probability 0.2.
    """
    # Loop through the parameters, creating GROUPROTATION information for any bonds described by an individual
    # parameter.
    with open(filename, "w") as output:
        for param in params.keys():
            for residue in molecule.residues.nodes():
                if residue.name == param[0]:
                    rotated_atoms = get_rotated_atoms((residue, param[1]))
                    # Create an identifiable group name consisting of [index]_[res_name]_[atom_1]_[atom_2]
                    group_name = "_".join((str(residue.index),
                                           residue.name.strip(),
                                           rotated_atoms[0].name.strip(),
                                           rotated_atoms[1].name.strip()))
                    output.write(" ".join(("GROUP",
                                           group_name,
                                           str(rotated_atoms[0].index + 1),
                                           str(rotated_atoms[1].index + 1),
                                           str(len(rotated_atoms[2])),
                                           str(params[param][0]),
                                           str(params[param][1]),
                                           "\n")))
                    for atom in rotated_atoms[2]:
                        output.write(str(atom.index + 1) + "\n")
                    output.write("\n")


def group_rotation_dict(molecule, params):
    groups = {}
    for param in params.keys():
        for residue in molecule.residues.nodes():
            if residue.name == param[0]:
                rotated_atoms = get_rotated_atoms((residue, param[1]))
                # Create an identifiable group name consisting of [index]_[res_name]_[atom_1]_[atom_2]
                group_name = "_".join((str(residue.index),
                                       residue.name.strip(),
                                       rotated_atoms[0].name.strip(),
                                       rotated_atoms[1].name.strip()))
                groups[group_name] = {}
                groups[group_name]["bond_atom_1"] = rotated_atoms[0].index
                groups[group_name]["bond_atom_2"] = rotated_atoms[1].index
                groups[group_name]["group_atoms"] = [atom.index for atom in rotated_atoms[2]]
                groups[group_name]["max_angle_magnitude"] = params[param][0]
                groups[group_name]["selection_probability"] = params[param][1]
    return groups


def get_rotated_atoms(bond):
    """
    Returns a list of the atoms in the rotating group for the GROUPROTATION script.
    
    bond has the format:

    bond = (residue, (atom_x, atom_y))
    """
    residue = bond[0]
    for residue_atom in residue.atoms:
        if residue_atom.name == bond[1][0]:
            atom_a = residue_atom
        if residue_atom.name == bond[1][1]:
            atom_b = residue_atom
    # Remove the edge, find the smallest subgraph and then add the edge back
    residue.molecule.atoms.remove_edge(atom_a, atom_b)
    rotating_atoms = sorted(list(nx.connected_components(residue.molecule.atoms))[-1])
    residue.molecule.atoms.add_edge(atom_a, atom_b)
    if atom_a in rotating_atoms:
        atom_1 = atom_b
        atom_2 = atom_a
        rotating_atoms.remove(atom_a)
    elif atom_b in rotating_atoms:
        atom_1 = atom_a
        atom_2 = atom_b
        rotating_atoms.remove(atom_b)
    return atom_1, atom_2, rotating_atoms


def parse_topology_file(topology_filename):
    topology_data = read_topology(topology_filename)
    mol = create_molecule(topology_data)
    return mol


def default_parameters(topology_filename):
    # should not import something from playground.  should this function
    # be moved to playground?
    import playground.group_rotation.amino_acids as amino

    topology_data = read_topology(topology_filename)
    parsed = create_molecule(topology_data)
    return group_rotation_dict(parsed, amino.def_parameters)


if __name__ == "__main__":
    import chirality

    topology_data = read_topology("/home/khs26/flu.prmtop")
    mol = create_molecule(topology_data)
    mol.read_coords("/home/khs26/flu.inpcrd")

    # for atom in mol.atoms:
    # if atom.name == 'CA' and atom.residue.name == 'ILE':
    #         ca = atom