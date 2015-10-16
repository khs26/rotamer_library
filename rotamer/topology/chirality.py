#!/usr/bin/python

import read_amber_prmtop as ra
import networkx as nx
import itertools


class GhostAtom(ra.Atom):
    pass


valences = {"H": 1,
            "C": 4,
            "N": 3,
            "O": 2,
            "S": 2,
            "Se": 2,
            "P": 5}


def tetravalent_atoms(atoms):
    """
    Identifies possible candidates (those with 4 bonds)

    :param atoms: Graph of atoms
    :return: List of atoms with 4 bonds
    """
    candidates = [atom for atom in atoms.nodes() if len(nx.edges(atoms, atom)) == 4]
    return candidates


def multi_bonds(atoms):
    multibonded = [atom for atom in atoms.nodes()
                   if len(nx.edges(atoms, atom)) < valences[atom.element]]
    for i, atom in enumerate(multibonded):
        paired = False
        for other in atoms.neighbors(atom):
            if isinstance(other, GhostAtom):
                paired = True
                continue
            if len(nx.edges(atoms, other)) < valences[other.element]:
                ghost_atom = GhostAtom(**(atom.__dict__))
                ghost_atom.name = atom.name + "*"
                ghost_other = GhostAtom(**(other.__dict__))
                ghost_other.name = other.name + "*"
                atoms.add_edge(other, ghost_atom)
                atoms.add_edge(atom, ghost_other)
                paired = True

def remove_ghost_atoms(atoms):
    """
    Removes ghost atoms from the atom graph (not needed after chirality checks).

    :param atoms: Atom graph
    """
    ghost_atoms = [atom for atom in atoms.nodes() if isinstance(atom, GhostAtom)]
    atoms.remove_nodes_from(ghost_atoms)

def rankable_neighbours(chiral_cands):
    """ Checks if the chiral atom candidates have rankable substituents on each site (i.e. discounting those whose
     neighbour list contains the same univalent atoms).

    :param chiral_cands: Atoms to test.
    :return: maybe_chiral, not_chiral: lists of possibly chiral and achiral atoms
    """

    maybe_chiral, not_chiral = [], []
    for chiral_cand in chiral_cands:
        atoms = chiral_cand.molecule.atoms
        neighbours = atoms.neighbors(chiral_cand)
        # Univalent atoms only have the original chiral_cand atom in their neighbour list. Possibly twice, because of
        # the multi-bond routine.
        univalent = [nb for nb in neighbours if all([nb2 == chiral_cand for nb2 in atoms.neighbors(nb)])]
        if len(univalent) > 1 and any([x.mass == y.mass for x, y in itertools.combinations(univalent, 2)]):
            not_chiral.append(chiral_cand)
        else:
            maybe_chiral.append(chiral_cand)
    return maybe_chiral, not_chiral


def chiral_order(atoms, chiral_atom, depth=6):
    # print "\n\nResidue:", chiral_atom.residue, "atom:", chiral_atom
    # print "Neighbours:", atoms.neighbors(chiral_atom)
    # Create a list of ordered atoms to be passed back
    ordered = []
    # Do a quick check whether there are multiple hydrogens
    neighbors = atoms.neighbors(chiral_atom)
    hydrogens = [atom for atom in neighbors if atom.element == "H"]
    if len(hydrogens) < 2:
        tree = nx.bfs_tree(atoms, chiral_atom)
        # Generate the list of shortest paths in the molecule, neglecting the trivial path [chiral_atom]
        paths = sorted(nx.single_source_shortest_path(tree, chiral_atom, depth).values(), reverse=True,
                       key=lambda x: map(lambda at: at.mass, x))[:-1]
        while paths:
            # Pop the first element (highest priority path) from the list of paths and remove any duplicates.
            path = paths.pop(0)
            # print "Path considered:", path
            paths_no_dups = [unpruned for unpruned in paths if unpruned != path]
            # print "Paths:", paths
            # print "Paths without dups:", paths_no_dups
            # If there are any duplicates, the paths list will be smaller and we can't resolve a highest priority
            if len(paths_no_dups) != len(paths):
                paths = paths_no_dups
            # Otherwise, the path is higher priority than all the other paths, so its second atom is the neighbour with
            # highest priority.
            else:
                # print "Best path:", path
                ranked_atom = path[1]
                # print "Ranked atom:", ranked_atom
                ordered.append(ranked_atom)
                # Drop all the paths containing our ranked atom.
                paths = [unpruned for unpruned in paths if unpruned[1] is not ranked_atom]
    else:
        ordered = []
        # ordered = [atom for atom in neighbors if atom.element != "H"]
        # ordered += [atom for atom in neighbors if atom.element == "H"]
    return ordered


def get_chiral_sets(atoms):
    """
    Driver routine for all the chirality stuff.

    :param atoms: Atom graph
    :return: Dictionary of chiral centres and CIP-ordered neighbours
    """
    chiral_cands = tetravalent_atoms(atoms)
    chiral_cands = rankable_neighbours(chiral_cands)[0]
    multi_bonds(atoms)
    chiral_centres = {}
    for i, chiral_atom in enumerate(chiral_cands):
        ordered = chiral_order(atoms, chiral_atom)
        if len(ordered) == 4:
            chiral_centres[chiral_atom] = ordered
    remove_ghost_atoms(atoms)
    return chiral_centres


def get_chiral_atoms(atoms):
    return get_chiral_sets(atoms).keys()


def write_chirality_file(input_filename, output_filename):
    molecule = ra.parse_topology_file(input_filename)
    atoms = molecule.atoms
    chiral_centres = get_chiral_sets(atoms)
    with open(output_filename, "w") as output_file:
        for atom in sorted(chiral_centres.keys(), cmp=lambda x, y: cmp(x.index, y.index)):
            # Write out the list of chiral atoms and their CIP-ranked neighbours.
            output_string = "{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4:>8d}\n".format(atom.index + 1,
                                                                           *[other_atom.index + 1 for other_atom in
                                                                             chiral_centres[atom]])
            output_file.write(output_string)

def calculate_chirality(coords, chiral_centres):
    import numpy as np

    # For centre atom C and atoms ordered I, J, K and L
    # Calculate dihedral of I-C-L-J
    for atom_list in chiral_centres:
        b1 = coords[atom_list[0]] - coords[atom_list[1]]
        b2 = coords[atom_list[4]] - coords[atom_list[0]]
        b3 = coords[atom_list[2]] - coords[atom_list[4]]
        b1xb2 = np.cross(b1, b2)
        b2xb3 = np.cross(b2, b3)
        b1xb2_x_b2xb3 = np.cross(b1xb2, b2xb3)
        b2_norm = b2 / np.linalg.norm(b2)
        angle = np.arctan2(np.dot(b1xb2_x_b2xb3, b2_norm), np.dot(b1xb2, b2xb3))
        print angle

if __name__ == "__main__":
    import rotamer.io.amber
    molecule = ra.parse_topology_file("../library/coords.prmtop")
    atoms = molecule.atoms
    chiral_centres = get_chiral_sets(atoms)
    chiral_centres_list = [[k.index] + [val.index for val in v] for k, v in chiral_centres.items()]
    coords = rotamer.io.amber.read_amber_restart("../library/coords.inpcrd")
    calculate_chirality(coords.reshape((-1, 3)), chiral_centres_list)
