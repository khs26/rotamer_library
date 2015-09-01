import networkx as nx

"""
Defines residue sidechains, atom graphs and dihedral definitions for residues.
"""


def add_hydrogens(graph, number_dict):
    """ Adds hydrogens to heavy atoms in a sidechain graph.

    :param graph: sidechain graph
    :param number_dict: dictionary of number of hydrogens per heavy atom
    """
    hydrogens = ["".join(["H", str(i)]) for i in range(0, sum(number_dict.values()))]
    for heavy_atom in number_dict:
        for i in range(0, number_dict[heavy_atom]):
            graph.add_edge(heavy_atom, hydrogens.pop(0))


def update_graph_dict(graph):
    """ Updates the node dictionary for a graph with elements, so that it can be used in isomorphism tests.

    :param graph: sidechain graph
    """
    for node in graph.nodes():
        graph.node[node]["element"] = node.rstrip("0123456789")

# Dihedral chain contains the list of heavy atoms which define the rotamer dihedrals.
dihedral_chain = {}

# A - Alanine
backbone = [("C0", "C1")]
hydrogen_count = {"C1": 3}
alanine = nx.Graph(backbone)
add_hydrogens(alanine, hydrogen_count)
dihedral_chain["ALA"] = None

# C - Cysteine
heavy_atoms = [("C0", "C1"),
               ("C1", "S0")]
hydrogen_count = {"C1": 2,
                  "S0": 1}
cysteine = nx.Graph(heavy_atoms)
add_hydrogens(cysteine, hydrogen_count)
dihedral_chain["CYS"] = ("N-1", "C0", "C1", "S0")

# Cx - Cysteine with salt bridge
# Not needed (can't really compare non-salt bridged rotamer states)

# D - Aspartate
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "O0"),
               ("C2", "O1")]
hydrogen_count = {"C1": 2}
aspartate = nx.Graph(heavy_atoms)
add_hydrogens(aspartate, hydrogen_count)
dihedral_chain["ASP"] = ("N-1", "C0", "C1", "C2", "O0")

# E - Glutamate
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "O0"),
               ("C3", "O1")]
hydrogen_count = {"C1": 2,
                  "C2": 2}
glutamate = nx.Graph(heavy_atoms)
add_hydrogens(glutamate, hydrogen_count)
dihedral_chain["GLU"] = ("N-1", "C0", "C1", "C2", "C3", "O0")

# F - Phenylalanine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "C4"),
               ("C4", "C5"),
               ("C5", "C6"),
               ("C6", "C7"),
               ("C7", "C2")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "C5": 1,
                  "C6": 1,
                  "C7": 1}
phenylalanine = nx.Graph(heavy_atoms)
add_hydrogens(phenylalanine, hydrogen_count)
dihedral_chain["PHE"] = ("N-1", "C0", "C1", "C2", "C3")

# G - Glycine
# Not needed
dihedral_chain["GLY"] = None

# Hd - Histidine (delta protonated)
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "N0"),
               ("N0", "C4"),
               ("C4", "N1"),
               ("N1", "C2")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "N1": 1}
histidine_d = nx.Graph(heavy_atoms)
add_hydrogens(histidine_d, hydrogen_count)
dihedral_chain["HID"] = ("N-1", "C0", "C1", "C2", "N1")

# He - Histidine (epsilon protonated)
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "N0"),
               ("N0", "C4"),
               ("C4", "N1"),
               ("N1", "C2")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "N0": 1}
histidine_e = nx.Graph(heavy_atoms)
add_hydrogens(histidine_e, hydrogen_count)
dihedral_chain["HIE"] = ("N-1", "C0", "C1", "C2", "N1")

# Hp - Histidine (both protonated)
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "N0"),
               ("N0", "C4"),
               ("C4", "N1"),
               ("N1", "C2")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "N0": 1,
                  "N1": 1}
histidine_p = nx.Graph(heavy_atoms)
add_hydrogens(histidine_p, hydrogen_count)
dihedral_chain["HIP"] = ("N-1", "C0", "C1", "C2", "N1")

# I - Isoleucine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C1", "C3"),
               ("C3", "C4")]
hydrogen_count = {"C1": 1,
                  "C2": 3,
                  "C3": 2,
                  "C4": 3}
isoleucine = nx.Graph(heavy_atoms)
add_hydrogens(isoleucine, hydrogen_count)
dihedrals = [("C0", "C1")]
dihedral_chain["ILE"] = ("N-1", "C0", "C1", "C3", "C4")

# K - Lysine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "C4"),
               ("C4", "N0")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "C3": 2,
                  "C4": 2,
                  "N0": 3}
lysine = nx.Graph(heavy_atoms)
add_hydrogens(lysine, hydrogen_count)
dihedral_chain["LYS"] = ("N-1", "C0", "C1", "C2", "C3", "C4", "N0")

# L - Leucine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C2", "C4")]
hydrogen_count = {"C1": 2,
                  "C2": 1,
                  "C3": 3,
                  "C4": 3}
leucine = nx.Graph(heavy_atoms)
add_hydrogens(leucine, hydrogen_count)
dihedral_chain["LEU"] = ("N-1", "C0", "C1", "C2", "C3")

# M - Methionine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "S0"),
               ("S0", "C3")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "C3": 3}
methionine = nx.Graph(heavy_atoms)
add_hydrogens(methionine, hydrogen_count)
dihedral_chain["MET"] = ("N-1", "C0", "C1", "C2", "S0", "C3")

# N - Asparagine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "O0"),
               ("C2", "N0")]
hydrogen_count = {"C1": 2,
                  "N0": 2}
asparagine = nx.Graph(heavy_atoms)
add_hydrogens(asparagine, hydrogen_count)
dihedral_chain["ASN"] = ("N-1", "C0", "C1", "C2", "O0")

# P - Proline
# Not needed

# Q - Glutamine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "O0"),
               ("C3", "N0")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "N0": 2}
glutamine = nx.Graph(heavy_atoms)
add_hydrogens(glutamine, hydrogen_count)
dihedral_chain["GLN"] = ("N-1", "C0", "C1", "C2", "C3", "O0")

# R - Arginine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "N0"),
               ("N0", "C4"),
               ("C4", "N1"),
               ("C4", "N2")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "C3": 2,
                  "N0": 1,
                  "N1": 2,
                  "N2": 2}
arginine = nx.Graph(heavy_atoms)
add_hydrogens(arginine, hydrogen_count)
dihedral_chain["ARG"] = ("N-1", "C0", "C1", "C2", "C3", "N0", "C4", "N1")

# S - Serine
heavy_atoms = [("C0", "C1"),
               ("C1", "O0")]
hydrogen_count = {"C1": 2,
                  "O0": 1}
serine = nx.Graph(heavy_atoms)
add_hydrogens(serine, hydrogen_count)
dihedral_chain["SER"] = ("N-1", "C0", "C1", "O0")

# T - Threonine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C1", "O0")]
hydrogen_count = {"C1": 1,
                  "C2": 3,
                  "O0": 1}
threonine = nx.Graph(heavy_atoms)
add_hydrogens(threonine, hydrogen_count)
dihedral_chain["THR"] = ("N-1", "C0", "C1", "O0")

# V - Valine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C1", "C3")]
hydrogen_count = {"C1": 1,
                  "C2": 3,
                  "C3": 3}
valine = nx.Graph(heavy_atoms)
add_hydrogens(valine, hydrogen_count)
dihedral_chain["VAL"] = ("N-1", "C0", "C1", "C2")

# W - Tryptophan
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "N0"),
               ("N0", "C4"),
               ("C4", "C5"),
               ("C5", "C6"),
               ("C6", "C7"),
               ("C7", "C8"),
               ("C8", "C9"),
               ("C9", "C4"),
               ("C9", "C2")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "N0": 1,
                  "C5": 1,
                  "C6": 1,
                  "C7": 1,
                  "C8": 1}
tryptophan = nx.Graph(heavy_atoms)
add_hydrogens(tryptophan, hydrogen_count)
dihedral_chain["TRP"] = ("N-1", "C0", "C1", "C2", "C3")

# Y - Tyrosine
heavy_atoms = [("C0", "C1"),
               ("C1", "C2"),
               ("C2", "C3"),
               ("C3", "C4"),
               ("C4", "C5"),
               ("C5", "C6"),
               ("C6", "C7"),
               ("C7", "C2"),
               ("C5", "O0")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "O0": 1,
                  "C6": 1,
                  "C7": 1}
tyrosine = nx.Graph(heavy_atoms)
add_hydrogens(tyrosine, hydrogen_count)
dihedral_chain["TYR"] = ("N-1", "C0", "C1", "C2", "C3")

amino_acids = {"ALA": alanine,
               "CYS": cysteine,
               "ASP": aspartate,
               "GLU": glutamate,
               "PHE": phenylalanine,
               "HID": histidine_d,
               "HIE": histidine_e,
               "HIP": histidine_p,
               "ILE": isoleucine,
               "LYS": lysine,
               "LEU": leucine,
               "MET": methionine,
               "ASN": asparagine,
               "GLN": glutamine,
               "ARG": arginine,
               "SER": serine,
               "THR": threonine,
               "VAL": valine,
               "TRP": tryptophan,
               "TYR": tyrosine}

map(update_graph_dict, amino_acids.values())
