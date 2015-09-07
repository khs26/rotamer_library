import residue_sidechains as res_scs
import read_amber_prmtop as amber
import chirality as chir
import networkx as nx


def res_atom_graph(molecule_graph, residues):
    """returns the graph of atoms in residues

    :param molecule: molecule containing residue and atom graphs
    :type molecule: amber.Molecule
    :param residues: residues
    :type residues: amber.Residue
    :return: graph containing only atoms in residues
    :rtype: nx.Graph
    """
    res_atoms = []
    for res in residues:
        res_atoms += res.atoms
    return molecule_graph.atoms.subgraph(res_atoms)


def find_sidechains(molecule_graph):
    # Identify chiral atoms
    atoms = molecule_graph.atoms
    chiral_centres = chir.get_chiral_sets(atoms)
    # Identify sidechains (Ca-Cb-X), apart from proline and glycine.
    sidechains = {}
    for k, v in chiral_centres.items():
        carbons = [atom for atom in v if atom.element == 'C']
        amides = [carbon for carbon in carbons
                  if any([type(nb) == chir.GhostAtom and nb.element == 'O' for nb in nx.neighbors(atoms, carbon)])
                  and any([nb.element == 'N' or nb.element == 'O' for nb in nx.neighbors(atoms, carbon)])]
        nbs_n = [nb for nb in v if nb.element == 'N']
        if amides and nbs_n:
            amide_bond = (k, amides[0])
            n_bond = (k, nbs_n[0])
            h_bond = (k, [h for h in nx.neighbors(atoms, k) if h.element == 'H'][0])
            # Now find sidechains by cutting the Ca-C, Ca-N and Ca-H bonds
            atoms.remove_edges_from([amide_bond, n_bond, h_bond])
            sidechain_atoms = [atom for atom in [comp for comp in nx.connected_components(atoms) if k in comp][0]
                               if type(atom) != chir.GhostAtom]
            atoms.add_edges_from([amide_bond, n_bond, h_bond])
            if not any([k in cycle for cycle in nx.cycle_basis(atoms.subgraph(sidechain_atoms))]):
                sidechains[k] = atoms.subgraph(sidechain_atoms)
    return sidechains


def residue_from_sidechain(sidechains):
    """ For each sidechain in sidechains, determines the type of residue and generates the mapping from the atom types
    in residue_sidechains to atoms in the molecule.

    :param sidechains: {Ca: sidechain atoms graph} dictionary
    :return: residues: residue type
    :return: mapping:  atom mapping dictionary
    """
    residues = {}
    mapping = {}
    for sc in sidechains:
        for res in res_scs.amino_acids:
            graph_match = nx.algorithms.isomorphism.GraphMatcher(res_scs.amino_acids[res], sidechains[sc],
                                                                 node_match=(lambda x, y: x['element'] == y['element']))
            if graph_match.is_isomorphic():
                residues[sc.residue] = res
                mapping[sc.residue] = graph_match.mapping
                break
    return residues, mapping


if __name__ == "__main__":
    import os.path

    topology_data = amber.read_topology(os.path.normpath("/home/khs26/flu.prmtop"))
    molecule = amber.create_molecule(topology_data)
    cands = chir.tetravalent_atoms(molecule.atoms)
    chir.multi_bonds(molecule.atoms)
    cands2 = chir.rankable_neighbours(cands)[0]
    print len(molecule.atoms), len(cands), len(cands2)
    # scs = find_sidechains(molecule, [res for res in molecule.residues.nodes()])
    # ress, maps = residue_from_sidechain(scs)
    # for k, v in sorted(ress.items()):
    # print k, v
    #      for i, j in res_scs.dihedrals:
    #         print i, j, chir.chiral_order(molecule.atoms, maps[k][i]), chir.chiral_order(molecule.atoms, maps[k][j])