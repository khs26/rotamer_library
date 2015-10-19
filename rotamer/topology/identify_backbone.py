import networkx as nx

from read_amber_prmtop import parse_topology_file


def find_backbone_carbons(molecule):
    atoms = molecule.atoms
    carbons = [atom for atom in atoms if atom.element == 'C']
    oxygens = [atom for atom in atoms if atom.element == 'O']
    nitrogens = [atom for atom in atoms if atom.element == 'N']
    # Carbons with three bonds, with at least one neighbouring oxygen.
    c_carbonyls = [atom for atom in carbons
                   if atoms.degree(atom) == 3
                   and list(set(nx.neighbors(atoms, atom)) & set(oxygens))]
    # Carbons with four bonds, with at least one neighbouring carbonyl and one neighbouring nitrogen.
    c_alphas = [atom for atom in carbons
                if atoms.degree(atom) == 4
                and list(set(nx.neighbors(atoms, atom)) & set(c_carbonyls))
                and list(set(nx.neighbors(atoms, atom)) & set(nitrogens))]
    backbone_atoms = {}
    for c_alpha in c_alphas:
        this_res = c_alpha.residue
        this_res.backbone_map = {}
        this_res.backbone_map["CA"] = c_alpha
        # Carbonyl neighbouring C-alpha
        this_res.backbone_map["C'"] = list(set(nx.neighbors(atoms, c_alpha)) & set(c_carbonyls))[0]
        # Nitrogen neighbouring C-alpha
        this_res.backbone_map["N"] = list(set(nx.neighbors(atoms, c_alpha)) & set(nitrogens))[0]
        backbone_atoms[this_res] = this_res.backbone_map.values()
    return backbone_atoms


if __name__ == "__main__":
    molecule = parse_topology_file("../library/coords.prmtop")
    find_backbone_carbons(molecule)
    for res in sorted(molecule.residues.nodes(), key=lambda x: x.index):
        print res, res.backbone_map.items()