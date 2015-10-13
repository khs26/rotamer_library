from read_amber_prmtop import parse_topology_file
import networkx as nx

def find_backbone_carbons(molecule):
    atoms = molecule.atoms
    carbons = [atom for atom in atoms if atom.element == 'C']
    oxygens = [atom for atom in atoms if atom.element == 'O']
    nitrogens = [atom for atom in atoms if atom.element == 'N']
    for c in carbons:
        print c, nx.neighbors(atoms, c)
    c_carbonyls = [atom for atom in carbons
                   if atoms.degree(atom) == 3
                   and list(set(nx.neighbors(atoms, atom)) & set(oxygens))]
    c_alphas = [atom for atom in carbons
                if atoms.degree(atom) == 4
                and list(set(nx.neighbors(atoms, atom)) & set(c_carbonyls))
                and list(set(nx.neighbors(atoms, atom)) & set(nitrogens))]
    backbone_atoms = {}
    for c_alpha in c_alphas:
        this_res = c_alpha.residue
        backbone_atoms[this_res] = (c_alpha, list(set(nx.neighbors(atoms, c_alpha)) & set(c_carbonyls))[0])
        this_res.backbone_map = {}
        this_res.backbone_map["CA"] = c_alpha
        this_res.backbone_map["C'"] = backbone_atoms[c_alpha.residue][1]
        this_res.backbone_map["N"] = list(set(nx.neighbors(atoms, c_alpha)) & set(nitrogens))[0]
    return backbone_atoms

if __name__ == "__main__":
    molecule = parse_topology_file("/home/khs26/flu.prmtop")