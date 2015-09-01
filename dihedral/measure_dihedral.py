import numpy as np

def dihedral_angle(coords, dihedral):
    """
    Measure dihedral angle for this set of atoms.

    :param coords: Coordinates for the entire molecule
    :param dihedral: List of atoms in the dihedral.
    :return angle: Measured dihedral angle
    """
    atom_coords = [coords[x.index] for x in dihedral]
    b1, b2, b3 = [atom_coords[i+1] - atom_coords[i] for i in range(0, 3)]
    b2_b3 = np.cross(b2, b3)
    b1_b2 = np.cross(b1, b2)
    angle = np.arctan2(np.linalg.norm(b2) * np.dot(b1, b2_b3), np.dot(b1_b2, b2_b3))
    return angle

def restrict_angle_value(angle, symmetry_order):
    """ Restricts an angle between -180 and 180 to be the smallest (closest to zero) value which is equivalent by
    rotational symmetry of order symmetry_order.

    :param angle: Angle to be normalised
    :param symmetry_order: Rotational symmetry
    :return angle: New value of angle
    """
    max_angle = np.pi / symmetry_order
    if angle < -max_angle:
        angle += 2.0 * max_angle
    elif angle > max_angle:
        angle -= 2.0 * max_angle
    return angle

def dihedrals_with_symmetry(coords, residue, residue_ids, dihedrals):
    """
    Measures dihedral angles for a given residue, accounting for symmetry wrt rotation (e.g. carboxylate oxygens or
    terminal methyls in leucine).

    :param coords:
    :param residue:
    :param residue_ids:
    :param dihedrals:
    :return dihedral_values:
    """
    dihedral_values = {}
    dihedral_values[residue] = [(dihedral, dihedral_angle(coords, dihedral)) for dihedral in dihedrals]
    # This compares against the list of residues which have symmetry.
    if residue_ids[residue] in ["ARG", "ASP", "GLU", "LEU", "PHE", "TYR", "VAL"]:
        print "Restricting:", dihedral_values[residue][-1]
        dihedral_values[residue][-1] = (dihedral_values[residue][-1][0], restrict_angle_value(dihedral_values[residue][-1][1], 2))
    else:
        dihedral_values[residue] = [(dihedral, dihedral_angle(coords, dihedral)) for dihedral in dihedrals]
    return dihedral_values


symmetric_atoms = {"ASP": ["O0", "O1"],
                   "GLU": ["O0", "O1"],
                   "ARG": ["N1", "N2"],
                   "LEU": ["C3", "C4"],
                   "PHE": ["C3", "C7"],
                   "TYR": ["C3", "C7"],
                   "VAL": ["C2", "C3"]}

def dihedral_with_symmetry(coords, dihedral):
    """
    Measures dihedral angles for a given dihedral, accounting for symmetry wrt rotation (e.g. carboxylate oxygens or
    terminal methyls in leucine).

    :param coords:
    :param residue:
    :param residue_ids:
    :param dihedrals:
    :return dihedral_values:
    """
    residue = dihedral.residue
    # This compares against the list of residues which have symmetry.
    if residue.identity in symmetric_atoms:
        if any([dihedral.atom_map[sym_atom] in dihedral.atoms for sym_atom in symmetric_atoms[residue.identity]]):
            # print "Restricting:", dihedral.residue, dihedral.atoms
            angle = restrict_angle_value(dihedral_angle(coords, dihedral.atoms), 2)
        else:
            # print "Not restricting:", dihedral.residue, dihedral.atoms
            angle = dihedral_angle(coords, dihedral.atoms)
    else:
        angle = dihedral_angle(coords, dihedral.atoms)
    return angle

