import numpy as np
import tables as ts
import playground.group_rotation.amino_acids as amino
import pele.amber.read_amber as ra
import playground.group_rotation.chirality as chir
import networkx as nx

class RotamerGroupTemplate(ts.IsDescription):
    """
    A compound data type for interoperating with ROTAMER_GROUP_TEMPLATE in GMIN.

    Fortran type definition:
        TYPE ROTAMER_GROUP_TEMPLATE
            CHAR(LEN=16)                :: GROUP_NAME
            CHAR(LEN=16)                :: RES_NAME
            CHAR(LEN=16)                :: DIHEDRAL_ATOM_NAMES(4)
            CHAR(LEN=16), ALLOCATABLE   :: MOVING_ATOM_NAMES(:)
        END TYPE ROTAMER_GROUP_TEMPLATE

    These are designed to be read from HDF5 files, and so this implementation uses
    the appropriate types from PyTables.

    As we do not know the length of the moving atoms array a priori, we instead
    use a CArray (HDF5 compressible dataset) to store moving atom names for each
    group. The individual row corresponding to a group then contains the name of the
    CArray with the relevant atom names in it.
    """
    group_name = ts.StringCol(itemsize=16)
    res_name = ts.StringCol(itemsize=16)
    dihedral_atom_names = ts.StringCol(itemsize=16, shape=(4))
    moving_atoms_carray = ts.StringCol(itemsize=24)

# Open a file in "w"rite mode
fileh = ts.open_file("amber_rotamer_groups.h5", mode="w")

# Get the HDF5 root group
root = fileh.root

# Create the groups for the templates themselves and one for storing the moving atoms arrays.
for groupname in ("RotamerGroupTemplates", "MovingAtomsArrays"):
    group = fileh.create_group(root, groupname)

# Create a filter, telling it to compress with zlib.
filters = ts.Filters(complib='zlib')

for amino_acid in (amino.amino_acids):
    table = fileh.create_table("/RotamerGroupTemplates",
                               amino_acid,
                               RotamerGroupTemplate,
                               "Template for {res}".format(res=amino_acid))

    # Get the record object associated with the table.
    group_template = table.row

    # Read in an appropriate topology file and create a molecular graph.
    filename = '/scratch/khs26/rotamer_lib_igb2/{res}/{res}/{res}/coords.prmtop'.format(res=amino_acid)
    topology = ra.read_topology(filename)
    mol_graph = ra.create_atoms_and_residues(topology)

    # Get the residue name for the first residue.
    res = next(residue for residue in mol_graph.residues.nodes() if residue.index == 1)

    # Get a list of dihedrals we are interested in for this residue.
    dihedrals = sorted([k[1] for k in amino.def_parameters if k[0] == amino_acid
                                                              and not ('C' in k[1] and 'CA' in k[1])])

    # For each pair of atoms in a dihedral, find their highest-ranked neighbours for defining the dihedral angle.
    dihedral_atoms = {}
    dihedral_moving_atoms = {}
    for atom_pair in dihedrals:
        atom0 = next(n for n in mol_graph.atoms.nodes() if n.name == atom_pair[0] and n.residue == res)
        atom1 = next(n for n in mol_graph.atoms.nodes() if n.name == atom_pair[1] and n.residue == res)
        atom_1 = next(atom for atom in chir.chiral_order(mol_graph.atoms, atom0, depth=2) if atom != atom1)
        atom2 = next(atom for atom in chir.chiral_order(mol_graph.atoms, atom1, depth=2) if atom != atom0)
        dihedral_atoms[(atom0.name, atom1.name)] = (atom_1, atom0, atom1, atom2)
        # Now find the moving atoms by breaking the dihedral bond and choosing the subgraph containing atom1.
        mol_graph.atoms.remove_edge(atom0, atom1)
        dihedral_moving_atoms[(atom0.name, atom1.name)] = nx.node_connected_component(mol_graph.atoms, atom1)
        mol_graph.atoms.add_edge(atom0, atom1)


    # Loop through the possible dihedral atom pairs for the amino acid.
    # i is going to form part of the CArray name
    for i, dihedral in enumerate(dihedrals):
        moving_atom_names = [atom.name for atom in dihedral_moving_atoms[dihedral]]
        carray_name = 'carray_{res}_{ind}'.format(res=amino_acid, ind=str(i))
        print amino_acid, i, dihedral, moving_atom_names, carray_name
        ca = fileh.create_carray(root.MovingAtomsArrays,
                                 name=carray_name,
                                 atom=ts.StringAtom(16),
                                 shape=(len(moving_atom_names),),
                                 filters=filters)
        ca[0:] = moving_atom_names

        group_template['group_name'] = "{res}_{dih0}_{dih1}".format(res=amino_acid, dih0=dihedral[0], dih1=dihedral[1])
        group_template['res_name'] = "{res}".format(res=amino_acid)
        group_template['dihedral_atom_names'] = np.array([x.name for x in dihedral_atoms[dihedral]])
        group_template['moving_atoms_carray'] = carray_name

        # Append this element to the row and move on.
        group_template.append()

    # Flush the table buffers
    table.flush()

# Read the records from table "/RotamerGroupTemplates/ARG" and select some
table = root.RotamerGroupTemplates.ARG
e = [(p['group_name'], p['res_name'], p['dihedral_atom_names'], p['moving_atoms_carray']) for p in table]
for elem in e:
    print("Selected values ==>", elem)
    print("Carray:", root.MovingAtomsArrays._v_children[elem[-1]][:])
print("Total selected records ==> ", len(e))

# Finally, close the file (this also will flush all the remaining buffers!)
fileh.close()