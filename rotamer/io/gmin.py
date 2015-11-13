import numpy as np

def read_xyz(xyz_filename, natoms):
    """
    Read coordinates in from a xyz file (possibly containing multiple structures).

    :param xyz_filename: filename to be read from
    :param natoms: number of atoms in the structure
    :return: 3d array of coordinates with shape (n_structures, n_atoms, 3)
    """
    all_coords = np.fromfile(xyz_filename, sep=" ")
    return np.reshape(all_coords, (-1, natoms, 3))

def read_one_structure(lowest_file):
    """
    Read a single structure from a lowest file.

    :param lowest_file: File object for the lowest file.
    :return: Dictionary containing the structure's attributes (e.g. coords and energy).
    """
    natoms = int(lowest_file.readline().strip())
    data_line = lowest_file.readline().split()
    index = int(data_line[3][:-1])
    energy = float(data_line[4])
    first_found = int(data_line[9])
    coords = np.zeros([natoms, 3])
    for i in range(0, natoms):
        coords[i] = map(float, lowest_file.readline().split()[1:])
    return {"index": index, "energy": energy, "first found": first_found, "coords": coords}


def read_lowest(lowest_filename):
    """
    Read all the structures from a lowest file.

    :param lowest_filename: File name of the lowest file to read.
    :return: List of dictionaries, each dictionary corresponding to a single structure from the lowest
             file.
    """
    structures = []
    with open(lowest_filename, "r") as lowest_file:
        while True:
            try:
                structures.append(read_one_structure(lowest_file))
            except ValueError:
                # We raise a ValueError if we go too far. Not sure why this indicates end-of-file.
                break
    return structures


def write_coords(coords, coords_filename):
    """
    Writes coordinates to xyz format.

    :param coords: Numpy array-like (ndarray or list) of coordinates to be printed.
    :param coords_filename: File name to be written to.
    """
    coords = np.reshape(coords, (-1, 3))
    with open(coords_filename, "w") as coords_file:
        coords_file.write("{: >10d}\n".format(coords.shape[0]))
        for coord in coords:
            coords_file.write(("{: >20.10f}" * 3 + "\n").format(*coord))


if __name__ == "__main__":
    lowest_structs = read_lowest("/scratch/khs26/ser_lys_fe_bh/ff03_igb2/temp_0.0/rep5/output/lowest")
    write_coords(np.array(lowest_structs[0][-1]), "test_output")
