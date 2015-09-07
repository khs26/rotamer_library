import numpy as np


def read_one_structure(lowest_file):
    natoms = int(lowest_file.readline().strip())
    data_line = lowest_file.readline().split()
    index = int(data_line[3][:-1])
    energy = float(data_line[4])
    first_found = int(data_line[-1])
    coords = np.zeros([natoms, 3])
    for i in range(0, natoms):
        coords[i] = map(float, lowest_file.readline().split()[1:])
    return index, energy, first_found, coords


def read_lowest(lowest_filename):
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
    coords = np.reshape(coords, (-1, 3))
    with open(coords_filename, "w") as coords_file:
        coords_file.write("{: >10d}\n".format(coords.shape[0]))
        for coord in coords:
            coords_file.write(("{: >20.10f}" * 3 + "\n").format(*coord))


if __name__ == "__main__":
    lowest_structs = read_lowest("/scratch/khs26/ser_lys_fe_bh/ff03_igb2/temp_0.0/rep5/output/lowest")
    write_coords(np.array(lowest_structs[0][-1]), "test_output")