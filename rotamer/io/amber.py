import numpy as np
import os.path


def read_amber_restart(filename):
    """
    Reads coordinates from an AMBER restart file.

    :param filename: File name to be read from.
    :return: Numpy array of coordinates with shape (natoms, 3).
    """
    field_length = 12
    coords = []
    with open(filename, "r") as coords_file:
        # Throw away the first line, since it just contains the name of the molecule.
        coords_file.readline()
        # The next line contains the number of atoms.
        number_of_atoms = int(coords_file.readline())
        # Later lines contain coordinates in 12-character wide fields.
        for line in coords_file:
            line = line.rstrip()
            coords += map(float, [line[i:i + field_length] for i in range(0, len(line), field_length)])
    # If the number of coordinates is not equal to 3 * number of atoms, raise a RuntimeError.
    if len(coords) != number_of_atoms * 3:
        raise RuntimeError("Number of coordinates in coords file and number of atoms are inconsistent.")
    return np.reshape(coords, (-1, 3))


def write_amber_restart(filename, coords):
    """
    Writes coordinates to a file in AMBER restart format.

    :param filename: File name to write to.
    :param coords: Numpy array-like (list or ndarray) of coordinates.
    """
    coords = np.reshape(coords, (-1, 3))
    with open(filename, "w") as restart_file:
        restart_file.write("{0}\n".format(os.path.basename(filename)))
        restart_file.write("{0: >5d}\n".format(coords.shape[0]))
        buff = []
        # Coordinates are printed in the format 6F12.7. If there are an odd number
        # of atoms, this means that we'll have one line that is 3F12.7 at the end.
        for i, coord in enumerate(coords):
            buff += coord.tolist()
            if i % 2 == 1:
                restart_file.write(("{: >12.7f}" * 6 + "\n").format(*buff))
                buff = []
        if buff:
            restart_file.write(("{: >12.7f}" * 3 + "\n").format(*buff))


if __name__ == "__main__":
    write_amber_restart("test.rst", [0, 12, 15, 2, 0.5, 3, 4, 5, 6])
