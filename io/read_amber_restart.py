__author__ = 'khs26'


def read_amber_coords(filename):
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
        raise ex.RuntimeError("Number of coordinates in coords file and number of atoms are inconsistent.")
    return coords