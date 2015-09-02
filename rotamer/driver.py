import os
import site
import sys
# Add
site.addsitedir(os.path.normpath(sys.argv[0] + "/../.."))
import rotamer.io.amber as amber_io
import rotamer.io.gmin as gmin_io

if __name__ == "__main__":
    if sys.argv[1] == "init":
        pass
    elif sys.argv[1] == "move":
        coords = amber_io.read_amber_restart(".coords_before_rotamer.rst")
        coords = coords * 100
        gmin_io.write_coords(coords, ".coords_after_rotamer.rst")