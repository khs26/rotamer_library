from rotamer.topology.residue_sidechains import amino_acids
import rotamer.dihedral.find_dihedrals as find_dihedrals
from rotamer.io.gmin import read_lowest
import rotamer.topology.read_amber_prmtop as ra
import numpy as np


class BaseRotamerState(object):
    """
    Base class for rotamer states for the GMIN rotamer library.
    """
    pass


class FullRotamerState(BaseRotamerState):
    """
    Defines a rotamer state for use with the GMIN rotamer library.

    A rotamer state corresponds to a single minimum from the rotamer library simulations. Each state has the values of
    all the dihedrals in the corresponding tripeptide and the energy (relative to the global minimum for that
    tripeptide).
    """
    pass


class RotamerStateFactory(object):
    """
    Creates new rotamer states, either from raw GMIN data or from a HDF5 database.
    """

    def from_lowest(self, prmtop_filename, lowest_filename, state_type):
        """
        Creates a new RotamerState of appropriate type from a GMIN lowest file.

        :param prmtop_filename: File name of the AMBER topology file
        :param lowest_filename: File name of the GMIN lowest file
        :param state_type: Type of rotamer state (i.e. what sort of dependence on neighbours, backbone dihedrals etc.)
        """
        molecule = ra.parse_topology_file(prmtop_filename)
        find_dihedrals.map_dihedrals(*(molecule.identify_residues()))
        self.configs = read_lowest(lowest_filename)
        find_dihedrals.phi_psi_dihedrals(molecule)
        find_dihedrals.sidechain_dihedrals(molecule)
        lowest0 = self.lowest_config_coords(0)
        lowest25 = self.lowest_config_coords(25)
        for res in sorted([res for res in molecule.residues.nodes() if res.name != "ACE" and res.name != "NME"],
                          key=lambda x: x.index):
            print "---------------"
            print res.name, res.index
            print "---------------"
            print "phi", res.dihedrals["phi"].measure_dihedral(lowest0) * 180.0 / np.pi
            print "psi", res.dihedrals["psi"].measure_dihedral(lowest0) * 180.0 / np.pi
            for dihedral, dihe_object in sorted([(k, v) for k, v in res.dihedrals.items() if 'chi' in k]):
                print dihedral, dihe_object.measure_dihedral(lowest0) * 180.0 / np.pi
            print ""
        for res in sorted([res for res in molecule.residues.nodes() if res.name != "ACE" and res.name != "NME"],
                          key=lambda x: x.index):
            print "---------------"
            print res.name, res.index
            print "---------------"
            print "phi", res.dihedrals["phi"].measure_dihedral(lowest25) * 180.0 / np.pi
            print "psi", res.dihedrals["psi"].measure_dihedral(lowest25) * 180.0 / np.pi
            for dihedral, dihe_object in sorted([(k, v) for k, v in res.dihedrals.items() if 'chi' in k]):
                print dihedral, dihe_object.measure_dihedral(lowest25) * 180.0 / np.pi
            print ""

    def lowest_config_coords(self, index):
        return self.configs[index]['coords']

def norm_lowest_configs(lowest_configs):
    global_min_energy = min([x["energy"] for x in lowest_configs])
    for config in lowest_configs:
        config["energy"] -= global_min_energy
    return global_min_energy

def lowest_to_dihedral_csv(prmtop_filename, lowest_filename, csv_filename):
    """
    Reads a topology and lowest filename and writes dihedral information to a csv file.

    The file has the following structure:

    <res_1>,<res_2>,<res_3>,...,<res_n>
    global minimum,<energy of global minimum>
    <res_1>,phi,psi,[chi_1],...,[chi_n]
    <res_2>,phi,psi,[chi_1],...,[chi_n]
    .
    .
    .
    <res_n>,phi,psi,[chi_1],...,[chi_n]
    <norm_energy_0>,<res_1_phi>,<res_1_psi>,[<res_1_chi_1>],...,<res_2_phi>,...,<res_n_phi>,...,[<res_n_chi_n>]
    <norm_energy_1>,<res_1_phi>,<res_1_psi>,[<res_1_chi_1>],...,<res_2_phi>,...,<res_n_phi>,...,[<res_n_chi_n>]
    .
    .
    .
    <norm_energy_N>,<res_1_phi>,<res_1_psi>,[<res_1_chi_1>],...,<res_2_phi>,...,<res_n_phi>,...,[<res_n_chi_n>]

    N.B. energies are given relative to the global minimum and angles are given in degrees in the order specified at the
    top of the file

    :param prmtop_filename: Amber topology filename
    :param lowest_filename: GMIN lowest filename
    :return:
    """
    import csv

    molecule = ra.parse_topology_file(prmtop_filename)
    find_dihedrals.map_dihedrals(*(molecule.identify_residues()))
    find_dihedrals.phi_psi_dihedrals(molecule)
    find_dihedrals.sidechain_dihedrals(molecule)
    lowest_configs = read_lowest(lowest_filename)
    res_no_termini = [res for res in sorted(molecule.residues.nodes(), key=lambda x: x.index)
                      if res.name != "ACE" and res.name != "NME"]
    with open(csv_filename, "w") as csv_file:
        dihe_writer = csv.writer(csv_file)
        dihe_writer.writerow(res_no_termini)
        global_min_energy = norm_lowest_configs(lowest_configs)
        dihe_writer.writerow(["global minimum", global_min_energy])
        ordered_dihe = []
        for res in res_no_termini:
            backbone = sorted([dihe for dihe in res.dihedrals if 'p' in dihe])
            sidechain = sorted([dihe for dihe in res.dihedrals if 'p' not in dihe])
            dihe_writer.writerow([res.name] + backbone + sidechain)
            ordered_dihe += [res.dihedrals[x] for x in backbone]
            ordered_dihe += [res.dihedrals[x] for x in sidechain]
        for config in lowest_configs:
            energy = config["energy"]
            angles = [x.measure_dihedral(config["coords"]) * 180.0 / np.pi for x in ordered_dihe]
            dihe_writer.writerow(["{:.8f}".format(energy)] + ["{: 9.4f}".format(angle) for angle in angles])


if __name__ == "__main__":
    from rotamer.topology.residue_sidechains import amino_acids
    import itertools
    import os.path
    import time
    now = time.time()
    for comb in itertools.product(amino_acids, repeat=3):
        print comb, "{:.1f}s".format(time.time() - now)
        directory = os.path.join(os.path.curdir, *comb)
        prmtop = os.path.join(directory, "coords.prmtop")
        lowest = os.path.join(directory, "lowest")
        dihedral_csv = os.path.join(os.path.curdir, "dihedrals", "".join(("_".join(comb), ".csv")))
        lowest_to_dihedral_csv(prmtop, lowest, dihedral_csv)