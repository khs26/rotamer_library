import os
import os.path
import glob
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster
import numpy as np
from rotamer.topology.residue_sidechains import amino_acids

def find_csvs(residue_string):
    """
    Returns a list of CSVs with the specified central residue.

    :param residue_string: 3-letter abbreviation of the residue
    :return: list of csv paths
    """
    root = os.path.abspath(os.sep)
    csv_directory = os.path.join(root, "home", "khs26", "rotamer_library", "dihedrals")
    csv_list = glob.glob(os.path.join(csv_directory, "*_" + residue_string + "_*.csv"))
    return csv_list

def csv_to_dataframe(csv_filename):
    """
    Generates a dataframe with appropriate column names from a rotamer library csv file.

    :param csv_filename: Name of CSV file to open
    :return: Dataframe with rotamer data
    """
    rotamer_data = pd.read_csv(csv_filename, skiprows=5, header=None)
    with open(csv_filename, "r") as csv_file:
        residues = csv_file.readline().split(",")
        global_min = float(csv_file.readline().split(",")[1])
        column_names = ["energy"]
        for i in range(0, len(residues)):
            split_line = csv_file.readline().split(",")
            column_names += ["_".join((residues[i].strip(), x.strip())) for x in split_line[1:]]
        rotamer_data.columns = column_names
    return rotamer_data

def get_columns_matching(dataframe, pattern, regex=False):
    """
    Returns the columns containing pattern (or satisfying the regex pattern, if regex=True).

    :param dataframe: Dataframe whose columns are being tested
    :param pattern: String to test columns for
    :param regex: If True, then check pattern as regex
    :return: Dataframe slice of matching columns
    """
    if not regex:
        matching = [col for col in dataframe.columns.values.tolist() if pattern in col]
    else:
        import re
        matching = [col for col in dataframe.columns.values.tolist() if re.match(pattern, col)]
    return dataframe.loc[:, matching]

def angle_distance2(angle1, angle2):
    """
    Returns the square of the Euclidean distance between two pairs of rotamer states, accounting for the wrap around of
    angle values.

    :param angle1, angle2: numpy arrays containing the angles (in degrees)
    :return: Euclidean distance squared, accounting for wrap
    """
    raw_deltas = np.absolute(angle1 - angle2)
    wrapped_deltas = np.fmin(raw_deltas, (360.0 - raw_deltas))
    return np.sum(np.square(wrapped_deltas))

if __name__ == "__main__":
    import time
    amino_acid = "ARG"
    arg_csvs = find_csvs(amino_acid)
    all_dfs = [csv_to_dataframe(csv_name) for csv_name in arg_csvs]
    only_centre = [get_columns_matching(df, r"(energy)|(2 " + amino_acid + ")", True) for df in all_dfs]
    joined = pd.concat(only_centre)
    # NDArray of rotamer angles
    as_ndarray = joined.iloc[:, 1:].values
    # Create and use a DBScan instance
    print "Started clustering:"
    now = time.time()
    dbscan = sklearn.cluster.DBSCAN(eps=200, metric=angle_distance2).fit_predict(as_ndarray)
    end = time.time()
    print end - now
    print dbscan
    n_clusters = len(set(dbscan)) - (1 if -1 in dbscan else 0)
    print n_clusters


# # print rotamer_data.head()
# correlations = rotamer_data.corr()
# print correlations[correlations < 1.0].describe()
# plot = scatter_matrix(rotamer_data[["1 SER_phi", "1 SER_psi", "1 SER_chi1"]], alpha=0.2)
# print plot
# for k, v in plot[0][1].__dict__.items():
#     print k, v
# plot[0][0].set_xlim(-180, 180)
# plot[0][1].set_xlim(-180, 180)
# plot[0][1].set_ylim(-180, 180)
# plot[1][0].set_xlim(-180, 180)
# plot[1][0].set_ylim(-180, 180)
# plot[1][1].set_xlim(-180, 180)
# plt.show()