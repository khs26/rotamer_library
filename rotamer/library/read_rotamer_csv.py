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

def wrapped_distance(point1, point2, wrap=360.0):
    """
    Returns the distance between two points on a circular surface, which wraps at `wrap`.

    :param point1, point2: coordinates of the points
    :param wrap: wrapping value
    :return: distance
    """
    delta = abs(point1 - point2)
    return np.sum(np.fmin(delta, wrap - delta))

def angle_distance2(angles1, angles2):
    """
    Returns the square of the Euclidean distance between two pairs of rotamer states, accounting for the wrap around of
    angle values.

    :param angles1, angles2: numpy arrays containing the angles (in degrees)
    :return: Euclidean distance squared, accounting for wrap
    """
    raw_deltas = np.absolute(angles1 - angles2)
    wrapped_deltas = np.fmin(raw_deltas, (360.0 - raw_deltas))
    return np.sum(np.square(wrapped_deltas))

def sum_delta_function(array1, array2):
    """
    Returns the number of identical array elements

    :param array1, array2: numpy arrays containing values for comparison
    :return: Number of identical array elements
    """
    return np.sum(array1 != array2)

def find_gap_in_data(data, width, min, max):
    """
    Finds gaps in data (where there are no values) of specified width, between min and max.

    :param data: array-like containing data
    :param width: width of empty window (and step size)
    :param min: bottom of search space
    :param max: top of search space
    :return: tuple containing lowest value (or None if none found)
    """
    ret_tuple = None
    nsteps = (max - min) / float(width)
    for bottom in np.linspace(min, max-width, num=nsteps):
        top = bottom + width
        # If our data is all either lower than bottom or greater than top, return (bottom, top)
        if np.all(np.logical_or(data < bottom, data > top)):
            ret_tuple = (bottom, top)
            break
    return ret_tuple


def cluster_wrapper(data, algorithm=None, **kwargs):
    """
    Wrapper for the clustering algorithms available in scikit-learn.

    :param data: Data to be clustered (in a format suitable to pass straight to the clustering algorithms).
    :param algorithm: Callable algorithm to use.
    :param kwargs: Keywords for the algorithm used.
    :return: - pandas DataFrame with data in columns corresponding to clusters.
             - pandas Series with labels in the place of angles.
             - sklearn.cluster object, which can be used to predict cluster identities of individual angles.
    """
    # Shift data so that there isn't a break in clusters around -180 and 180
    gap = find_gap_in_data(data, width=0.5, min=-180, max=180)
    data[data < gap[0]] += 360.0
    # Predict clusters
    algo_object = algorithm(**kwargs)
    labels = algo_object.fit_predict(data)
    # Shift data back to (-180, 180)
    data[data > 180.0] -= 360.0
    # Convert labels to a DataFrame
    data_1d = data.transpose()[0]
    label_series = {}
    for label in set(labels):
        label_series[label] = pd.Series(data_1d[labels == label])
    return pd.DataFrame(label_series), pd.Series(labels), algo_object


def cluster_angles_kmeans(dataframe, n_clusters):
    """
    Assigns each dihedral angle to a cluster using the mini-batch K-means algorithm. The number of clusters for each
    angle must be determined beforehand (e.g. from inspection).

    :param dataframe: dataframe to cluster (this can also be a slice of a larger dataframe)
    :param n_clusters: number of clusters to generate
    :return: dataframe in which angles have been replaced by cluster ids
    """
    as_ndarray = dataframe.values
    labelled_states = {}
    clusterers = {}
    for i, col in enumerate(as_ndarray.transpose()[:]):
        clustered, labels, cluster_object = cluster_wrapper(data=col[np.newaxis, :].transpose(),
                                                            algorithm=sklearn.cluster.MiniBatchKMeans,
                                                            n_clusters=n_clusters[i])
        labelled_states[dataframe.columns[i]] = labels
        clusterers[dataframe.columns[i]] = cluster_object
    #     clustered.plot(kind='hist', bins=360, edgecolor='none', normed=1)
    # plt.show()
    labelled = pd.DataFrame(labelled_states)
    return labelled, clusterers

if __name__ == "__main__":
    """
    TODO:
        - Dump list of cluster centroids for each dihedral.
        - Reweight probabilities with Boltzmann distribution (at 298K = 0.59kcal/mol).
        - Tabulate the probabilities.
    """
    import cPickle
    import copy
    import itertools as it
    for amino_acid in amino_acids:
        if amino_acid in ("ALA", "GLY"):
            continue
        if not amino_acid in ("LEU"):
            continue
        ok = False
        while not ok:
            arg_csvs = find_csvs(amino_acid)
            all_dfs = [csv_to_dataframe(csv_name) for csv_name in arg_csvs]
            only_centre = [get_columns_matching(df, r"(energy)|(2 " + amino_acid + ")", True) for df in all_dfs]
            joined = pd.concat(only_centre)
            joined.iloc[:, 3:].hist(bins=360, weights=np.exp(-joined.iloc[:, 0] / 0.59).reshape(-1, 1), normed=1, edgecolor='none')
            # joined.iloc[:, 3:].hist(bins=360, normed=1)
            # plt.show()
            clusters = map(int, raw_input("Cluster proposals: ").split())
            labelled, clusterers = cluster_angles_kmeans(joined.iloc[:, 3:], clusters)
            # ok = "y" in raw_input("OK? ")
            ok = True
        normalisation = np.sum(np.exp(-joined.iloc[:, 0] / 0.59))
        weights = np.exp(-joined.iloc[:, 0] / 0.59).values
        labelled["weight"] = weights
        states = labelled.iloc[:, :-1].drop_duplicates().values
        as_ndarray = labelled.iloc[:, :-1].values
        prob = {}
        for state in states:
            angles = []
            for i in range(len(state)):
                angles.append(clusterers[labelled.columns.values[i]].cluster_centers_[state[i]][0])
            prob[tuple(state)] = np.sum(labelled[np.all(state == as_ndarray, axis=1)].iloc[:, -1]) / normalisation
            print " ".join(["{0: >8.3f}".format(angle) for angle in angles]) + "{0: >7.2f}%".format(100 * prob[tuple(state)])
        print sum(prob.values())
            # print labelled[np.all(state == as_ndarray, axis=1)]
        # print labelled

    # # NDArray of rotamer angles
    # as_ndarray = labelled.values
    # labels = labelled.drop_duplicates().values
    # # Join it up with the phi/psi angles and the energies
    # all_with_labels = pd.concat([pd.DataFrame(joined.iloc[:, :3].values), labelled], axis=1)
    # all_with_labels.columns.values[:3] = joined.columns.values[:3]
    # # print all_with_labels[np.logical_and((all_with_labels["2 HIP_chi1"] == 0), (all_with_labels["2 HIP_chi2"] == 0))]
    # print labels
    # for label in labels:
    #     print "state:", label, "count:", np.sum(np.all(label == as_ndarray, axis=1))