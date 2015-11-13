import os
import os.path
import glob
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster
import numpy as np
from pandas.tools.plotting import scatter_matrix, andrews_curves
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

def cluster_angles(dataframe, eps, metric=angle_distance2, max_unclustered=0.00):
    """
    Assigns each dihedral angle to a cluster using the DBScan algorithm. If the clustering algorithm returns a higher
    proportion of unclustered points than max_unclustered, repeat the clustering with a doubled eps until a suitable
    number have been clustered.

    :param dataframe: dataframe to cluster (this can also be a slice of a larger dataframe)
    :param eps: epsilon value to cutoff nearest neighbour analysis
    :param metric: callable, which return a distance between two angles
    :param max_unclustered: proportion of the data which can be unclustered
    :return: dataframe in which angles have been replaced by cluster ids
    """
    as_ndarray = dataframe.values
    # Create a dataframe containing cluster labels, instead of angles
    label_df = pd.DataFrame()
    # Run the clustering until we satisfy the max_unclustered condition.
    pd.set_option('display.width', 1000)
    for i, col in enumerate(as_ndarray.transpose()[:]):
        unclustered = col.shape[0]
        this_eps = eps
        while unclustered/float(col.shape[0]) > max_unclustered:
            dbscan = sklearn.cluster.DBSCAN(eps=this_eps, metric=metric).fit(col[np.newaxis, :].transpose())
            labels = dbscan.labels_
            cluster_dict = {}
            for label in set(labels):
                 if label == -1:
                     continue
                 cluster_dict[str(label)] = pd.Series(col[labels == label])
            clusters = pd.DataFrame(cluster_dict)
            # pd.DataFrame(col).hist(bins=36)
            clusters.plot(kind='kde')
            plt.show()
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            unclustered = np.sum(labels == -1)
            this_eps = 1.1 * this_eps
            # clustered = pd.DataFrame()
            # with open('clusters', "w") as clusters:
            #     for j in range(-1, n_clusters):
            #         clustered = pd.concat([clustered, pd.DataFrame([col[l] for l, k in enumerate(labels) if k == j])], axis=1)
            #         clusters.write("Cluster {: d}\n".format(j))
            #         clusters.write("-----------------\n")
            #         for element in clustered.iloc[:, j][clustered.iloc[:, j].notnull()].tolist():
            #             clusters.write("{: 8.3f}\n".format(element))
            #         clusters.write("-----------------\n")
            #         # mean = np.sum(clustered.iloc[:, j]) / clustered.iloc[:, j].size
            #         # stdev = np.sqrt(np.sum(np.square(clustered.iloc[:, j] - mean)) / clustered.iloc[:, j].size)
            #         # print "{: d} mean: {: 8.3f} stdev: {: 8.3f}".format(j, mean, stdev)
            #         # outliers = (clustered.iloc[:, j] - mean) > (3 * stdev)
            #         # print "-----------"
            #         # print "Outliers:", clustered.iloc[:, j][outliers]
            #         # print "-----------"
            #         clustered.columns[j] == str(j)
            #     clusters.write("------------------------------------------------\n")
            #     print clustered.transpose()
            #     # plot(kind='kde')
            #     plt.show()
            #     exit()
            # print "".join(["{: 7.2f}{: 7.2f}\n".format(mu, sigma) for (mu, sigma) in sorted(zip(clustered.mean(0), clustered.std(0)))])
            # clustered.hist(bins=36)
            # plt.show()
        # Add this column to the label dataframe
        label_df = pd.concat([label_df, pd.DataFrame(labels)], axis=1)
        print "==================="
        print dataframe.columns[i]
        print "n_clusters:", n_clusters, "unclustered:", unclustered
        print "final eps:", this_eps / 1.1
        for j in range(0, n_clusters):
            print j, ":", np.sum(labels == j)
    # Copy the column names over
    label_df.columns = dataframe.columns
    return label_df

if __name__ == "__main__":
    from scipy import stats
    amino_acid = "ARG"
    arg_csvs = find_csvs(amino_acid)
    all_dfs = [csv_to_dataframe(csv_name) for csv_name in arg_csvs]
    only_centre = [get_columns_matching(df, r"(energy)|(2 " + amino_acid + ")", True) for df in all_dfs]
    joined = pd.concat(only_centre)
    labelled = cluster_angles(joined.iloc[:10000, 3:], 0.1, metric=wrapped_distance)
    print labelled
    # NDArray of rotamer angles
    as_ndarray = labelled.values.tolist()
    as_ndarray = [tuple(ls) for ls in as_ndarray]
    print len(set(as_ndarray))
    joined.hist(bins=36)
    plt.show()
    # print sorted(as_ndarray)
    # print "Started clustering:"
    # dbscan = sklearn.cluster.DBSCAN(eps=1, metric=sum_delta_function).fit(as_ndarray)
    # print dbscan.labels_
    # # scatter_matrix(joined.iloc[:1000, 3:], diagonal='kde')
    # # plt.show()
    # # for i in np.linspace(1, 200, 10):
    # #     cluster(i)
    # print "Epsilon:", 200
    # now = time.time()
    # # dbscan = sklearn.cluster.DBSCAN(eps=200, metric=angle_distance2).fit_predict(as_ndarray)
    # labels = dbscan.labels_
    # end = time.time()
    # print end - now, "s"
    # # print dbscan
    # n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    # print "Clusters:", n_clusters, "Unclustered:", np.sum(labels == -1)
    # # print "====================="
    # core_samples_mask = np.zeros_like(dbscan.labels_, dtype=bool)
    # core_samples_mask[dbscan.core_sample_indices_] = True
    # unique_labels = set(labels)
    # np.set_printoptions(precision=2, suppress=True)
    # cluster_reps = []
    # for i in unique_labels:
    #     class_member_mask = (labels == i)
    #     for j, val in enumerate(class_member_mask):
    #         if val:
    #             print i, as_ndarray[j]
    #             cluster_reps.append(as_ndarray[j])
    #             break
    # for j, i in enumerate(sorted(cluster_reps, key=lambda x: (x[0], x[1], x[2], x[3], x[4]))):
    #     print j, i
    # colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    # for k, col in zip(unique_labels, colors):
    #     # if not (k == 3 or k == 15 or k == 32 or k == 54):
    #     #     continue
    #     if k == -1:
    #         # Black used for noise.
    #         col = 'k'
    #     class_member_mask = (labels == k)
    #
    #     xy = as_ndarray[class_member_mask & core_samples_mask]
    #     plt.plot(xy[:, 0], xy[:, 2], 'o', markerfacecolor=col, markeredgecolor='k', markersize=14)
    #
    #     xy = as_ndarray[class_member_mask & ~core_samples_mask]
    #     plt.plot(xy[:, 0], xy[:, 2], 'o', markerfacecolor=col, markeredgecolor='k', markersize=6)
    # plt.xlim([-180, 180])
    # plt.ylim([-180, 180])
    # plt.title('Estimated number of clusters: %d' % n_clusters)
    # plt.show()
    # # Create and use a DBScan instance
    # # print "Started clustering:"
    # # now = time.time()
    # # dbscan = sklearn.cluster.DBSCAN(eps=300, metric=angle_distance2).fit_predict(as_ndarray)
    # # end = time.time()
    # # print end - now, "s"
    # # print dbscan
    # # n_clusters = len(set(dbscan)) - (1 if -1 in dbscan else 0)
    # # print n_clusters