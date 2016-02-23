import rotamer.library.read_rotamer_csv as read_csv

# Restrict range of symmetric dihedrals to -90 to 90
# Get rid of arginine chi5

symmetric_dihedrals = {"ASP": "chi2",
                       "GLU": "chi3",
                       "LEU": "chi2",
                       "PHE": "chi2",
                       "TYR": "chi2",
                       "VAL": "chi1"}

def get_csv_list(directory):
    import os.path
    import glob
    all_csvs = glob.glob(os.path.join(directory, "*.csv"))
    matching_csvs = [csv_file for csv_file in all_csvs if any([k in csv_file for k in symmetric_dihedrals])]
    return matching_csvs

def columns_to_fix(df):
    """
    Returns the list of columns to be fixed.

    :param df:
    :return: list of column names
    """
    return [col for col in df.columns.values if any([k in col and v in col for k, v in symmetric_dihedrals.items()])]

def fix_values(df, columns):
    """
    Changes the values of the dihedrals in the columns that need to be fixed.

    :param df: Dataframe
    :param columns: Columns to be fixed
    """
    df[df.loc[:, columns] > 90] -= 180
    df[df.loc[:, columns] < -90] += 180
    arg_chi5s = [col for col in df.columns.values if "ARG" in col and "chi5" in col]
    return df.drop(arg_chi5s, axis=1)


if __name__ == "__main__":
    import re
    matching_csvs = get_csv_list("/home/khs26/rotamer_library/dihedrals")
    for i, csv in enumerate(matching_csvs):
        print csv, "{0}/{1}".format(i, len(matching_csvs))
        df = read_csv.csv_to_dataframe(csv)
        new_df = fix_values(df, columns_to_fix(df))
        formatter = {col: lambda x: "{0: 9.4f},".format(x) for col in new_df.columns.values}
        formatter["energy"] = lambda x: "{0:<.8f},".format(x)
        as_string = new_df.to_string(header=False, index=False, formatters=formatter, justify="left")
        fix_string = re.sub(",$", "", re.sub(",\n", "\n", re.sub("\n ", "\n", re.sub("^ ", "", re.sub(", ", ",", as_string)))))
        with open(csv, "r") as csv_in:
            with open("_".join([csv, "fixed"]), "w") as csv_out:
                for i in range(0, 5):
                    line_in = csv_in.readline()
                    if "ARG" in line_in:
                        line_out = re.sub(",chi5", "", line_in)
                    else:
                        line_out = line_in
                    csv_out.write(line_out)
                csv_out.write(fix_string)
