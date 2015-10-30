import numpy as np

# What the selected rotamer move should take into account when selecting
# a new state.
depends = {"backbone": False,
           "neighbours": False,
           "neighbours_sidechain": False}

# Proposal acceptance temperature
temperature = None

# Energy distribution function
distribution = None

def read_settings(filename):
    """ Reads settings for rotamer library moves from the specified file.

    :param filename: Filename to read settings from
    :return: True if settings read correctly
    """
    with open(filename, "r") as settings_file:
        pass
    return check_settings()

def check_settings():
    """ Validates the current settings.

    :return: True if all necessary settings are set.
    """
    return temperature is not None and distribution is not None

if __name__ == "__main__":
    pass