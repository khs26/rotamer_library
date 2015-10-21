import numpy as np
import pandas as pd
import sys

rotamer_data = pd.read_csv(sys.argv[1], skiprows=5)
print rotamer_data