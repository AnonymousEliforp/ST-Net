import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# initialize dataframe
n = 200
ngroup = 3
df = pd.DataFrame({'data': np.random.rand(n), 'group': map(np.floor, np.random.rand(n) * ngroup)})

group = 'group'
column = 'data'
grouped = df.groupby(group)

for key, item in grouped:
    print(grouped.get_group(key), "\n\n")