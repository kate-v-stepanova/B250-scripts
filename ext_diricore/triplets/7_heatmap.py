import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

import os
import sys
import glob

sns.set(rc={'figure.figsize':(10, 16)})

infile = sys.argv[1]
outfile = sys.argv[2]


data = pd.read_csv(infile, sep="\t")
genes = data['gene'].tolist()
data.index = genes
data = data.drop('gene', axis=1)
import pdb; pdb.set_trace()
data = data[:80]
data.columns = [col.replace('norm_', '') for col in data.columns]

g = sns.clustermap(data, cmap="vlag", robust=True, col_cluster=False, xticklabels=True, yticklabels=True, figsize=(10,16), cbar_kws={"ticks":[0,3]})
g.ax_heatmap.set_title('Heatmap for the subset of genes: {}'.format(os.path.basename(infile).replace('.tsv', '')))
plot = g.fig
print('Writing file: {}'.format(outfile))
plot.savefig(outfile)
