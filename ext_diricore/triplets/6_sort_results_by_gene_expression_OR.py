import sys
import os

import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep="\t")
df = df.drop('aa', axis=1)

import pdb; pdb.set_trace()
df.index = df['gene'].tolist()
df = df.drop('gene', axis=1)
df = df.drop('norm_CAF_Ctrl_si_Pro', axis=1)
df = df.loc[(df!=0).any(axis=1)]
df['diff1'] = df['norm_CAF_PYCR1_si'] - df['norm_CAF_Ctrl_si']
df['diff2'] = df['norm_CAF_PYCR1_si_Pro'] - df['norm_CAF_PYCR1_si']
df = df.sort_values(by='diff1', ascending=False)

df = df.loc[(df['diff1'] > 0) | (df['diff2'] < 0)]

df['gene'] = df.index
df = df[['gene', 'norm_CAF_Ctrl_si', 'norm_CAF_PYCR1_si', 'norm_CAF_PYCR1_si_Pro']]
print('Writing {}'.format(outfile))
df.to_csv(outfile, sep="\t", index=False)

