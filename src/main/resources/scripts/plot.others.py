# -*- coding: utf-8 -*-

import matplotlib as mpl
from _socket import NI_DGRAM
# mpl.use('Agg')
import matplotlib.pyplot as p
import numpy as np
import pylab
import scipy.stats as stats
import networkx as nwx
import glob
import builtins
from matplotlib.pyplot import margins
import os.path
import json, sys

import pandas as pd

from collections import OrderedDict

RED = 'darkred'
G1 = 'slategrey'
G2 = 'darkgrey'
         
df = pd.read_csv('./src/main/resources/data/motif-experiments.csv')

# print(df)
print(df.keys())

types = ['ours', 'full motif', 'approximate motif', 'full census', 'approximate census']
markers = ['s', 'o', 'd', 'x', '+']

fig = p.figure(figsize=(8, 6))
norm = mpl.colors.Normalize(vmin=df['year'].min(), vmax=df['year'].max())

for k, m in zip(types, markers):
    print(k)
    dfs = df.loc[df['type'] == k]
    p.scatter(dfs['nodes'].values.astype('float'), dfs['links'].values.astype('float'),
                marker=m, c=dfs['year'].values.astype('float'), cmap='copper_r', label=k, norm=norm)
    
    print(len(dfs.index))
#      
#     if k != 'ours':
#         for i in range(len(dfs.index)):
#             xy = (dfs['nodes'].values[i], dfs['links'].values[i])
#             p.annotate(int(dfs['source'].values[i]), xy, xycoords='data', xytext=(5, -1),
#                 textcoords='offset points')


ax = p.gca()
ax.set_xlim([10, 1e9])
ax.set_ylim([10, 1e10])

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_aspect('equal')

ax.set_xlabel('number of nodes')
ax.set_ylabel('number of links')
    
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
 
p.colorbar()
leg = p.legend(frameon=False)

for h in leg.legendHandles:
    h.set_color('k')

p.savefig('others.pdf')

