# -*- coding: utf-8 -*-

import matplotlib as mpl
from _socket import NI_DGRAM
mpl.use('Agg')
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

RED = 'darkred'
G1 = 'slategrey'
G2 = 'darkgrey'

def plot(df, color, marker):
    
    mean = df.groupby([0]).mean()
    sem = df.groupby([0]).sem()
    
    ax = p.subplot(321)
    
    handles = ax.scatter(mean[2].values, mean[6].values-mean[5].values, color=color, marker=marker)
    ax.errorbar(mean[2].values, mean[6].values-mean[5].values, yerr=1.96*sem[6], fmt='none', color=RED, elinewidth=3)

    z, _, _, _, _ = np.polyfit(np.log10(mean[2].values), np.log10(mean[6].values-mean[5].values), 1, full=True)
    xs = np.logspace(np.log10(mean[2].min()), np.log10(mean[2].max()), num=500)
    hn, = ax.plot(xs, 10**np.poly1d(z)(np.log10(xs)), linestyle='dotted', zorder=-10, color=color, label='%.03f' % z[0])
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)   
    ax.spines["bottom"].set_visible(False)
    ax.get_xaxis().set_tick_params(labelbottom='off')
    
    ax.set_ylabel('motif tests (s)')

    ax1 = ax
    ax = p.subplot(322, sharey=ax1)

    ax.scatter(mean[3].values, mean[6].values-mean[5].values, color=color, marker=marker)
    ax.errorbar(mean[3].values, mean[6].values-mean[5].values, yerr=1.96*sem[6], fmt='none', color=RED, elinewidth=3)
    
    z, _, _, _, _ = np.polyfit(np.log10(mean[3].values), np.log10(mean[6].values-mean[5].values), 1, full=True)
    xs = np.logspace(np.log10(mean[3].min()), np.log10(mean[3].max()), num=500)
    hl, = ax.plot(xs, 10**np.poly1d(z)(np.log10(xs)), linestyle='dotted', zorder=-10, color=color, label='%.03f' % z[0])
        
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([0.1, 7000])
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.get_xaxis().set_tick_params(labelbottom='off')
    ax.get_yaxis().set_tick_params(labelleft='off')

    
    ax2 = ax
    
    ax = p.subplot(323, sharex=ax1)
    
    ax.scatter(mean[2].values, mean[5].values, color=color, marker=marker)
    ax.errorbar(mean[2].values, mean[5].values, yerr=1.96*sem[5], fmt='none', color=RED, elinewidth=3)

    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)   
    ax.spines["bottom"].set_visible(False)
    ax.get_xaxis().set_tick_params(labelbottom='off')
    
    ax.set_ylabel('sampling (s)')

    ax3 = ax
    ax = p.subplot(324, sharey=ax3, sharex=ax2)

    ax.scatter(mean[3].values, mean[5].values, color=color, marker=marker)
    ax.errorbar(mean[3].values, mean[5].values, yerr=1.96*sem[5], fmt='none', color=RED, elinewidth=3)

    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.get_xaxis().set_tick_params(labelbottom='off')
    ax.get_yaxis().set_tick_params(labelleft='off')

    ax4 = ax

    ax = p.subplot(325, sharex=ax1)
    min = df.groupby([0]).min()
    max = df.groupby([0]).max()
    ax.scatter(mean[2].values, mean[4].values, color=color, marker=marker)
    ax.vlines(x=mean[2].values, ymin=min[4], ymax=max[4], color=RED, linewidth=3)
    ax.set_ylim([-5, 105])
    ax.set_xscale('log')
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    
    ax.set_xlabel('nr. of nodes')
    ax.set_ylabel('nr. of motifs')
    
    ax5 = ax
    
    ax = p.subplot(326, sharey=ax5, sharex=ax2)
    min = df.groupby([0]).min()
    max = df.groupby([0]).max()
    ax.scatter(mean[3].values, mean[4].values, color=color, marker=marker)
    ax.vlines(x=mean[3].values, ymin=min[4], ymax=max[4], color=RED, linewidth=3)
    
    ax.hlines(xmin=[], xmax=[], y=[0, 100], color='lightgrey', linewidth=3)

    
    ax.set_xscale('log')
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.get_yaxis().set_tick_params(labelleft='off')

    
    ax.set_xlabel('nr. of links')
    
    return handles, hn, hl
         
df = pd.read_csv('output.csv', header=None, dtype={1:np.bool}, true_values=[' true'], false_values=[' false'])


print(df.groupby([0,1]).mean().sort_values([5]))
#print(df.groupby([0,1]).count())
print('max nodes', df[2].max())
print('max links', df[3].max())

directed = df.loc[~ df[1]].drop([1], axis=1)
undirected = df.loc[df[1]].drop([1], axis=1)

fig = p.figure(figsize=(10, 6))

hd, hdn, hdl = plot(directed, G1, 'd')
hu, hun, hul = plot(undirected, G2, 'o')

hd.set_label('directed')
hu.set_label('undirected')

ax = fig.axes[0]
l1 = ax.legend(handles=[hd, hu], loc='upper left', frameon=False)
ax.add_artist(l1)
l2 = ax.legend(handles=[hdn, hun], loc='lower right', frameon=False)
ax = fig.axes[1]
l3 = ax.legend(handles=[hdl, hul], loc='lower right', frameon=False)

p.tight_layout()
p.savefig('konect.pdf')

print('Maximum SEM for the total time ', df.groupby([0]).sem()[6].max())
print('Maximum SEM for the sampling ', df.groupby([0]).sem()[5].max())

