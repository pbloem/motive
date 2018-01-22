# -*- coding: utf-8 -*-

import matplotlib as mpl
from _socket import NI_DGRAM
mpl.use('Agg')
import matplotlib.pyplot as p
import numpy as n
import pylab
import scipy.stats as stats
import networkx as nwx
import glob
import builtins
from matplotlib.pyplot import margins
import os.path
import json

RED = 'darkred'
G1 = 'lightgrey'
G2 = 'silver'
G3 = 'darkgrey'

mpl.style.use('classic')

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

mpl.rc('font', **font)

margin = 0.05
extra = 0.05

row1height = 0.6
row2height = 0.2
row3height = 0.2

# To be run in the workspace dir of the UCompareBeta module
barwidth = 0.9
pluswidth = 0.45

# Load experiment metadata
with open('metadata.json') as mdfile:
    metadata = json.load(mdfile)

sub_index = metadata["subindex"]
nums_instances = metadata["nums instances"]
motif_size = metadata["motif size"]
directed = False

ni = len(nums_instances)

# Load the frequencies and factors
frequencies = n.genfromtxt('frequencies.csv', delimiter=',')
(nummotifs, width) = frequencies.shape
factors = n.genfromtxt('factors.csv', delimiter=',')
means = n.genfromtxt('means.csv', delimiter=',')

runs = width // ni

freqMeans = means[:,0:ni]
factMeans = means[:,ni:2*ni]

fig = p.figure(figsize=(16,7))

### 1) Plot the factors
ax1 = fig.add_axes([0.0 + margin + extra, row3height + row2height + margin, 1.0 - 2.0 * margin- extra, row1height - 2.0 * margin]); 

ind = n.arange(nummotifs)

bw = barwidth/ni
for i in range(ni):
    color = G1
    label = u'$k = 0$'
    if i == 1:
        color = G2
        label = u'$k = 10$'
    if i == 2:
        color = G3
        label = u'$k = 100$'
    
    # the means as bars
    print(color)
    bars = ax1.bar(ind - barwidth/2.0 + i * bw, factMeans[:, i], bw, color=color, zorder=1, linewidth=0)
    bars.set_label(label)
    
for i in range(ni):
    # the data as scatter
    for s in range(nummotifs):
        
        min = n.min(factors[s,i*runs:(i+1)*runs])
        max = n.max(factors[s,i*runs:(i+1)*runs])
        ax1.vlines((ind[s] - barwidth/2.0 + (i+0.5) * bw),min, max, colors=RED, linewidths=2, zorder=3)
    
ax1.set_xlim([0 - pluswidth, nummotifs - 1 + pluswidth])

ax1.hlines(0, - pluswidth, nummotifs - 1 + pluswidth)

yloc = p.MaxNLocator(7)
ax1.get_yaxis().set_major_locator(yloc)

ax1.get_yaxis().set_tick_params(which='both', direction='out')

ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)
ax1.spines["bottom"].set_visible(False)
ax1.spines["left"].set_visible(False)

ax1.get_xaxis().set_tick_params(which='both', top='off', bottom='off', labelbottom='off')
ax1.get_yaxis().set_tick_params(which='both', left='off', right='off')

# top = n.max(factor)
# if n.min(factor) < - top and top > 0:
#   ax1.set_ylim(bottom=-top)
   
# negative grid (white lines over the bars)   
ticks = ax1.get_yaxis().get_majorticklocs()   
ticks = n.delete(ticks, n.where(n.logical_and(ticks < 0.00001, ticks > -0.00001)))
ax1.hlines(ticks, - pluswidth, nummotifs - 1 + pluswidth, color='w', zorder=2)

ax1.legend()
ax1.set_ylabel('log-factor (bits)')

### 2) Plot the motifs

bottom = margin
height = row2height - margin

side = pluswidth - 0.5
width = (1.0 - 2.0 * margin - extra) / (nummotifs + 2.0 * side)

i = 0
for path in sorted(glob.glob('motif.*.edgelist'))[:nummotifs]:
    axsmall = fig.add_axes([margin + extra + side*width + width * i, bottom, width, height])
    axsmall.axis('off')
    
    graph = nwx.read_edgelist(path,create_using=(nwx.DiGraph() if directed else nwx.Graph()))
    ng = nwx.number_of_nodes(graph)
    
    pos = nwx.spring_layout(graph)
    nodes = nwx.draw_networkx_nodes(graph, pos, ax=axsmall, node_size=12)
    if nodes != None:
        nodes.set_edgecolor(RED)
        nodes.set_color(RED)
    color = RED if i == sub_index else 'k'
    edges = nwx.draw_networkx_edges(graph, pos, alpha=0 if directed else 1, fc=color, edge_color=color)
    if nodes == None or ng < motif_size:
        (minx, maxx) = axsmall.get_xlim()
        ran = maxx - minx
        rem = motif_size if (nodes == None) else motif_size - ng
        axsmall.scatter((n.arange(rem) * (0.333/rem) + 0.666) * ran + minx, 0 * n.ones(rem), s=12, color=RED)   

    i = i + 1
    
### 3)  Frequency graph

ax3 = fig.add_axes([0.0 + margin + extra, row2height + margin, 1.0 - 2.0 * margin - extra, row3height - margin]) 

# ax3.bar(ind - barwidth/2.0, freq, barwidth, color='k')
for i in range(ni):
    color = G1
    if i == 1:
        color = G2
    if i == 2:
        color = G3
    
    # the means as bars
    ax3.bar(ind - barwidth/2.0 + i * bw, freqMeans[:, i], bw, color=color, zorder=1, linewidth=0)
# the data as scatter
for i in range(ni):
    for s in range(nummotifs):
        min = n.min(frequencies[s,i*runs:(i+1)*runs])
        max = n.max(frequencies[s,i*runs:(i+1)*runs])
        ax3.vlines((ind[s] - barwidth/2.0 + (i+0.5) * bw),min, max, colors=RED, linewidths=2, zorder=3)
        
ax3.get_yaxis().set_tick_params(which='both', direction='out')

ax3.set_xlim([0 - pluswidth, nummotifs - 1 + pluswidth])

# reduce the number of ticks
yloc = p.MaxNLocator(4)
ax3.yaxis.set_major_locator(yloc)

ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)
ax3.spines["left"].set_visible(False)

ax3.get_xaxis().tick_bottom()
ax3.get_xaxis().set_tick_params(which='both', top='off', bottom='off', right='off', labelbottom='off')
ax3.get_yaxis().set_tick_params(which='both', left='off', right='off')

ax3.set_ylim([0, ax3.get_ylim()[1]])

ticks = ax3.get_yaxis().get_majorticklocs()   
ticks = n.delete(ticks, n.where(n.logical_and(ticks < 0.00001, ticks > -0.00001)))
ax3.hlines(ticks, - pluswidth, nummotifs - 1 + pluswidth, color='w', zorder=2)

ax3.set_ylabel('freq.')

p.savefig('synthetic-plot.png')
p.savefig('synthetic-plot.pdf')
