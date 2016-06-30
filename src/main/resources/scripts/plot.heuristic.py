import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as p
import numpy as n
import pylab
import scipy.stats as stats
import networkx as nwx
import glob
import builtins

def clean(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.get_xaxis().set_tick_params(which='both', top='off')
    ax.get_yaxis().set_tick_params(which='both', right='off')
    ax.axhline(color='k')

data = n.genfromtxt('out.csv', delimiter=',')
exDegrees = data[:,0]
numNodes = data[:,1]
numLinks = data[:,2]
profits = data[:,3]

alpha = 0.7
size = 20

ax = p.subplot(311)
p.scatter(numNodes, profits, c=exDegrees, s=size, linewidth=0, alpha=alpha)
ax.text(.5,.9,'number of nodes (color=)',horizontalalignment='center',
        transform=ax.transAxes)
clean(ax)
p.colorbar()

ax = p.subplot(312)
p.scatter(numLinks, profits, c=exDegrees, s=size, linewidth=0, alpha=alpha)
ax.text(.5,.9,'number of links (color=exdegree)',horizontalalignment='center',
        transform=ax.transAxes)
clean(ax)
p.colorbar()

ax = p.subplot(313) 
p.scatter(exDegrees, profits, c=numLinks, s=size, linewidth=0, alpha=alpha)
ax.text(.5,.9,'exdegree (color=num links)',horizontalalignment='center',
        transform=ax.transAxes)
clean(ax)
p.colorbar()

p.savefig('heuristic.png')

p.figure()

data = n.genfromtxt('cutoffs.csv', delimiter=',')
p.scatter(data[:,0], data[:,1], linewidth=0, alpha=0.05)

p.ylim([0, 500])
p.savefig('cutoffs.png')
