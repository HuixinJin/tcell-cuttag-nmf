from collections import defaultdict, Counter
import urllib
import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch

import nimfa
import fastcluster
from scipy.cluster.hierarchy import fcluster
import sys
sys.setrecursionlimit(100000)

for m in (['H3K9me3']):
#for m in ("H3K27ac","H3K4me1",'H3K4me3','H3K9me3'):
	#if m in ("H3K27ac", "H3K4me1"):
		#rank = 15
	#else:
		#rank = 9
	rank = 3 
	data = pd.read_table('../matrix/' + m + '_variable_peak_matrix.tsv',
        	             sep = '\t', header=0, index_col=0)
	data.columns = data.columns.str[:-5]
	data_avg = data.groupby(data.columns, axis=1).mean()
	data = np.array(data_avg)
	data = preprocessing.Normalizer().fit_transform(data)


	nmf = nimfa.Nmf(data.T, rank=rank, seed='nndsvd', max_iter=100, 
        	        update='divergence', objective='div', n_run=50, track_factor=True)
	nmf_fit = nmf()

	def clean_axis(ax):
	    ax.get_xaxis().set_ticks([])
	    ax.get_yaxis().set_ticks([])
	    for sp in ax.spines.values():
	        sp.set_visible(False)

	#fig = plt.figure(figsize=(13.9, 10))
	#heatmapGS = gridspec.GridSpec(1, 2, wspace=.01, hspace=0., width_ratios=[0.25,1])

	C = 1 - nmf_fit.fit.consensus()
	with open('t_rank3/result/'+str(m)+'/C.pkl', 'wb') as f:
		pickle.dump(C, f)
	Y = sch.linkage(C, method='average')
	#Y = fastcluster.linkage(C, method='average')
	with open('t_rank3/result/'+str(m)+'/Y_linkage.pkl', 'wb') as f:
		pickle.dump(Y, f)

	#denAX = fig.add_subplot(heatmapGS[0,0])


	#denD = sch.dendrogram(Y, orientation='left', link_color_func=lambda k: 'black')
	#clean_axis(denAX)

	#heatmapAX = fig.add_subplot(heatmapGS[0,1])
	#D = C[denD['leaves'], :][:, denD['leaves']]
	#axi = heatmapAX.imshow(D, interpolation='nearest', aspect='equal', origin='lower', cmap='RdBu') 
	#clean_axis(heatmapAX)

	#cb = fig.colorbar(axi, fraction=0.046, pad=0.04, aspect=10) 
	#cb.set_label('Distance', fontsize=20)

	#plt.savefig("t/result/" +str(m)+"/consensus_heatmap.pdf", dpi = 300, bbox_inches='tight')


	# Extract flat clusters at a specific height
	clusters = fcluster(Y, t=0.5, criterion='distance')

	# Save cluster assignments to a text file
	np.savetxt('t_rank3/result/' + str(m) + '/cluster_assignments.txt', clusters, fmt='%d', header='Cluster', comments='')

