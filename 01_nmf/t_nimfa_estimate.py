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

for m in ("H3K27ac", "H3K4me1",'H3K4me3','H3K9me3'):
	print(str(m) + ' start!')
	data = pd.read_table('../matrix/'+m+'_variable_peak_matrix.tsv',sep = '\t', header=0, index_col=0)
	data.columns = data.columns.str[:-5]
	data_avg = data.groupby(data.columns, axis=1).mean()
	data = np.array(data_avg).T
	data = preprocessing.Normalizer().fit_transform(data)


	rank_cands = range(2, 6, 1)
	snmf = nimfa.Nmf(data, seed='nndsvd', max_iter=100)
	summary = snmf.estimate_rank(rank_range=rank_cands, n_run=10, what='all')
	
	rss = [summary[rank]['rss'] for rank in rank_cands]
	coph = [summary[rank]['cophenetic'] for rank in rank_cands]
	disp = [summary[rank]['dispersion'] for rank in rank_cands]
	spar = [summary[rank]['sparseness'] for rank in rank_cands]
	spar_w, spar_h = zip(*spar)
	evar = [summary[rank]['evar'] for rank in rank_cands]

	# 创建多个子图
	fig, axes = plt.subplots(3, 2, figsize=(12, 12))
	fig.tight_layout(pad=5.0)

	# 绘制 RSS
	axes[0, 0].plot(rank_cands, rss, 'o-', color='blue', linewidth=2)
	axes[0, 0].set_title('RSS')
	axes[0, 0].set_xlabel('Rank')
	axes[0, 0].set_ylabel('RSS Value')
	axes[0, 0].grid(True)

	# 绘制 Cophenetic correlation
	axes[0, 1].plot(rank_cands, coph, 'o-', color='green', linewidth=2)
	axes[0, 1].set_title('Cophenetic Correlation')
	axes[0, 1].set_xlabel('Rank')
	axes[0, 1].set_ylabel('Cophenetic Value')
	axes[0, 1].grid(True)

	# 绘制 Dispersion
	axes[1, 0].plot(rank_cands, disp, 'o-', color='red', linewidth=2)
	axes[1, 0].set_title('Dispersion')
	axes[1, 0].set_xlabel('Rank')
	axes[1, 0].set_ylabel('Dispersion Value')
	axes[1, 0].grid(True)

	# 绘制 Sparsity (Basis)
	axes[1, 1].plot(rank_cands, spar_w, 'o-', color='purple', linewidth=2)
	axes[1, 1].set_title('Sparsity (Basis)')
	axes[1, 1].set_xlabel('Rank')
	axes[1, 1].set_ylabel('Sparsity Value')
	axes[1, 1].grid(True)

	# 绘制 Sparsity (Mixture)
	axes[2, 0].plot(rank_cands, spar_h, 'o-', color='orange', linewidth=2)
	axes[2, 0].set_title('Sparsity (Mixture)')
	axes[2, 0].set_xlabel('Rank')
	axes[2, 0].set_ylabel('Sparsity Value')
	axes[2, 0].grid(True)
	
	# 绘制 Explained variance
	axes[2, 1].plot(rank_cands, evar, 'o-', color='brown', linewidth=2)
	axes[2, 1].set_title('Explained Variance')
	axes[2, 1].set_xlabel('Rank')
	axes[2, 1].set_ylabel('Explained Variance Value')
	axes[2, 1].grid(True)
	
	# 显示图像
	plt.savefig('t/result/' + str(m) + '/rank_estimate_2-6.pdf', dpi=300, bbox_inches='tight')
	#plt.show()

	with open('t/result/'+str(m)+'/estimate_2-6.pkl', 'wb') as f:
		 pickle.dump(snmf, f)
	
	print(str(m) + ' finish!')


