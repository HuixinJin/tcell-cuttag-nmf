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
#sys.setrecursionlimit(13)

#rank = 3
#data = pd.read_table('../matrix/H3K9me3_variable_peak_matrix.tsv', sep = '\t', header=0, index_col=0)
#data.columns = data.columns.str[:-5]
#data_avg = data.groupby(data.columns, axis=1).mean()
#data = np.array(data_avg)
#data = preprocessing.Normalizer().fit_transform(data)

#nmf = nimfa.Nmf(data.T, rank=rank, seed='nndsvd', max_iter=100, update='divergence', objective='div', n_run=50, track_factor=True)
#nmf_fit = nmf()
#C = 1 - nmf_fit.fit.consensus()
#with open('t_rank3/result/H3K9me3/C.pkl', 'wb') as f:
	#C = pickle.load(f)
	#pickle.dump(C, f)
with open('t/result/H3K9me3/C.pkl', 'rb') as f:
	C = pickle.load(f)
#Y = sch.linkage(C, method='average')
with open('t/result/H3K9me3/Y_linkage.pkl', 'rb') as f:
	Y = pickle.load(f)

print(f"Y shape: {Y.shape}")  
print(f"Y data type: {Y.dtype}")  
print(f"First few rows of Y:\n{Y[:5]}")  


from scipy.cluster.hierarchy import leaves_list  
leaf_order = leaves_list(Y)  
print(f"Leaf order shape: {leaf_order.shape}")  

 
import matplotlib.pyplot as plt  
import seaborn as sns  
import numpy as np  
# 获取重排后的距离矩阵  
C_ordered = C[leaf_order, :][:, leaf_order]  

# 创建聚类热图  
# 降低热图渲染质量  

plt.figure(figsize=(12, 10), dpi=72)  # 降低DPI  
plt.rcParams['figure.facecolor'] = 'w'  
plt.rcParams['savefig.dpi'] = 72  
plt.rcParams['savefig.bbox'] = 'tight'  
plt.rcParams['savefig.pad_inches'] = 0.1  
sns.heatmap(C_ordered, cmap='coolwarm_r',   
            xticklabels=False, yticklabels=False,  
            rasterized=True)  # 将矢量图转为栅格图  
plt.title('NMF Consensus Distance Matrix (Clustered)', fontsize=14)  
plt.tight_layout()  

# 使用较低的文件质量保存  
plt.savefig('t/result/H3K9me3/consensus_heatmap.png', dpi=72)  
plt.close()  

#denD = sch.dendrogram(Y, orientation='left', link_color_func=lambda k: 'black', no_plot=True)
	#truncate_mode='level',  # 添加树状图截断  
	#p=5)  # 只显示前5层（可以根据需要调整）

