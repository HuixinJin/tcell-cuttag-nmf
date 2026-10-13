import pandas as pd
import pickle

# H3K9me3_rank3

peaks = pd.read_table('../../matrix/H3K9me3_variable_peak_matrix.tsv', sep = '\\t')
print('peaks shape:\n')
print(peaks.shape)

cluster = pd.read_table('result/H3K9me3/cluster_assignments.txt')
print('cluster shape:\n')
print(cluster.shape)

cluster['chr'] = [index.split(':')[0] for index in peaks.index.values]
cluster['start'] = [index.split(':')[1].split('-')[0] for index in peaks.index.values]
cluster['end'] = [index.split(':')[1].split('-')[1] for index in peaks.index.values]
cluster.head()

cluster['Cluster'].value_counts()

for i in set(cluster['Cluster']):
	cluster_df = cluster.loc[cluster['Cluster']==i,['chr','start','end']]
	cluster_df.to_csv('result/H3K9me3/cluster_bed/cluster_' + str(i) + '_peaks.bed', index=None, header=None, sep = '\t')
