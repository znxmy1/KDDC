import pandas as pd
import scanpy as sc

shuju=sc.read_h5ad('COVID19_ALL.h5ad')
shuju
#查看某一列中的细胞类型
print(shuju.obs['majorType'].unique())
#提取两类细胞亚群
T_cell = shuju[shuju.obs['majorType'].isin(['CD4', 'CD8'])]  # 根据细胞类型提取
T_cell
#pipei
tcr=pd.read_csv('GSE158055_covid19_tcr_vdjnt_pclone.tsv',sep='\t')
tcr
#匹配，取出共有交集
T_tcr = T_cell[T_cell.obs.index.isin(tcr['cellBarcode'])]
tcr_1=tcr[tcr['cellBarcode'].isin(T_tcr.obs.index)]
tcr_1.to_csv('tcr_clean.csv')
T_tcr.write_h5ad('T_tcr.h5ad')

import pandas as pd
import scanpy as sc
T_tcr=sc.read_h5ad('T_tcr.h5ad')
T_tcr
# 确保计算邻居图
sc.pp.neighbors(T_tcr)
# 计算UMAP
sc.tl.umap(T_tcr)
sc.pl.umap(T_tcr, color='majorType', palette=['#4ABBAD', '#E585BB'], save='T_tcr_umap_majorType.pdf')
yanse = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', 
    '#7f7f7f', '#bcbd22', '#17becf', '#a0aec0', '#f6e05e', '#f56565', '#ed8936', 
    '#fc8181', '#fbbf24', '#38b2ac', '#63b3ed', '#6b46c1', '#d53f8c', '#ed64a6',
    '#9b79b0', '#f4a261', '#e9d8a6', '#e53e3e', '#dd6b20', '#f6ad55', '#fbd38d',
    '#fcf6e4', '#e2e8f0', '#d6bcfa', '#c6f6d5', '#b2f5ea'
]
sc.pl.umap(T_tcr, color='celltype', palette=yanse, save='T_tcr_umap_celltype.pdf')
sc.pl.umap(T_tcr, color='CoVID-19 severity', palette=['#80A6E2', '#FBDD85', '#F46F43'], save='T_tcr_umap_severity.pdf')
tcr_c = T_tcr[T_tcr.obs['CoVID-19 severity'] == 'control']
tcr_m = T_tcr[T_tcr.obs['CoVID-19 severity'] == 'mild/moderate']
tcr_s = T_tcr[T_tcr.obs['CoVID-19 severity'] == 'severe/critical']
# 显示 mild 的细胞
sc.pl.umap(tcr_c, color='CoVID-19 severity', palette=['#80A6E2'],title='control', save='tcr_umap_control.pdf')
# 显示 moderate 的细胞
sc.pl.umap(tcr_m, color='CoVID-19 severity', palette=['#FBDD85'],title='mild/moderate', save='tcr_umap_moderate.pdf')
# 显示 severe 的细胞
sc.pl.umap(tcr_s, color='CoVID-19 severity',palette=['#F46F43'], title='severe/critical', save='tcr_umap_severe.pdf')

# 创建一个包含状态和细胞类型的数据框
data1 = T_tcr.obs[['CoVID-19 severity', 'celltype']]
# 计算每个状态下每种细胞类型的计数
celltype_counts = data1.groupby(['CoVID-19 severity', 'celltype']).size().unstack(fill_value=0)
# 计算占比
celltype_proportions1 = celltype_counts.div(celltype_counts.sum(axis=1), axis=0)
import matplotlib.pyplot as plt
# 绘制堆叠条形图

celltype_proportions1.plot(kind='bar', stacked=True, figsize=(12, 8), color=yanse)
plt.savefig('tcr_celltype_proportions_by_severity.pdf')
celltype_proportions1.to_csv('T_tcr_bili.csv')
T_tcr.write_h5ad('T_tcr1.h5ad')


T_tcr=sc.read_h5ad('T_tcr1.h5ad')
T_tcr
tcr_clean=pd.read_csv('tcr_clean.csv')
print(tcr_clean.head())

columns_to_extract = ['CoVID-19 severity', 'celltype','majorType']
data_frame = T_tcr.obs[columns_to_extract].copy()
data_frame['cellBarcode']=data_frame.index
print(data_frame.head())
#
extract = ['cellBarcode', 'sampleID','TCRA_cgene','TCRA_vgene','TCRA_dgene','TCRA_jgene','TCRA_cdr3aa','TCRA_cdr3nt']
clean1= tcr_clean[extract].copy()
print(clean1.head())
#
tcr_he=pd.merge(data_frame,clean1)
tcr_he=tcr_he.set_index('cellBarcode')
tcr_he。to_csv('tcr_he.csv')
# #检查是否有重复值
tcr_he.index.duplicated().sum()
#
# #计算氨基酸序列的长度
aa_sequence='TCRA_cdr3aa'
tcr_he['cdr3aa_length']=tcr_he[aa_sequence].apply(len)
print(tcr_he.head())

group_length=tcr_he.groupby(['sampleID','CoVID-19 severity','cdr3aa_length']).size().reset_index(name='count')
print(group_length)
import matplotlib.pyplot as plt
import seaborn as sns
#
# ###展示不同状态下每个样本的氨基酸序列长度分布
colors=['#0094ff', '#008d00', 'orange', 'red']
g = sns.FacetGrid(group_length, col='CoVID-19 severity', col_wrap=3, height=4, aspect=1.5)
# 在每个分面中绘制不同样本的折线图
g.map_dataframe(sns.lineplot, x='cdr3aa_length', y='count', hue='CoVID-19 severity',style='sampleID', marker='o',palette=colors)
plt.show()
g.savefig("t_umap_severity_sampleID.pdf")
#
# #计算每个人不同状态下vdj基因的频率
tcr_he = tcr_he.apply(lambda col: col.cat.codes if col.dtype.name == 'category' else col)
tcr_he1=tcr_he.fillna(0)
sample_vgene=tcr_he1.groupby(['sampleID','TCRA_vgene']).size().reset_index(name='count')
sample_dgene=tcr_he1.groupby(['sampleID','TCRA_dgene']).size().reset_index(name='count')
sample_jgene=tcr_he1.groupby(['sampleID','TCRA_jgene']).size().reset_index(name='count')
vgene_matrix =sample_vgene.pivot_table(index='sampleID', columns=['TCRA_vgene'], values='count', fill_value=0)
dgene_matrix =sample_dgene.pivot_table(index='sampleID', columns=['TCRA_dgene'], values='count', fill_value=0)
jgene_matrix =sample_jgene.pivot_table(index='sampleID', columns=['TCRA_jgene'], values='count', fill_value=0)
he1=pd.merge(vgene_matrix,dgene_matrix,left_index=True, right_index=True)
he2=pd.merge(he1,jgene_matrix,left_index=True, right_index=True)

df=pd.DataFrame({'sampleID':tcr_he['sampleID'],'CoVID-19 severity':tcr_he['CoVID-19 severity']})
df1=df.drop_duplicates()
df2 = df1.set_index('sampleID')

he3 = pd.merge(he2, df2, left_index=True, right_index=True)
he3.to_csv('tcr_vdj_label.csv')



