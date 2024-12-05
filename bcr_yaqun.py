import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import pandas as pd
import scanpy as sc

shuju=sc.read_h5ad('COVID19_ALL.h5ad')
shuju
#查看某一列中的细胞类型
print(shuju.obs['majorType'].unique())
#提取一类细胞亚群
B_cell = shuju[shuju.obs['majorType'] == 'B']  # 提取标注为 "T_cells" 的细胞
B_cell

bcr=pd.read_csv('GSE158055_covid19_bcr_vdjnt_pclone.tsv',sep='\t')
bcr
#匹配，取出共有交集
B_bcr = B_cell[B_cell.obs.index.isin(bcr['cellBarcode'])]
bcr_1=bcr[bcr['cellBarcode'].isin(B_bcr.obs.index)]
bcr_1.to_csv('bcr_clean.csv')
B_bcr.write_h5ad('B_bcr.h5ad')

B_bcr=sc.read_h5ad('B_bcr.h5ad')
B_bcr
# 确保计算邻居图
sc.pp.neighbors(B_bcr)
# 计算UMAP
sc.tl.umap(B_bcr)
sc.pl.umap(B_bcr, color='celltype', palette=['#0094ff', '#008d00', 'orange', 'red'], save='B_bcr_umap_celltype.pdf')
sc.pl.umap(B_bcr, color='CoVID-19 severity', palette=['#80A6E2', '#FBDD85', '#F46F43'], save='B_bcr_umap_severity.pdf')
bcr_c = B_bcr[B_bcr.obs['CoVID-19 severity'] == 'control']
bcr_m = B_bcr[B_bcr.obs['CoVID-19 severity'] == 'mild/moderate']
bcr_s = B_bcr[B_bcr.obs['CoVID-19 severity'] == 'severe/critical']
# 显示 mild 的细胞
sc.pl.umap(bcr_c, color='CoVID-19 severity', palette=['#80A6E2'],title='control', save='bcr_umap_control.pdf')
# 显示 moderate 的细胞
sc.pl.umap(bcr_m, color='CoVID-19 severity', palette=['#FBDD85'],title='mild/moderate', save='bcr_umap_moderate.pdf')
# 显示 severe 的细胞
sc.pl.umap(bcr_s, color='CoVID-19 severity',palette=['#F46F43'], title='severe/critical', save='bcr_umap_severe.pdf')

# 创建一个包含状态和细胞类型的数据框
data = B_bcr.obs[['CoVID-19 severity', 'celltype']]
# 计算每个状态下每种细胞类型的计数
celltype_counts = data.groupby(['CoVID-19 severity', 'celltype']).size().unstack(fill_value=0)
# 计算占比
celltype_proportions = celltype_counts.div(celltype_counts.sum(axis=1), axis=0)
import matplotlib.pyplot as plt
# 绘制堆叠条形图
colors=['#0094ff', '#008d00', 'orange', 'red']
celltype_proportions.plot(kind='bar', stacked=True, figsize=(12, 8), color=colors)
plt.savefig('celltype_proportions_by_severity.pdf')
celltype_proportions.to_csv('B_bcr_bili.csv')
B_bcr.write_h5ad('bcr1.h5ad')

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

bcr=sc.read_h5ad('bcr1.h5ad')
bcr
sc.tl.umap(bcr)
sc.pl.umap(bcr)
sc.tl.leiden(bcr,resolution=0.4)
sc.pl.umap(bcr,color='leiden',show=False)
plt.savefig('bcr_cluster.pdf', bbox_inches='tight')
plt.close()
sc.tl.rank_genes_groups(bcr, "leiden", method="t-test")
sc.pl.rank_genes_groups(bcr, n_genes=10, sharey=False,show=False)
plt.savefig('bcr_topgene.pdf', bbox_inches='tight')
plt.close()

result = bcr.uns['rank_genes_groups']
groups = result["names"].dtype.names
topgene=pd.DataFrame(
     {
         group + "_" + key[:1]: result[key][group]
         for group in groups
         for key in ["names"]
  }
 ).head(50)
topgene.to_csv('bcr_topgene.csv')


marker_genes=['CD79A', 'CD79B', 'MS4A1','MZB1', 'XBP1']
sc.pl.violin(bcr, marker_genes, groupby="leiden",stripplot=False,show=False)
plt.savefig('bcr_canonicalgenes.pdf',dpi=300,bbox_inches='tight')
plt.close()

#0
markers=['FOS','CXCR4','BTG1','IGHD','BACH2','IL4R','CD69','TCL1A','FOSB','JUND','DUSP1']
sc.pl.violin(bcr, markers, groupby="leiden",stripplot=False)
plt.savefig('bcr_canonicalgenes.pdf',dpi=300,bbox_inches='tight')
plt.close()
#1
markers=['CD37','CD79B','TCL1A','FCER2','IGHV7-4-1','IGHV5-10-1','TNFSF13B','IGLL5']
sc.pl.violin(bcr, markers, groupby="leiden",stripplot=False)
#2
markers=['TCL1A','IL4R','FCER2','IGHD']
sc.pl.violin(bcr, markers, groupby="leiden",stripplot=False)

markers=['HLA-DQA2','HLA-DRB5','HLA-DRB1','IGHV4-34','IGHV7-4-1']
sc.pl.dotplot(bcr, markers, groupby="leiden")
#
markers=['CD27','TCL1A','MS4A1','AIM2','SOX5','TNFRSF1B']
markers=['CD27','TCL1A']
sc.pl.violin(bcr, markers, groupby="leiden",stripplot=False)

markers=['CXCR4','TCL1A','CD27','HLA-B','HLA-DQA2','LY6E','TNFSF13B']

markers=['TCL1A']
markers=['CD27']
markers=['MS4A1','AIM2','SOX5','TNFRSF1B']
sc.pl.violin(bcr,markers,groupby='leiden',stripplot=False)

new_cluster_names = {'0':"B_c01_TCL1A",#0
    '1':"B_c01_TCL1A",#1
    '2':"B_c02_CD27",#2
    '3':"B_c01_TCL1A",#3
    '4':"B_c02_CD27",#4
    '5':"B_c02_CD27",#5
    '6':"B_c023_TNFRSF1B",#6
    '7':"B_c01_TCL1A",#7
    '8':"B_c01_TCL1A",#8
    '9':"B_c01_TCL1A"#9
                     }
bcr.obs['cate'] = bcr.obs['leiden'].replace(new_cluster_names)
sc.pl.umap(bcr,color='cate',show=False)
plt.savefig('bcr_cell.pdf',dpi=300,bbox_inches='tight')
plt.close()

#输出标记基因图
sc.pl.violin(bcr,'TCL1A',groupby='cate',stripplot=False,show=False,rotation=15)
plt.savefig('bcr_TCL1A.pdf',dpi=300,bbox_inches='tight')
plt.close()

sc.pl.violin(bcr,'CD27',groupby='cate',stripplot=False,show=False,rotation=15)
plt.savefig('bcr_CD27.pdf',dpi=300,bbox_inches='tight')
plt.close()

sc.pl.violin(bcr,'TNFRSF1B',groupby='cate',stripplot=False,show=False,rotation=15)
plt.savefig('bcr_TNFRSF1B.pdf',dpi=300,bbox_inches='tight')
plt.close()
#细胞比例
data = bcr.obs[['CoVID-19 severity', 'cate']]
# 计算每个状态下每种细胞类型的计数
celltype_counts = data.groupby(['CoVID-19 severity', 'cate']).size().unstack(fill_value=0)
# 计算占比
celltype_proportions = celltype_counts.div(celltype_counts.sum(axis=1), axis=0)
celltype_proportions.to_csv('bcr_ratio.csv')

# 绘制堆叠条形图
fig,ax=plt.subplots(figsize=(8,6))
celltype_proportions.plot(kind='bar',stacked=True,ax=ax)
ax.set_xticklabels(ax.get_xticklabels(),rotation=15,ha='right')
plt.legend(loc='right',bbox_to_anchor=(1.3,0.8),borderaxespad=0)
plt.savefig('bcr_ratio.pdf',bbox_inches='tight')
#####计算基因集得分
from sklearn.metrics import roc_auc_score
#耗竭
haojie=["PDCD1","CTLA4","TIGIT","TCF7","NR4A1","TOX","TOX2","NFAT5","LAG3","HAVCR2","NR4A2"]
cuyan=["TNF","IL6","IL18","IL1B","IL2","IL7", "CCL2", "CXCL8"]
kangyan=["IL10","IL11","IL22","IL37","IL4","IL5","IL13"]
huohua=["TYROBP","MYD88","IKBKB","TNF","NFKB1","NOS2","IRF5","NFKB2","IL1B","IL1A","TNFSF11","SLAMF1","LAG3"]

gene_expression=bcr[:,haojie].X
scores=np.mean(gene_expression,axis=1)
bcr.obs['haojie_score']=score

from scipy import stats

b1=bcr[bcr.obs['cate']=='B_c01_TCL1A']
grouped_data = b1.obs.groupby('CoVID-19 severity')['gene_score'].apply(list)
kruskal_result = stats.kruskal(*grouped_data)
print(f"Kruskal-Wallis H 检验统计量: {kruskal_result.statistic}, p值: {kruskal_result.pvalue}")
sc.pl.violin(b1,'gene_score',groupby='CoVID-19 severity',stripplot=False,rotation=15,show=False)
plt.text(0.5, 0.95, f"Kruskal-Wallis H 检验: p = {kruskal_result.pvalue:.4f}",
         ha='center', va='top', transform=plt.gca().transAxes)
plt.show()

b2=bcr[bcr.obs['cate']=='B_c02_CD27']
rouped_data = b2.obs.groupby('CoVID-19 severity')['gene_score'].apply(list)
kruskal_result = stats.kruskal(*grouped_data)
print(f"Kruskal-Wallis H 检验统计量: {kruskal_result.statistic}, p值: {kruskal_result.pvalue}")
sc.pl.violin(b2,'gene_score',groupby='CoVID-19 severity',stripplot=False,rotation=15,show=False)
plt.text(0.5, 0.95, f"Kruskal-Wallis H 检验: p = {kruskal_result.pvalue:.4f}",
         ha='center', va='top', transform=plt.gca().transAxes)
plt.show()

b3=bcr[bcr.obs['cate']=='B_c03_TNFRSF1B']
sc.pl.violin(b3,haojie,groupby='CoVID-19 severity',stripplot=False,rotation=15)
rouped_data = b3.obs.groupby('CoVID-19 severity')['gene_score'].apply(list)
kruskal_result = stats.kruskal(*grouped_data)
print(f"Kruskal-Wallis H 检验统计量: {kruskal_result.statistic}, p值: {kruskal_result.pvalue}")
sc.pl.violin(b3,'gene_score',groupby='CoVID-19 severity',stripplot=False,rotation=15,show=False)
plt.text(0.5, 0.95, f"Kruskal-Wallis H 检验: p = {kruskal_result.pvalue:.4f}",
         ha='center', va='top', transform=plt.gca().transAxes)
plt.show()

sc.pl.violin(b1,cuyan,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,cuyan,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,cuyan,groupby='CoVID-19 severity',stripplot=False,rotation=15)

kangyan=["IL10","IL11","IL22","IL37","IL4","IL5","IL13"]
sc.pl.violin(b1,kangyan,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,kangyan,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,kangyan,groupby='CoVID-19 severity',stripplot=False,rotation=15)

huohua=["TYROBP","MYD88","IKBKB","TNF","NFKB1","NOS2","IRF5","NFKB2","IL1B","IL1A","TNFSF11","SLAMF1","LAG3"]
sc.pl.violin(b1,huohua,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,huohua,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,huohua,groupby='CoVID-19 severity',stripplot=False,rotation=15)

###差异基因的比较
selected_statuses=['control','mild/moderate']

b1_subset = b1[b1.obs['CoVID-19 severity'].isin(selected_statuses)]
# b1_subset.obs['selected_status'] = b1_subset.obs['status']
sc.tl.rank_genes_groups(b1_subset, groupby='CoVID-19 severity', method='t-test', n_genes=100)
sc.pl.rank_genes_groups(b1_subset, n_genes=100, sharey=False)
result = b1_subset.uns['rank_genes_groups']
groups = result['names'].dtype.names

df_results = pd.DataFrame()

# 提取每组的基因名、p值等信息
for group in groups:
    df_results[group] = result['names'][group]  # 基因名
    df_results[f'{group}_pval'] = result['pvals'][group]  # p值
    df_results[f'{group}_pval_adj'] = result['pvals_adj'][group]  # 校正后的 p值
    df_results[f'{group}_logfoldchanges'] = result['logfoldchanges'][group]  # 对数倍数变化

# 重置索引
df_results.reset_index(drop=True, inplace=True)

df_results.to_csv('differential_expression_results.csv', index=False)



b3_subset = b3[b3.obs['CoVID-19 severity'].isin(selected_statuses)]
# b1_subset.obs['selected_status'] = b1_subset.obs['status']
sc.tl.rank_genes_groups(b3_subset, groupby='CoVID-19 severity', method='t-test', n_genes=100)
sc.pl.rank_genes_groups(b3_subset, n_genes=100, sharey=False)
result = b3_subset.uns['rank_genes_groups']
groups = result['names'].dtype.names
df_results = pd.DataFrame()

# 提取每组的基因名、p值等信息
for group in groups:
    df_results[group] = result['names'][group]  # 基因名
    df_results[f'{group}_pval'] = result['pvals'][group]  # p值
    df_results[f'{group}_pval_adj'] = result['pvals_adj'][group]  # 校正后的 p值
    df_results[f'{group}_logfoldchanges'] = result['logfoldchanges'][group]  # 对数倍数变化

# 重置索引
df_results.reset_index(drop=True, inplace=True)
df_results.to_csv('b3_differential_expression_results.csv', index=False)


#tonglu
mapk=['JUN', 'JUND', 'GADD45B', 'DUSP1', 'NF1', 'PDGFA', 'FOS', 'RAC1', 'FGF23']
sc.pl.violin(b1,mapk,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,mapk,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,mapk,groupby='CoVID-19 severity',stripplot=False,rotation=15)

Oxide=['ATP5F1A', 'MT-ND4', 'MT-CO2', 'MT-ND2', 'MT-ND3', 'MT-ND1']
sc.pl.violin(b1,Oxide,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,Oxide,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,Oxide,groupby='CoVID-19 severity',stripplot=False,rotation=15)

tgf=['SKI', 'CREBBP', 'ID3', 'SKIL', 'TNF']
sc.pl.violin(b1,tgf,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,tgf,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,tgf,groupby='CoVID-19 severity',stripplot=False,rotation=15)

il17=['JUN', 'JUND', 'FOSB', 'FOS']
sc.pl.violin(b1,il17,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b2,il17,groupby='CoVID-19 severity',stripplot=False,rotation=15)
sc.pl.violin(b3,il17,groupby='CoVID-19 severity',stripplot=False,rotation=15)

####yaqun

b1=bcr[bcr.obs['cate']=='B_c01_TCL1A']
b2=bcr[bcr.obs['cate']=='B_c02_CD27']
b3=bcr[bcr.obs['cate']=='B_c03_TNFRSF1B']

len(b1.obs['sampleID'].value_counts())
b1_1=b1.obs[['CoVID-19 severity','sampleID']].value_counts()
b1_1.groupby('CoVID-19 severity').size()

len(b2.obs['sampleID'].value_counts())
b2_1=b2.obs[['CoVID-19 severity','sampleID']].value_counts()
b2_1.groupby('CoVID-19 severity').size()

len(b3.obs['sampleID'].value_counts())
b3_1=b3.obs[['CoVID-19 severity','sampleID']].value_counts()
b3_1.groupby('CoVID-19 severity').size()

# tgf=['SKI', 'CREBBP', 'ID3', 'SKIL', 'TNF','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','TGFBR3','SMAD2','SMAD3','SMAD4','SMAD5','SMAD6','SMAD7']
bcr.obs['type_state']=bcr.obs['cate'].astype(str)+'_'+bcr.obs['CoVID-19 severity'].astype(str)
# sc.pl.dotplot(bcr,tgf,groupby='type_state')

il17=['JUN', 'JUND', 'FOSB', 'FOS']
sc.pl.dotplot(bcr,il17,groupby='type_state')

tgf=['SKI', 'CREBBP', 'ID3', 'SKIL', 'TNF','TGFB1','SMAD7']
sc.pl.dotplot(bcr,tgf,groupby='type_state')

Oxide=['ATP5F1A', 'MT-ND4', 'MT-CO2', 'MT-ND2', 'MT-ND3', 'MT-ND1']
sc.pl.dotplot(bcr,Oxide,groupby='type_state')
mapk=['JUN', 'JUND', 'GADD45B', 'DUSP1', 'NF1', 'PDGFA', 'FOS', 'RAC1', 'FGF23']
sc.pl.dotplot(bcr,mapk,groupby='type_state')


haojie=["PDCD1","CTLA4","TIGIT","TCF7","NR4A1","TOX","TOX2","NFAT5","LAG3","HAVCR2","NR4A2"]
sc.pl.dotplot(bcr,haojie,groupby='type_state')
cuyan=["TNF","IL6","IL18","IL1B","IL2","IL7", "CCL2", "CXCL8"]
sc.pl.dotplot(bcr,cuyan,groupby='type_state')
kangyan=["IL10","IL11","IL22","IL37","IL4","IL5","IL13"]
sc.pl.dotplot(bcr,kangyan,groupby='type_state')
huohua=["TYROBP","MYD88","IKBKB","TNF","NFKB1","NOS2","IRF5","NFKB2","IL1B","IL1A","TNFSF11","SLAMF1","LAG3"]
sc.pl.dotplot(bcr,huohua,groupby='type_state')

genes=['MT-CO2','MT-ND3','JUN','DUSP1','FOS','FGF23',"NFKB1","NFKB2",'ID3','SKIL','TNF']
sc.pl.dotplot(bcr,genes,groupby='type_state',show=False,swap_axes=True)
plt.savefig('bcr_genes.pdf',dpi=300,bbox_inches='tight')
plt.close()

markers=['TCL1A','CD27','TNFRSF1B']
sc.pl.heatmap(bcr,markers,groupby='cate')
sc.pl.stacked_violin(bcr,markers,groupby='cate',stripplot=False,show=False)
plt.savefig('bcr_markers.pdf',dpi=300,bbox_inches='tight')
plt.close()



