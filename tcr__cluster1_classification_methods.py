import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score
from sklearn.feature_selection import f_classif
from sklearn.model_selection import train_test_split

# kmer2
# cluster_kmer = pd.read_csv('tcr_cluster1_kmer_2.csv')
# cluster_kmer = pd.read_csv('tcr_cluster1_kmer_3.csv')
cluster_kmer = pd.read_csv('tcr_cluster1_kmer_4.csv')
cluster_kmer = cluster_kmer.set_index('sample')
cluster_kmer_1 = cluster_kmer.T
cluster_kmer_1 = cluster_kmer_1.div(cluster_kmer_1.sum(axis=1), axis=0)
cluster_kmer_1['class'] = cluster_kmer_1.index
cluster_kmer_1['state'] = cluster_kmer_1['class'].str.split('-').str[1].str[0]
cluster_kmer_1 = cluster_kmer_1.set_index('class')
x = cluster_kmer_1.drop(columns=['state'])
zero_count = (x == 0).sum(axis=0)
y = cluster_kmer_1['state']
y_class=pd.factorize(y)[0]
print(y)
print(y_class)
#kmer2_kmer3
# cluster_kmer=pd.read_csv('tcr_cluster1_kmer_2.csv')
# cluster_kmer2 = pd.read_csv('tcr_cluster1_kmer_3.csv')
# cluster_kmer=cluster_kmer.set_index('sample')
# cluster_kmer_1=cluster_kmer.T
# cluster_kmer_1=cluster_kmer_1.div(cluster_kmer_1.sum(axis=1),axis=0)
# cluster_kmer_1['class']=cluster_kmer_1.index
# cluster_kmer2= cluster_kmer2.set_index('sample')
# cluster_kmer2_1 = cluster_kmer2.T
# cluster_kmer2_1 = cluster_kmer2_1.div(cluster_kmer2_1.sum(axis=1), axis=0)
# cluster_kmer2_1['class'] = cluster_kmer2_1.index
#
# he1 = pd.merge(cluster_kmer_1, cluster_kmer2_1, on='class')
# he1['state'] = he1['class'].str.split('-').str[1].str[0]
# he1 = he1.set_index('class')
# #
# x=he1.drop(columns=['state'])
# zero_count = (x == 0).sum(axis=0)
# y = he1['state']
# y_class=pd.factorize(y)[0]

#kmer2_kmer4
# cluster_kmer=pd.read_csv('tcr_cluster1_kmer_2.csv')
# cluster_kmer2 = pd.read_csv('tcr_cluster1_kmer_4.csv')
# cluster_kmer=cluster_kmer.set_index('sample')
# cluster_kmer_1=cluster_kmer.T
# cluster_kmer_1=cluster_kmer_1.div(cluster_kmer_1.sum(axis=1),axis=0)
# cluster_kmer_1['class']=cluster_kmer_1.index
# cluster_kmer2= cluster_kmer2.set_index('sample')
# cluster_kmer2_1 = cluster_kmer2.T
# cluster_kmer2_1 = cluster_kmer2_1.div(cluster_kmer2_1.sum(axis=1), axis=0)
# cluster_kmer2_1['class'] = cluster_kmer2_1.index
#
# he1 = pd.merge(cluster_kmer_1, cluster_kmer2_1, on='class')
# he1['state'] = he1['class'].str.split('-').str[1].str[0]
# he1 = he1.set_index('class')
#
# x=he1.drop(columns=['state'])
# zero_count = (x == 0).sum(axis=0)
# y = he1['state']
# y_class=pd.factorize(y)[0]

###kmer3_kmer4
# cluster_kmer=pd.read_csv('tcr_cluster1_kmer_3.csv')
# cluster_kmer2 = pd.read_csv('tcr_cluster1_kmer_4.csv')
# cluster_kmer=cluster_kmer.set_index('sample')
# cluster_kmer_1=cluster_kmer.T
# cluster_kmer_1=cluster_kmer_1.div(cluster_kmer_1.sum(axis=1),axis=0)
# cluster_kmer_1['class']=cluster_kmer_1.index
# cluster_kmer2= cluster_kmer2.set_index('sample')
# cluster_kmer2_1 = cluster_kmer2.T
# cluster_kmer2_1 = cluster_kmer2_1.div(cluster_kmer2_1.sum(axis=1), axis=0)
# cluster_kmer2_1['class'] = cluster_kmer2_1.index
#
# he1 = pd.merge(cluster_kmer_1, cluster_kmer2_1, on='class')
# he1['state'] = he1['class'].str.split('-').str[1].str[0]
# he1 = he1.set_index('class')
#
# x=he1.drop(columns=['state'])
# zero_count = (x == 0).sum(axis=0)
# y = he1['state']
# y_class=pd.factorize(y)[0]



#kmer2_kmer3_kmer4
# cluster_kmer=pd.read_csv('tcr_cluster1_kmer_2.csv')
# cluster_kmer2 = pd.read_csv('tcr_cluster1_kmer_3.csv')
# cluster_kmer3 = pd.read_csv('tcr_cluster1_kmer_4.csv')
#
# cluster_kmer=cluster_kmer.set_index('sample')
# cluster_kmer_1=cluster_kmer.T
# cluster_kmer_1=cluster_kmer_1.div(cluster_kmer_1.sum(axis=1),axis=0)
# cluster_kmer_1['class']=cluster_kmer_1.index
# cluster_kmer2= cluster_kmer2.set_index('sample')
# cluster_kmer2_1 = cluster_kmer2.T
# cluster_kmer2_1 = cluster_kmer2_1.div(cluster_kmer2_1.sum(axis=1), axis=0)
# cluster_kmer2_1['class'] = cluster_kmer2_1.index
# cluster_kmer3 = cluster_kmer3.set_index('sample')
# cluster_kmer3_1 = cluster_kmer3.T
# cluster_kmer3_1 = cluster_kmer3_1.div(cluster_kmer3_1.sum(axis=1), axis=0)
# cluster_kmer3_1['class'] = cluster_kmer3_1.index
#
# he1 = pd.merge(cluster_kmer_1, cluster_kmer2_1, on='class')
# he2 = pd.merge(he1, cluster_kmer3_1, on='class')
# he2['state'] = he2['class'].str.split('-').str[1].str[0]
# he2 = he2.set_index('class')
# x=he2.drop(columns=['state'])
# zero_count=(x==0).sum(axis=0)
# y = he2['state']
# y_class=pd.factorize(y)[0]

######
######
######
models = {
    "logisticRegression": LogisticRegression(max_iter=200),
    "svm": SVC(probability=True),
    'mlp': MLPClassifier(),
    'Knn': KNeighborsClassifier(),
    'Decision': DecisionTreeClassifier(),
    'RandomForest': RandomForestClassifier(),
    'XGBoost': XGBClassifier(eval_metric='logloss')
}
# 初始化用于存储最佳结果的字典
best_results = {}

for threshold in range(139, 140, 1):  # 修改为101，以包括100
    columns_keep = zero_count[zero_count <= threshold].index
    x_filter = x[columns_keep]
    feature_names = x_filter.columns
    f_scores, p_values = f_classif(x_filter, y_class)
    # 创建 DataFrame 来存储 ANOVA 分数和 p 值
    anova_results = pd.DataFrame({
        'Feature': feature_names,
        'F-Score': f_scores,
        'P-Value': p_values
    })
    for name, model in models.items():
        model_values = []
        auc_score_list = []
        for i in np.arange(0.001, 0.05, 0.001):
            selected_features = anova_results[anova_results['P-Value'] < i]['Feature']
            x_selected = x_filter[selected_features]
            x_train, x_test, y_train, y_test = train_test_split(x_selected, y_class, test_size=0.2, random_state=42)
            model.fit(x_train, y_train)
            y_score = model.predict_proba(x_test)

            auc_score = []
            for j in range(3):  # 针对每个类
                auc = roc_auc_score(y_test == j, y_score[:, j])
                auc_score.append(auc)

            value = sum(auc_score) / 3
            model_values.append(value)
            auc_score_list.append(auc_score)

        max_value = max(model_values)
        max_index = model_values.index(max_value)
        max_auc_score = auc_score_list[max_index]
        # 存储最佳结果
        best_results[(name, threshold)] = {
            'Max AUC': max_value,
            'Optimal P-Value': (max_index + 1) / 1000,
            'AUC Scores': max_auc_score
        }
# 输出最佳结果
for (model_name, threshold), result in best_results.items():
    print(
        f'Model: {model_name}, Threshold: {threshold}, '
        f'Max AUC: {result["Max AUC"]:.4f}, '
        f'Optimal P-Value: {result["Optimal P-Value"]:.3f}, '
        f'AUC Scores: {result["AUC Scores"]}')