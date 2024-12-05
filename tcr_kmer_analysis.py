# -*- coding: utf-8 -*-

import pandas as pd #数据分析
import itertools #迭代器
import pandas as pd #数据分析
import os#读取目录与多个文件操作

####数据集拆分
# data=pd.read_csv('tcr_he.csv')
# print(data.head())
# print(data.shape)
# print(data['sampleID'].value_counts())
# grouped=data.groupby('sampleID')
# subsets = {sampleID: group for sampleID, group in grouped}
# for sampleID, subset in subsets.items():
#     subset.to_csv(f"{sampleID}_subset.csv", index=False)



amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
# 生成所有可能的k-mers
def generate_kmers(k):
    return [''.join(p) for p in itertools.product(amino_acids, repeat=k)]
#itertools.product产生多个列表和迭代器的积

# 从序列中提取k-mers
def extract_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmers.append(sequence[i:i + k])
    return kmers

# 生成k-mer矩阵
def kmer_matrix(sequences, k, sequence_freq):
    kmers = generate_kmers(k)
    matrix = []

    for seq, freq in zip(sequences, sequence_freq):
        kmer_counts = {kmer: 0 for kmer in kmers}
        seq_kmers = extract_kmers(seq, k)
        for kmer in seq_kmers:
            kmer_counts[kmer] += 1
        matrix.append([kmer_counts[kmer] * freq/len(seq) for kmer in kmers])

    return pd.DataFrame(matrix, columns=kmers)

folder_path='D:/immune repertoire/TCR/tcr_cell/cluster5'
result_df = pd.DataFrame()
# 读取每个样本文件并处理
for filename in os.listdir(folder_path):
    # 检查文件是否为 CSV 文件
    if filename.endswith(".csv"):
        # 构建文件的完整路径
        file_path = os.path.join(folder_path, filename)
        # 读取 CSV 文件并将第一行作为列名
        sample_data = pd.read_csv(file_path)
        # 提取样本的氨基酸序列和频率
        sequences = sample_data.iloc[:,11]
        length = sample_data['cdr3aa_length']
        # 设置k值
        k=4
        # 生成k-mer矩阵
        matrix = kmer_matrix(sequences, k, length)
        #print(matrix)
        column_sums = matrix.sum(axis=0)#axis=0沿着行的方向求和，会计算每一列的和
        col_name = filename.split('.', 1)[0]
        column_sums_df = pd.DataFrame(column_sums, columns=['{}'.format(col_name)])
        # 将列求和添加到总和矩阵中
        result_df = pd.concat([result_df, column_sums_df], axis=1)
result_df.to_csv('tcr_cluster5_kmer_4.csv',index=True)

