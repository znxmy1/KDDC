#KDDC
KDDC：a new framework that integrates K-mers, Dataset filtering, Dimension reduction and Classification algorithms to achieve immune cell heterogeneity classification

###Data source

step1:immune repertoire data is obtained from GSE158055 of GEO Database.

step2:single cell RNA sequencing data is sourced from the website http://covid19.cancer-pku.cn/#/dimensional-reduction, which can provide a complete h5 file of COVID19.

###Description

bcr/tcr_yaqun.py: includes data input, data cleaning, cell clustering, cell annotation, and visualization sections.

bcr/tcr_kmer_analysis.py: extract patient sample information and cut and integrate cdr3aa sequences according to different lengths.

bcr/tcr_cluster1_classification_methods.py: reduce the dimensionality of specific cell subpopulations before classification, compare the performance of different classification algorithms, and evaluate them.
bcr_metabolism.py：alculate pathway scores for cell subpopulations in B cells.
bcr_protein.R: calculate the intercellular surface protein score of cell subpopulations in B cells.
In addition, we also provided some intermediate process data to verify the completeness of our analysis.

###contact
Feel free to submit an issue or contact us at 2264311577@qq.com for problems about the tool.
