# genomicsTDA

Code:

**compute_t_score_GMM.ipynb** - python jupyter notebook to compute the T scores for a Gaussian Mixture Model describing RNA-seq FPKM data.

**deseq_two_groups.R** - R code to conduct the differential gene expression analysis on the raw FPKM data, producing the log2 fold changes, p-values, and adjusted p-values.

**heat_kernel_signature_sensitivity_analysis.R** - R code to produce the confidence intervals from the empirical distribution of the GSSs of each subject in the data set.

**main_mapper_graph.ipynb** - python jupyter notebook to construct a Mapper graph using T0 scores for GMM. This code produces an image of the mapper graph and saves matrices storing node and edge information, to be used in the Heat Kernel Signature.

**mapper_PI.ipynb** - python jupyter notebook to compute the position index (PI) for each subject in a Mapper graph. This code allows the user to loop through a range of parameter values and sum the PI for each subject over all mapper graphs.

**mapper_tumor_utils.py** - python code containing functions used in the three jupyter notebooks above. 

**tsne_t0.R** - R code to produce the clusters using t-distributed stochastic neighbor embedding (tsne) of the T0 scores for the GMM.

Data:

**data_molecular.zip** - a folder containing the molecular expression data, to be used in **deseq_two_groups.R**, **tsne_t0.R**, **main_mapper_graph.ipynb**, and **mapper_PI.ipynb**.

**hks_matrices.zip** - a folder containing matrices describing Mapper graphs, to be used in **heat_kernel_signature_sensitivity_analysis.R**.

**Lung_FPKM.csv** - FPKM (Fragments Per Kilobase Million) expression data from 19648 distinct genes, for 314 healthy subjects (from GTEx) and 500 lung tumor subjects (from TCGA). This is used in **compute_t_score_GMM.ipynb**.

