# genomicsTDA

Code:

**compute_t_score_GMM.ipynb** - python jupyter notebook to compute the T scores for a Gaussian Mixture Model describing RNA-seq FPKM data.

**main_mapper_graph.ipynb** - python jupyter notebook to construct a Mapper graph using T0 scores for GMM. This code produces an image of the mapper graph and saves matrices storing node and edge information, to be used in the Heat Kernel Signature.

**mapper_PI.ipynb** - python jupyter notebook to compute the position index (PI) for each subject in a Mapper graph. This code allows the user to loop through a range of parameter values and sum the PI for each subject over all mapper graphs.

**mapper_tumor_utils.py** - python code containing functions used in the three jupyter notebooks above. 

Data:

**Lung_FPKM.csv** - FPKM (Fragments Per Kilobase Million) expression data from 19648 distinct genes, for 314 healthy subjects (from GTEx) and 500 lung tumor subjects (from TCGA). This is used in **compute_t_score_GMM.ipynb**.

Lung_T0.csv - T0 scores for the FPKM expression data above, to be used in **main_mapper_graph.ipynb** and **mapper_PI.ipynb**.
