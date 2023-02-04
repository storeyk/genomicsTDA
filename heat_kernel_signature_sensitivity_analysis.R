## file: heat_kernel_signature_sensitivity_analysis
## author: Farzana_Nasrin
## date: 8/1/2022
## This code produces the confidence intervals from the empirical distribution 
## of the GSSs of each subject in the data set. 
## Input: a matrix (W) for the weight of edges, a diagonal matrix (M) for weights of vertices, and a matrix (V) for the vertex set
## In particular W is generated from mapper graphs where the diagonal entries are the number of subject at each vertex and non-diagonal 
## entries (i,j) are the number of overlapping subjects between vertex i and j.
## the "mat.vertices.data" encodes the proportion of non-healthy subjects present at each vertex
## Output: a vector of GSSs and CIs of those estimations


library(matrixcalc)
library(asymmetry)
require(pracma)
require(gmodels)
library(Rmisc)



mat.edge.data = read.table(paste("Lung_T0_meancorr_eps600_r60_g80_edge_mat.csv",sep = ","),sep = ",",row.names = NULL)
mat.vertices.data = read.table(paste("lung_T0_meancorr_eps600_r60_g80_prop_mat.csv",sep = ","),sep = ",",row.names = NULL)
mat.subject.data = read.table(paste("Lung_T0_meancorr_eps600_r60_g80_subj_node_mat.csv",sep = ","),sep = ",",row.names = NULL)

matrix.M.0 = as.matrix(mat.vertices.data)
matrix.M.0.diag = diag(matrix.M.0)
matrix.W.0 = as.matrix(mat.edge.data)
sum_D.0 = colSums(matrix.W.0)
matrix.D.0 = diag(sum_D.0, dim(mat.edge.data)[1], dim(mat.edge.data)[1])
matrix.S.0 = as.matrix(mat.subject.data)
matrix.S.0 = matrix.S.0[,colSums(matrix.S.0) > 0]
subject.number = dim(matrix.S.0)[1]

vertices_order.0 = matrix.M.0.diag[order(-matrix.M.0.diag)]


idx =  ! diag(matrix.M.0) ==0
matrix.M = matrix.M.0[ idx , idx ]
node.number = dim(matrix.M)[1]
matrix.W= matrix.W.0[ idx , idx ]
matrix.S = matrix.S.0[  , idx ]
matrix.M.diag = diag(matrix.M)
matrix.D= matrix.D.0[ idx , idx ]
matrix.M.inv = diag(1/sqrt(matrix.M.diag), node.number,node.number)
weighted_Laplacian = matrix.M.inv %*% (matrix.D-matrix.W) %*% matrix.M.inv

asymm_index = function(m) {
  which(m != t(m), arr.ind = TRUE)
}

asymm_index(weighted_Laplacian)

eigen_list = eigen(weighted_Laplacian)
eigen_values = as.matrix(eigen_list$values, ncol =1)
eigen_vectors = eigen_list$vectors

t.number = 20
t = seq(from = 0, to = 1,length.out = t.number)
s_matrix =matrix(0, nrow = t.number, ncol = node.number)
## heat kernel signatures

for (i in (1:t.number)){
  for (j  in (1:node.number)){
    to_h = apply(eigen_values, 1, function(x){exp(-t[i]*x)})  
    phi = eigen_vectors[j,]
    s_matrix[i,j] =  dot(to_h,phi)
  }
}
  
matplot((s_matrix), type = "l")
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
HKS = rep.row(s_matrix[t.number,],node.number)

GSS = matrix.S %*% apply(t(s_matrix[t.number,]),2,c)

GSS_df = rbind(t(GSS_600_60_30), t(GSS_600_60_40), t(GSS_600_60_50), t(GSS_600_60_60), t(GSS_600_60_70), t(GSS_600_60_80))
GSS_mat = matrix(GSS_df, ncol = subject.number, byrow = TRUE)
GSS_CI = sapply(as.data.frame(GSS_df),function(x){CI(x, ci= 0.98)})
boxplot(as.data.frame(GSS_df))
matplot(t(GSS_CI), type = "l")
