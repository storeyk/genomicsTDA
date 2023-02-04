## file: tsne_t0
## author: Farzana_Nasrin
## date: 8/28/2022
## This code produces the clusters using t-distributed stochastic neighbor embedding (tsne) 
## of the T0 scores obtained by standardizing Gaussian mixture feeting 



library(tidyverse) ## ggplot2 is a package of tidyverse
library(Rtsne)  
cts <-  read.csv('Lung_T0.csv')
ctsmatrix = as.matrix(cts[,-1])
rownames(ctsmatrix) = cts[,1]

n_healthy = 314
n_tumor = 500
transpose_cts1 = as.matrix(t(ctsmatrix))
transpose_cts2 = as.matrix(c(rep("healthy", n_healthy), rep("tumor",n_tumor)))
transpose_cts = cbind(transpose_cts1,transpose_cts2)  



set.seed(142)
tSNE_fit <- Rtsne(transpose_cts1)
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2")

tSNE_df <- cbind(tSNE_df,transpose_cts2)

tSNE_df %>% head()

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = transpose_cts2))+
  geom_point(alpha = 1,size=2)+
  theme_classic(16)+
  theme(axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)),
        axis.title.y = element_text(size=rel(1.75)),
        axis.title.x = element_text(size=rel(1.75)))+
  scale_color_manual(breaks = c("healthy", "tumor"),
                     values=c("purple", "yellow"))