library(tidyverse)
library(normtools)
library(randomForest)
library(corrplot)
library(doParallel)

# Load data
dat <- read.csv("/hpf/largeprojects/MICe/yyee/dev/normtools/normtools/data/mouse_brain_ABI_coronal_gene_expression.csv", row.names = 1, check.names = F)

# Setup params
N_cores <- 40
#N_genes <- 12
N_genes <- 4345

# Sampling
#rand_genes <- sample(1:dim(dat)[1], size=N_genes, replace=F)
rand_genes <- 1:N_genes
rand_gene_names <- rownames(dat)[rand_genes]

# Set up output objects
dat_rf <- normtools:::construct_like_matrix(dat)[rand_genes,]
dat_rf_importance <- matrix(nrow=N_genes, ncol=dim(dat)[1], dimnames = list(rand_gene_names, rownames(dat)))
dat_rf_rsq <- numeric(N_genes) %>% `names<-`(rand_gene_names)

# Parallelization
cl <- makeCluster(N_cores, outfile="", rscript_args="--vanilla")
registerDoParallel(cl)
prog <- txtProgressBar(min=1, max=N_genes, style=3)
export_functions <- c("randomForest")
export_variables <- c("dat", "rand_genes")

# Loop
output <- foreach(i=1:N_genes, .export = c(export_functions), .verbose = FALSE) %dopar% {

  x <- t(dat)[,-c(rand_genes[i])]
  y <- t(dat)[,rand_genes[i]]
  rf <- randomForest(x = x, y=y)

  normalized_value <- (rf$y - rf$predicted)
  importance <- append(rf$importance, 0, after = rand_genes[i]-1)
  rsq <- rf$rsq[500]

  out <- list(normalized=normalized_value, importance=importance, rsq=rsq, gene=rand_gene_names[i])
  setTxtProgressBar(prog, i)
  return(out)
  #dat_rf[i, ] <- (rf$y - rf$predicted)
  #dat_rf_importance[i,] <- append(rf$importance, 0, after = rand_genes[i]-1)
  #dat_rf_rsq[i] <- rf$rsq[500]
}
close(prog)
registerDoSEQ()

save(list=c("output"), file="/hpf/largeprojects/MICe/yyee/dev/normtools/normtools/data/rf_cleaned.RData")
load("/hpf/largeprojects/MICe/yyee/dev/normtools/normtools/data/rf_cleaned.RData")
