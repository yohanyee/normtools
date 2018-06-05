library(tidyverse)
library(normtools)
library(randomForest)
library(corrplot)

# Load data
dat <- read.csv("data/mouse_brain_ABI_coronal_gene_expression.csv", row.names = 1, check.names = F)

dat_rel <- norm_row_sum(dat)
dat_reg <- norm_row_regression(dat)
dat_nbm <- norm_by_model(dat, ~ total, df=data.frame(total=rowSums(dat)))
dat_suc <- norm_successive(dat)


o <- dat_rel %>% cor %>% dist %>% hclust(method="complete") %>% .$order



N_genes <- 10
rand_genes <- sample(1:dim(dat)[1], size=N_genes, replace=F)
rand_gene_names <- rownames(dat)[rand_genes]




dat_rf <- normtools:::construct_like_matrix(dat)[rand_genes,]
dat_rf_importance <- matrix(nrow=N_genes, ncol=dim(dat)[1], dimnames = list(rand_gene_names, rownames(dat)))
dat_rf_rsq <- numeric(N_genes) %>% `names<-`(rand_gene_names)

for (i in 1:N_genes) {
  print(i)
  x <- t(dat)[,-c(rand_genes[i])]
  y <- t(dat)[,rand_genes[i]]
  rf <- randomForest(x = x, y=y)
  dat_rf[i, ] <- (rf$y - rf$predicted)
  dat_rf_importance[i,] <- append(rf$importance, 0, after = rand_genes[i]-1)
  dat_rf_rsq[i] <- rf$rsq[500]
}

corrplot(cor(dat_rel)[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_reg)[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_nbm)[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_suc)[o,o], method="color", tl.col="black", tl.cex=0.08)

qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_nbm)[o,o]))
qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_suc)[o,o]))
qplot(as.numeric(cor(dat_nbm)[o,o]), as.numeric(cor(dat_suc)[o,o]))

corrplot(cor(dat_rel)[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_rel[rand_genes,])[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_reg[rand_genes,])[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_suc[rand_genes,])[o,o], method="color", tl.col="black", tl.cex=0.08)
corrplot(cor(dat_rf)[o,o], method="color", tl.col="black", tl.cex=0.08)


qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_rel[rand_genes,])[o,o]))
qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_reg[rand_genes,])[o,o]))
qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_suc[rand_genes,])[o,o]))
qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_rf)[o,o]))

qplot(as.numeric(cor(dat_rel)[o,o]), as.numeric(cor(dat_rel[rand_genes,])[o,o]))
qplot(as.numeric(cor(dat_reg)[o,o]), as.numeric(cor(dat_reg[rand_genes,])[o,o]))

dat_rf_rsq
