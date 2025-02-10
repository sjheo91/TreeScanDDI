# Simulation Example
rm(list=ls())

source('code/function.R')
library(pROC)
library(MASS)
library(openxlsx)
library(doSNOW)
library(parallel)

numCores<-detectCores()
cl <- makeCluster(numCores)
registerDoSNOW(cl)

TS <- read.csv('data/Tree_Structure.csv')
set.seed(1)
tree <- sample(TS$TS,100,replace=F)
tree_str <- tree_split(tree)

## Scenario - No signal
n_iter <- 1000
n_signal <- c(1, 3, 5)
H1 <- c('1/1/1/1','2/2/1/1','3/2/2/1')
n_size <- c(10, 5, 1)
AE_bias <- c(0, 0.05, 0.1, 0.2)
drug_bias <- c(0, 0.05, 0.1, 0.2)
cut <- c(0)

grid <- expand.grid(n_signal=n_signal, H1=H1, n_size=n_size, AE_bias=AE_bias, drug_bias=drug_bias, 
                    cut=cut, stringsAsFactors=F)
res_list <- list()
grid <- grid[grid$AE_bias==0.05|grid$drug_bias==0.05,]

for(j in 1:nrow(grid)){
  
  pb <- txtProgressBar(max=n_iter, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  start <- proc.time()
  res <- foreach(i=1:n_iter, .combine=comb, .multicombine=T, .packages='pROC', .options.snow=opts) %dopar% {
    run_simulation(iter=i, tree_str=tree_str, n_AE=100, H0=c(1,1,1,1), num_AE_bias=100, 
                   H1=as.numeric(unlist(strsplit(grid$H1[j],'/',))), 
                   margin_drug=c(20000, 200000, 200000, 5000000)/grid$n_size[j],
                   drug_bias=c(grid$drug_bias[j]^2,grid$drug_bias[j],grid$drug_bias[j],0), 
                   n_signal=grid$n_signal[j],
                   AE_bias=grid$AE_bias[j],
                   cut=grid$cut[j])
  }
  proc.time()-start
  close(pb)
  print(j)
  
  tmp <- do.call(rbind,lapply(res, summary_perf))
  res_list[[j]] <- cbind(method=row.names(tmp), 
                         sapply(grid[j,],rep,nrow(tmp)), 
                         type_I_error=tmp[,1])
}

res_data <- data.frame(do.call(rbind,res_list))
row.names(res_data) <- 1:nrow(res_data)
write.xlsx(res_data, 'type_I_error.xlsx', row.names=F, overwrite=T)
