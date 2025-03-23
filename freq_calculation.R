library(tidyverse)
library(ggpubr)
library(data.table)
library(dplyr)
library(Hmisc)
library(glue)
library(binom)


df = rbind(comp_others4, modern_df)
summary(df)

df$T <- as.numeric(df$T) #total reads
df$R <- as.numeric(df$R) #risk allele containing reads


df_subset = split(df, list(df$Group, df$SNP_ID, df$Phenotype))
df_subset_v1 = df_subset[sapply(df_subset, function(x) dim(x)[1]) > 0]
df_subset_renamed = setNames(df_subset_v1, as.vector(1:length(df_subset_v1)))


### compute frequency

#r=number of risk allele copies
#t=number of total allele copies
#p=allele frequency -> we would like to find this
#error rate epsilon -> 0.01 

freq_loop <-  function(dataframe) {
  function_for_freq <- function( r  ,t ,p) { 
    (((p^2*dbinom(r,t,1-0.01,log=T)) + ( 2*p*(1-p)*dbinom(r,t,0.5,log=T) ) + ( (1-p)^2*dbinom(r,t,0.01,log=T))))
  } 
  
  ps=seq(0,1,by=0.01)
  logLike=matrix(NA,length(ps),1)
  rownames(logLike)=ps
  for (p in ps) { 
    x=0
    for (i in 1:nrow(dataframe)) {
      x=x+(function_for_freq(dataframe$R[i],dataframe$T[i],p))
      logLike[rownames(logLike)==p,]=x} 
  }
  period_name <- as.character(unique(dataframe$Group))
  snp_name <- as.character(unique(dataframe$SNP_ID))
  phenotype_name <- as.character(unique(dataframe$Phenotype))
  pop_size <- as.numeric(nrow(dataframe))
  dlogLike <- logLike - max(logLike)
  pHat <- ps[logLike == max(logLike)]     #Maximum likelihood estimate
  pHat_max <- as.numeric(pHat)
  
  # Wilson score interval
  alpha <- 0.05
  x <- pHat_max * pop_size / pop_size
  CI <- binom.confint(x * pop_size, pop_size, method = "wilson", type = "central")

  
  #put the result into a data frame
  likeResults <- data.frame(frequency = ps, diffLogLike = dlogLike, Loglikelihood = logLike,  Period= period_name, 
                            SNPid= snp_name, Phenotype=phenotype_name, POPsize= pop_size ,  pHat= pHat_max, CI_Lower= CI[,5] ,   CI_Upper = CI[,6])
  likeResults$Loglikelihood <- as.numeric(likeResults$Loglikelihood)
  maxlikelihood_result <- likeResults[which.max(likeResults$Loglikelihood),]
}


df_empty = data.frame()

for (i in 1:length(df_subset_renamed)) { 
  output = freq_loop(df_subset_renamed[[i]])
  df_empty = rbind(df_empty, output)
  maxll_results <- df_empty
}





comp_others_freq = maxll_results
all_freq = maxll_results
write.csv(all_freq, "/Users/dilanur/Desktop/all_freq.csv")

