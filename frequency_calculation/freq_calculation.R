library(tidyverse)
library(ggpubr)
library(data.table)
library(Hmisc)
library(glue)
library(binom)


df <- rbind(comp_others4, modern_df)
summary(df)

# Convert read count columns to numeric format
# T = total number of reads
# R = number of reads carrying the risk allele
df$T <- as.numeric(df$T)
df$R <- as.numeric(df$R)


# Split the dataset by Group, SNP_ID, and Phenotype
# Each subset represents one population/SNP/phenotype combination
df_subset <- split(df, list(df$Group, df$SNP_ID, df$Phenotype))

# Remove empty subsets
df_subset_v1 <- df_subset[sapply(df_subset, function(x) dim(x)[1]) > 0]

# Rename subsets with simple numeric names
df_subset_renamed <- setNames(
  df_subset_v1,
  as.vector(1:length(df_subset_v1))
)


# ============================================================
# Allele frequency estimation using maximum likelihood
# ============================================================

# r   = number of reads carrying the risk allele
# t   = total number of reads
# p   = population allele frequency to be estimated
# eps = sequencing error rate
#
# The model assumes three possible genotypes:
# RR: homozygous risk genotype
# Rr: heterozygous genotype
# rr: homozygous non-risk genotype
#
# For each individual, the likelihood is computed as:
#
# P(data | p) =
#   p^2        * P(reads | RR) +
#   2p(1 - p) * P(reads | Rr) +
#   (1 - p)^2 * P(reads | rr)
#
# Since dbinom(..., log = TRUE) returns log-probabilities,
# the likelihood is calculated in log-space using the log-sum-exp trick.

freq_loop <- function(dataframe) {
  
  # Compute the per-individual log-likelihood contribution
  function_for_freq <- function(r, t, p, eps = 0.01) {
    
    log_terms <- c(
      log(p^2)         + dbinom(r, t, 1 - eps, log = TRUE),  # RR genotype
      log(2*p*(1 - p)) + dbinom(r, t, 0.5,     log = TRUE),  # Rr genotype
      log((1 - p)^2)   + dbinom(r, t, eps,     log = TRUE)   # rr genotype
    )
    
    # Log-sum-exp trick for numerical stability.
    # This prevents underflow when likelihood values are very small.
    m <- max(log_terms)
    m + log(sum(exp(log_terms - m)))
  }
  
  
  # Candidate allele frequencies evaluated on a fixed grid.
  # A smaller step size gives more precise estimates but takes longer.
  ps <- seq(0, 1, by = 0.01)
  
  # Store log-likelihood values for each candidate allele frequency
  logLike <- matrix(NA, length(ps), 1)
  rownames(logLike) <- ps
  
  
  # For each candidate allele frequency, sum the log-likelihood
  # contributions across all individuals in the current subgroup
  for (p in ps) {
    
    total_log_likelihood <- 0
    
    for (i in 1:nrow(dataframe)) {
      total_log_likelihood <- total_log_likelihood + function_for_freq(
        r = dataframe$R[i],
        t = dataframe$T[i],
        p = p,
        eps = 0.01
      )
    }
    
    logLike[rownames(logLike) == p, ] <- total_log_likelihood
  }
  
  
  # Extract subgroup metadata
  period_name <- as.character(unique(dataframe$Group))
  snp_name <- as.character(unique(dataframe$SNP_ID))
  phenotype_name <- as.character(unique(dataframe$Phenotype))
  pop_size <- as.numeric(nrow(dataframe))
  
  
  # Difference from the maximum log-likelihood.
  # This is useful for comparing likelihood support across frequencies.
  dlogLike <- logLike - max(logLike)
  
  
  # Maximum likelihood estimate of allele frequency
  pHat <- ps[logLike == max(logLike)]
  pHat_max <- as.numeric(pHat)
  
  
  # Wilson score confidence interval
  #
  # Here, pHat is converted into an approximate count within the subgroup.
  # This provides a binomial confidence interval around the estimated frequency.
  x <- pHat_max * pop_size / pop_size
  
  CI <- binom.confint(
    x * pop_size,
    pop_size,
    method = "wilson",
    type = "central"
  )
  
  
  # Store likelihood results and subgroup-level metadata
  likeResults <- data.frame(
    frequency = ps,
    diffLogLike = dlogLike,
    Loglikelihood = logLike,
    Period = period_name,
    SNPid = snp_name,
    Phenotype = phenotype_name,
    POPsize = pop_size,
    pHat = pHat_max,
    CI_Lower = CI[, 5],
    CI_Upper = CI[, 6]
  )
  
  likeResults$Loglikelihood <- as.numeric(likeResults$Loglikelihood)
  
  
  # Return only the row corresponding to the maximum likelihood estimate
  maxlikelihood_result <- likeResults[which.max(likeResults$Loglikelihood), ]
  
  return(maxlikelihood_result)
}


# ============================================================
# Run frequency estimation for all subgroups
# ============================================================

# Empty data frame to collect maximum likelihood results
df_empty <- data.frame()

# Apply the frequency estimation function to each subgroup
for (i in 1:length(df_subset_renamed)) {
  
  output <- freq_loop(df_subset_renamed[[i]])
  
  df_empty <- rbind(df_empty, output)
  
  maxll_results <- df_empty
}


# Final frequency result table
comp_others_freq <- maxll_results
all_freq <- maxll_results


# Export results as CSV
write.csv(
  all_freq,
  "/Users/dilanur/Desktop/all_freq.csv",
  row.names = FALSE
)
