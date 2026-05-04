freq_loop <- function(dataframe) {
  
  # Compute the per-individual log-likelihood contribution
  # r   = number of reads carrying the risk allele
  # t   = total number of reads
  # p   = population allele frequency to be estimated
  # eps = sequencing error rate
  function_for_freq <- function(r, t, p, eps = 0.01) {
    
    # Genotype likelihood model:
    # RR: homozygous risk genotype
    # Rr: heterozygous genotype
    # rr: homozygous non-risk genotype
    #
    # The likelihood is:
    # P(data | p) =
    #   p^2        * P(reads | RR) +
    #   2p(1 - p) * P(reads | Rr) +
    #   (1 - p)^2 * P(reads | rr)
    #
    # Since dbinom(..., log = TRUE) returns log-probabilities,
    # all terms are computed in log-space.
    
    log_terms <- c(
      log(p^2)         + dbinom(r, t, 1 - eps, log = TRUE),  # RR genotype
      log(2*p*(1 - p)) + dbinom(r, t, 0.5,     log = TRUE),  # Rr genotype
      log((1 - p)^2)   + dbinom(r, t, eps,     log = TRUE)   # rr genotype
    )
    
    # Use the log-sum-exp trick for numerical stability.
    # This avoids underflow when likelihood values are very small.
    m <- max(log_terms)
    m + log(sum(exp(log_terms - m)))
  }
  
  # Candidate allele frequencies evaluated on a fixed grid
  ps <- seq(0, 1, by = 0.01)
  
  # Store log-likelihood values for each candidate frequency
  logLike <- matrix(NA, length(ps), 1)
  rownames(logLike) <- ps
  
  # For each candidate frequency, sum the log-likelihood contribution
  # across all individuals in the current subgroup
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
  
  # Difference from the maximum log-likelihood
  # Useful for comparing likelihood support across frequencies
  dlogLike <- logLike - max(logLike)
  
  # Maximum likelihood estimate of allele frequency
  pHat <- ps[logLike == max(logLike)]
  pHat_max <- as.numeric(pHat)
  
  # Wilson score confidence interval
  #
  # Here, pHat is converted into an approximate count of risk alleles
  # within the subgroup. This gives a binomial confidence interval
  # around the estimated allele frequency.
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
  
  # Return only the row with the maximum likelihood estimate
  maxlikelihood_result <- likeResults[which.max(likeResults$Loglikelihood), ]
  
  return(maxlikelihood_result)
}
