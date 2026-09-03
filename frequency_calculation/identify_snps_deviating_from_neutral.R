######################## Store SNPs showing significant deviations from neutral expectations ########################

all_data_fig = all_data_df


# 1) Initial dataset
ds <- all_data_fig
head(ds)
tail(ds)
dim(ds)


### Neolithic (P1) – Late Chalcolithic (P2) ANALYSIS ###
### Confidence intervals taken from Late Chalcolithic (P2) ###

# 2) Retain only SNPs with data available in both P1 and P2
#    P3 may or may not be present; only P1 and P2 are required for this comparison
ds2 <- ds %>%
  select(frequency, Period, SNPid, Phenotype, POPsize, CI_Lower, CI_Upper) %>%
  group_by(SNPid) %>%
  filter(all(c("P1", "P2") %in% Period)) %>%  # both P1 and P2 are required
  ungroup()

ds2
dim(ds2)


# 3) Calculate Neolithic - Late Chalcolithic frequency differences
ds3 <- ds2 %>%
  group_by(SNPid) %>%
  mutate(
    Neol_minus_LateChal = frequency[Period == "P1"] - frequency[Period == "P2"]
  ) %>%
  ungroup()


# Use ds4 in the SNP-based loop
ds4 <- ds3
head(ds4)


# 4) SNP list
SNPids <- unique(ds4$SNPid)


#### Neolithic - Late Chalcolithic comparison using Late Chalcolithic confidence intervals ####

# Results data frame
results2_n_lc_filtered <- data.frame(
  SNPid     = character(),
  Status    = character(),
  Phenotype = character(),
  stringsAsFactors = FALSE
)


# 5) SNP-based loop
for (i in seq_along(SNPids)) {

  current_snp <- SNPids[[i]]

  # Obtain the confidence interval for the current SNP from P2 (Late Chalcolithic)
  lowerCI <- ds4$CI_Lower[ds4$SNPid == current_snp & ds4$Period == "P2"][1]
  upperCI <- ds4$CI_Upper[ds4$SNPid == current_snp & ds4$Period == "P2"][1]

  # Select neutral SNPs whose P2 frequencies fall within the
  # confidence interval of the current SNP
  ds_rs <- ds4 %>%
    select(
      frequency, Period, Neol_minus_LateChal,
      SNPid, Phenotype, POPsize, CI_Lower, CI_Upper
    ) %>%
    group_by(SNPid) %>%
    filter(
      (Phenotype == "Neutral" &
         frequency >= lowerCI &
         frequency <= upperCI &
         Period == "P2") |
        SNPid == current_snp
    ) %>%
    ungroup()

  # Construct the reference distribution from matched neutral SNPs
  neutral_vals <- ds_rs$Neol_minus_LateChal[ds_rs$Phenotype == "Neutral"]

  # Skip the current SNP if no matched neutral SNPs are available
  if (length(neutral_vals) == 0) next

  # Calculate the 2.5th and 97.5th percentiles of the neutral distribution
  p97_5 <- quantile(neutral_vals, 0.975, na.rm = TRUE)
  p2_5  <- quantile(neutral_vals, 0.025, na.rm = TRUE)

  # Frequency difference and phenotype of the current SNP
  snp_value <- ds_rs$Neol_minus_LateChal[ds_rs$SNPid == current_snp][1]
  phenotype <- unique(ds_rs$Phenotype[ds_rs$SNPid == current_snp])[1]

  # Identify SNPs falling outside the neutral distribution
  if (snp_value > p97_5) {
    results2_n_lc_filtered <- rbind(
      results2_n_lc_filtered,
      data.frame(
        SNPid = current_snp,
        Status = "decrease",
        Phenotype = phenotype
      )
    )
  } else if (snp_value < p2_5) {
    results2_n_lc_filtered <- rbind(
      results2_n_lc_filtered,
      data.frame(
        SNPid = current_snp,
        Status = "increase",
        Phenotype = phenotype
      )
    )
  }
}

results2_n_lc_filtered

write.csv(
  results2_n_lc_filtered,
  "/Users/dilanur/Desktop/results_n_lc_filtered_son.csv"
)




### Late Chalcolithic (P2) – Modern (P3) ANALYSIS ###
### Confidence intervals taken from Modern (P3) ###

# 1) Retain SNPs with data available in both P2 and P3
#    P1 may or may not be present
ds2_lc_mod <- ds %>%
  select(frequency, Period, SNPid, Phenotype, POPsize, CI_Lower, CI_Upper) %>%
  group_by(SNPid) %>%
  filter(all(c("P2", "P3") %in% Period)) %>%  # both P2 and P3 are required
  ungroup()

ds2_lc_mod
dim(ds2_lc_mod)


# 2) Calculate Late Chalcolithic - Modern frequency differences
ds3_lc_mod <- ds2_lc_mod %>%
  group_by(SNPid) %>%
  mutate(
    LateChal_minus_Modern = frequency[Period == "P2"] - frequency[Period == "P3"]
  ) %>%
  ungroup()


# Dataset used in the SNP-based loop
ds4_lc_mod <- ds3_lc_mod
head(ds4_lc_mod)


# 3) SNP list
SNPids_lc_mod <- unique(ds4_lc_mod$SNPid)


# Results data frame
results2_lc_m_filtered <- data.frame(
  SNPid     = character(),
  Status    = character(),
  Phenotype = character(),
  stringsAsFactors = FALSE
)


# 4) SNP-based loop
# Late Chalcolithic - Modern comparison using Modern confidence intervals
for (i in seq_along(SNPids_lc_mod)) {

  current_snp <- SNPids_lc_mod[[i]]

  # Obtain the confidence interval for the current SNP from P3 (Modern)
  lowerCI <- ds4_lc_mod$CI_Lower[
    ds4_lc_mod$SNPid == current_snp &
      ds4_lc_mod$Period == "P3"
  ][1]

  upperCI <- ds4_lc_mod$CI_Upper[
    ds4_lc_mod$SNPid == current_snp &
      ds4_lc_mod$Period == "P3"
  ][1]

  # Select neutral SNPs whose P3 frequencies fall within the
  # confidence interval of the current SNP
  ds_rs_lc_mod <- ds4_lc_mod %>%
    select(
      frequency, Period, LateChal_minus_Modern,
      SNPid, Phenotype, POPsize, CI_Lower, CI_Upper
    ) %>%
    group_by(SNPid) %>%
    filter(
      (Phenotype == "Neutral" &
         frequency >= lowerCI &
         frequency <= upperCI &
         Period == "P3") |
        SNPid == current_snp
    ) %>%
    ungroup()

  # Matched neutral reference distribution
  neutral_vals <-
    ds_rs_lc_mod$LateChal_minus_Modern[
      ds_rs_lc_mod$Phenotype == "Neutral"
    ]

  # Skip the current SNP if no matched neutral SNPs are available
  if (length(neutral_vals) == 0) next

  p97_5 <- quantile(neutral_vals, 0.975, na.rm = TRUE)
  p2_5  <- quantile(neutral_vals, 0.025, na.rm = TRUE)

  # Frequency difference and phenotype of the current SNP
  snp_value <-
    ds_rs_lc_mod$LateChal_minus_Modern[
      ds_rs_lc_mod$SNPid == current_snp
    ][1]

  phenotype <-
    unique(
      ds_rs_lc_mod$Phenotype[
        ds_rs_lc_mod$SNPid == current_snp
      ]
    )[1]

  # Identify SNPs falling outside the neutral distribution
  # Positive difference -> decrease toward the Modern period
  # Negative difference -> increase toward the Modern period
  if (snp_value > p97_5) {
    results2_lc_m_filtered <- rbind(
      results2_lc_m_filtered,
      data.frame(
        SNPid = current_snp,
        Status = "decrease",
        Phenotype = phenotype
      )
    )
  } else if (snp_value < p2_5) {
    results2_lc_m_filtered <- rbind(
      results2_lc_m_filtered,
      data.frame(
        SNPid = current_snp,
        Status = "increase",
        Phenotype = phenotype
      )
    )
  }
}

results2_lc_m_filtered

write.csv(
  results2_lc_m_filtered,
  "/Users/dilanur/Desktop/results_lc_m_filtered_son.csv"
)




### Neolithic (P1) – Modern (P3) ANALYSIS ###
### Confidence intervals taken from Modern (P3) ###

# 1) Retain SNPs with data available in both P1 and P3
#    P2 may or may not be present
ds2_n_mod <- ds %>%
  select(frequency, Period, SNPid, Phenotype, POPsize, CI_Lower, CI_Upper) %>%
  group_by(SNPid) %>%
  filter(all(c("P1", "P3") %in% Period)) %>%  # both P1 and P3 are required
  ungroup()

ds2_n_mod
dim(ds2_n_mod)


# 2) Calculate Neolithic - Modern frequency differences
ds3_n_mod <- ds2_n_mod %>%
  group_by(SNPid) %>%
  mutate(
    Neol_minus_Modern = frequency[Period == "P1"] - frequency[Period == "P3"]
  ) %>%
  ungroup()


# Dataset used in the SNP-based loop
ds4_n_mod <- ds3_n_mod
head(ds4_n_mod)


# 3) SNP list
SNPids_n_mod <- unique(ds4_n_mod$SNPid)


# Results data frame
results2_n_m_filtered <- data.frame(
  SNPid     = character(),
  Status    = character(),
  Phenotype = character(),
  stringsAsFactors = FALSE
)


# 4) SNP-based loop
# Neolithic - Modern comparison using Modern confidence intervals
for (i in seq_along(SNPids_n_mod)) {

  current_snp <- SNPids_n_mod[[i]]

  # Obtain the confidence interval for the current SNP from P3 (Modern)
  lowerCI <- ds4_n_mod$CI_Lower[
    ds4_n_mod$SNPid == current_snp &
      ds4_n_mod$Period == "P3"
  ][1]

  upperCI <- ds4_n_mod$CI_Upper[
    ds4_n_mod$SNPid == current_snp &
      ds4_n_mod$Period == "P3"
  ][1]

  # Select neutral SNPs whose P3 frequencies fall within the
  # confidence interval of the current SNP
  ds_rs_n_mod <- ds4_n_mod %>%
    select(
      frequency, Period, Neol_minus_Modern,
      SNPid, Phenotype, POPsize, CI_Lower, CI_Upper
    ) %>%
    group_by(SNPid) %>%
    filter(
      (Phenotype == "Neutral" &
         frequency >= lowerCI &
         frequency <= upperCI &
         Period == "P3") |
        SNPid == current_snp
    ) %>%
    ungroup()

  # Matched neutral reference distribution
  neutral_vals <-
    ds_rs_n_mod$Neol_minus_Modern[
      ds_rs_n_mod$Phenotype == "Neutral"
    ]

  # Skip the current SNP if no matched neutral SNPs are available
  if (length(neutral_vals) == 0) next

  p97_5 <- quantile(neutral_vals, 0.975, na.rm = TRUE)
  p2_5  <- quantile(neutral_vals, 0.025, na.rm = TRUE)

  # Frequency difference and phenotype of the current SNP
  snp_value <-
    ds_rs_n_mod$Neol_minus_Modern[
      ds_rs_n_mod$SNPid == current_snp
    ][1]

  phenotype <-
    unique(
      ds_rs_n_mod$Phenotype[
        ds_rs_n_mod$SNPid == current_snp
      ]
    )[1]

  # Identify SNPs falling outside the neutral distribution
  # Positive difference -> decrease toward the Modern period
  # Negative difference -> increase toward the Modern period
  if (snp_value > p97_5) {
    results2_n_m_filtered <- rbind(
      results2_n_m_filtered,
      data.frame(
        SNPid = current_snp,
        Status = "decrease",
        Phenotype = phenotype
      )
    )
  } else if (snp_value < p2_5) {
    results2_n_m_filtered <- rbind(
      results2_n_m_filtered,
      data.frame(
        SNPid = current_snp,
        Status = "increase",
        Phenotype = phenotype
      )
    )
  }
}

results2_n_m_filtered

write.csv(
  results2_n_m_filtered,
  "/Users/dilanur/Desktop/results_n_m_filtered_son.csv"
)
