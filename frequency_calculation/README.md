# Variant Frequency Calculation Workflow

This folder contains the workflow used to generate variant-level allele-frequency estimates, annotate variants with phenotype information, and identify phenotype-associated SNPs whose temporal allele-frequency changes deviate from matched neutral expectations.

## Workflow Overview

1. **`variant_call.sh`** → Generates per-sample VCFs from BAM files.
2. **`process_variant_table.R`** → Adds phenotype, risk/minor allele, and group information to the combined variant table.
3. **`freq_calculation.R`** → Calculates allele frequencies and Wilson 95% confidence intervals for each Group × SNP × Phenotype combination.
4. **`identify_snps_deviating_from_neutral.R`** → Compares temporal allele-frequency changes of phenotype-associated SNPs with frequency-matched neutral SNP distributions.

---

## 1 | Prerequisites

| Software | Version / Notes |
| --- | --- |
| **bcftools** | ≥ 1.18 — binary path defined inside `variant_call.sh`. |
| **R** | ≥ 4.2 — required packages include `tidyverse`, `data.table`, `dplyr`, `ggpubr`, `Hmisc`, `glue`, and `binom`. |
| **SLURM** | Optional — `variant_call.sh` contains SBATCH directives for cluster execution. |

Additional inputs:

- Reference genome FASTA
- Target BED file
- Sample BAM list (`alkan_samples.txt`)
- Risk/minor allele annotation table
- Phenotype and population/group metadata

---

## 2 | Folder Layout

```text
frequency_calculation/
├── variant_call.sh
├── process_variant_table.R
├── freq_calculation.R
├── identify_snps_deviating_from_neutral.R
├── all_data_df_example.csv
└── vcf_files/
```

`all_data_df_example.csv` contains the first 50 rows of the dataset used in the neutral-comparison analysis. It is provided only to illustrate the expected input structure and column names and does not represent the complete dataset.

---

## 3 | Step-by-Step Execution

### 3.1 Variant Calling — `variant_call.sh`

```bash
sbatch variant_call.sh
```

or, after removing the SLURM directives:

```bash
bash variant_call.sh
```

The script:

1. Defines the `bcftools`, reference genome, BED file, sample list, and output directory.
2. Runs `bcftools mpileup` for each sample within the target regions.
3. Calls variants using `bcftools call`.
4. Produces individual VCF files.
5. Adds sample identifiers to the VCF headers.
6. Merges the individual VCF files into `combined.vcf`.

**Outputs:**

```text
vcf_files/*.vcf
combined.vcf
```

---

### 3.2 Variant Table Processing — `process_variant_table.R`

```bash
Rscript process_variant_table.R
```

Main steps:

1. Reads the variant table derived from `combined.vcf`.
2. Creates a genomic coordinate identifier.
3. Joins the variant table with risk/minor allele and phenotype annotations.
4. Splits the `DP4` field into reference and alternative read counts.
5. Calculates risk/minor allele read counts.
6. Adds population or temporal-group information.
7. Creates the processed data table used for downstream frequency calculations.

Some input paths and object names are specific to the original analysis environment and may need to be changed before reuse.

---

### 3.3 Allele-Frequency Calculation — `freq_calculation.R`

```bash
Rscript freq_calculation.R
```

For each population/group, SNP, and phenotype combination, the script:

1. Estimates allele frequency.
2. Calculates Wilson 95% confidence intervals.
3. Generates the frequency table used in downstream temporal comparisons.

The resulting data include variables such as:

```text
Period
SNPid
Phenotype
POPsize
frequency / pHat
CI_Lower
CI_Upper
```

The final frequency table was subsequently used to compare phenotype-associated variants with neutral SNPs.

---

### 3.4 Identification of SNPs Deviating from Neutral Expectations

`identify_snps_deviating_from_neutral.R` identifies phenotype-associated SNPs whose temporal allele-frequency changes fall outside the distribution observed among frequency-matched neutral SNPs.

Three temporal comparisons are performed:

- Neolithic (P1) vs. Late Chalcolithic (P2)
- Late Chalcolithic (P2) vs. Modern (P3)
- Neolithic (P1) vs. Modern (P3)

For each phenotype-associated SNP:

1. Neutral SNPs with allele frequencies falling within the corresponding confidence interval are selected.
2. The temporal allele-frequency differences of these matched neutral SNPs are used to construct a reference distribution.
3. The 2.5th and 97.5th percentiles of the neutral distribution are calculated.
4. Phenotype-associated SNPs falling outside these limits are classified as showing an increase or decrease relative to neutral expectations.

Example outputs include:

```text
results_n_lc_filtered_son.csv
results_lc_m_filtered_son.csv
results_n_m_filtered_son.csv
```

---

## 4 | Example Input Data

`all_data_df_example.csv` contains the first 50 rows of the input dataset used by `identify_snps_deviating_from_neutral.R`.

The example file is included to demonstrate the expected structure of the input table, including columns such as:

```text
frequency
Period
SNPid
Phenotype
POPsize
CI_Lower
CI_Upper
```

It is not the complete analysis dataset.

---

## 5 | Key Output Files

| File | Description |
| --- | --- |
| `combined.vcf` | Combined variant calls across samples |
| Processed variant table | Variant table enriched with phenotype, allele, and group information |
| `all_freq.csv` | Allele-frequency estimates and confidence intervals |
| `results_n_lc_filtered_son.csv` | SNPs deviating from neutral expectations between P1 and P2 |
| `results_lc_m_filtered_son.csv` | SNPs deviating from neutral expectations between P2 and P3 |
| `results_n_m_filtered_son.csv` | SNPs deviating from neutral expectations between P1 and P3 |

---

## 6 | Notes

Some file paths, object names, and population labels in the scripts are specific to the original analysis environment and may need to be modified before running the workflow on another system.

The example dataset included in this repository is provided only to illustrate the input format. The complete study dataset is not included.

Visualization scripts used to generate manuscript figures are not included in this folder; the repository focuses on the analytical steps used to generate the underlying results.

---

## 7 | References

- **bcftools** — Danecek et al. 2011. *Bioinformatics* 27(21): 2987–2993.
- **Wilson confidence interval** — Wilson, E.B. 1927. *Journal of the American Statistical Association* 22: 209–212.
