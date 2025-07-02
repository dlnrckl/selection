# Variant Frequency Calculation Workflow

This folder hosts a **three‑step pipeline** that starts with variant calling, enriches the resulting table with phenotype & grouping metadata, and finally computes per‑group risk/minor allele frequencies together with Wilson 95 % confidence intervals.

> **Workflow Overview**
>
> 1. **`variant_call.sh`** → Generates per‑sample VCFs from BAM files.
> 2. **`process_variant_table.R`** → Adds phenotype, risk/minor allele, and group columns to the combined variant table.
> 3. **`freq_calculation.R`** → Calculates allele frequencies for every *Group × SNP × Phenotype* combination.

---

## 1 | Prerequisites

| Software     | Version / Notes                                                                                     |
| ------------ | --------------------------------------------------------------------------------------------------- |
| **bcftools** |  ≥ 1.18 — binary path defined by the variable `bcftools` inside `variant_call.sh`.                  |
| **R**        |  ≥ 4.2 — required packages: `tidyverse`, `data.table`, `dplyr`, `ggpubr`, `Hmisc`, `glue`, `binom`. |
| **SLURM**    |  (optional) — `variant_call.sh` ships with SBATCH lines; remove them for local execution.           |

Additional inputs

* **Reference genome** (e.g. hg38 FASTA)
* **Target BED file** specifying the regions used in mpileup
* **Sample BAM list** (`alkan_samples.txt`, one BAM name per line)
* **Risk\_allele table**  (columns: *Chromosome, POS, Risk\_Allele, Minor\_Allele, …, Phenotype*)

---

## 2 | Folder Layout

```
frequency_calculation/
├── variant_call.sh
├── process_variant_table.R
├── freq_calculation.R
└── vcf_files/          # generated automatically; holds the VCF outputs
```

---

## 3 | Step‑by‑Step Execution

\### 3.1  Variant Calling — `variant_call.sh`

```bash
# Submit to SLURM (cluster)
sbatch variant_call.sh

# OR run locally (after removing SBATCH lines)
bash variant_call.sh
```

What the script does:

1. Defines key variables: `bcftools`, `ref`, `bedfile`, `INFILE`, `OUTDIR`.
2. Creates `OUTDIR` (`vcf_files/`) if it doesn’t exist.
3. For every sample listed in `INFILE`:

   * Runs `bcftools mpileup` (restricted to the BED regions, including DP & AD tags).
   * Pipes to `bcftools call -mv` → produces `${sample}.vcf`.
4. Adds sample name to each VCF header via `awk`.
5. Merges all modified VCFs into **`combined.vcf`**.

> **Outputs**: `vcf_files/*.vcf` and `combined.vcf`

\### 3.2  Table Enrichment — `process_variant_table.R`

```R
Rscript process_variant_table.R
```

Main steps:

1. Reads the raw CSV exported from `combined.vcf` (columns such as `CHROM, POS, DP, DP4, REF, ALT, BamID`).
2. Creates a unique coordinate key (`chr_pos`).
3. Joins with the *Risk\_allele* table to append `SNP_ID`, `Risk_Allele`, `Minor_Allele`, `Phenotype`.
4. Splits the `DP4` field into four depth columns (`Ref1, Ref2, Alt1, Alt2`).
5. Calculates **risk/minor allele read count** (`Total_RAC`).
6. Adds **Group** information (default set to `"Modern"`; change to suit your study design).
7. Writes an R data frame object called `modern_df` (columns: `SampleName, Group, POS, SNP_ID, Phenotype, R, T`).

> **Note**: Paths and object names (e.g. `comperative_alkan`, `Risk_allele`) are hard‑coded; modify them if your filenames differ.

\### 3.3  Frequency Calculation — `freq_calculation.R`

```R
Rscript freq_calculation.R
```

1. Loads `modern_df` (and optionally other group tables such as `comp_others4`).
2. For each *Group × SNP\_ID × Phenotype* subset:

   * Estimates allele frequency **p̂** by Maximum Likelihood.
   * Computes Wilson 95 % confidence interval using `binom.confint`.
3. Saves the results to **`all_freq.csv`** (default location: Desktop). Columns include:

   * `Period`, `SNPid`, `Phenotype`, `POPsize`, `pHat`, `CI_Lower`, `CI_Upper`, …

---

## 4 | Key Output Files

| File           | Description                                                    |
| -------------- | -------------------------------------------------------------- |
| `combined.vcf` | Multi‑sample VCF containing all called variants                |
| `modern_df`    |  R object — enriched variant table with phenotype & group info |
| `all_freq.csv` | Final per‑group allele frequencies with Wilson CIs             |

---

## 5 | Customization Tips

* **SLURM directives** — adjust `--partition`, `--time`, and resource flags to match your cluster.
* **Region selection** — point `bedfile` to a different BED to focus on other genomic intervals.
* **Group assignment** — edit the line `modern$Group = "Modern"` in `process_variant_table.R` to reflect archaeological periods, populations, etc.
* **Parameterization** — the R scripts currently rely on hard‑coded paths; consider adding argument parsing via `commandArgs(trailingOnly = TRUE)`.

---

## 6 | References

* **bcftools** — Danecek et al. 2011. *Bioinformatics* 27(21): 2987‑2993.
* **Wilson CI** — Wilson, E.B. 1927. *J. Amer. Stat. Assoc.* 22: 209‑212.
