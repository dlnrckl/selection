# Summary Statistics Pipeline

This repository contains the scripts used to preprocess 50-kb genomic windows and calculate population-genetic summary statistics from haplotype data.

## Workflow

For each target SNP, a 50-kb genomic window was defined around the target position (±25 kb). The analysis workflow was run from:

```text
/mnt/NEOGENE1/projects/selection_2023/dilanur/summ_stats_polarizasyon
```

The main workflow is implemented in `full_pipeline.sh` and consists of the following steps:

1. Intersect the target BED regions with the input VCF using BEDTools.
2. Polarize alleles according to the ancestral allele state using the Ensembl ancestral genome FASTA.
3. Retain biallelic variants using VCFtools.
4. Convert VCF files to MS format using the modified `vcf2ms.pl` script.
5. Remove the extra empty line after the `positions:` field in MS files.
6. Calculate summary statistics using `estimate_summary.py`.

## Files

### `full_pipeline.sh`

Main SLURM workflow that runs all preprocessing and summary-statistics steps sequentially.

### `make_beds.sh`

Creates one BED file for each target SNP/window from the input region file.

### `flip_chimp.py`

Polarizes VCF alleles using the ancestral allele state from the Ensembl ancestral genome FASTA. An earlier version of the workflow used 1000 Genomes data for allele-orientation harmonization; the final analysis used the chimpanzee-based ancestral FASTA, and this final script is provided as `flip_chimp.py`.

The ancestral FASTA directory used in the analysis was:

```text
/mnt/NEOGENE1/projects/selection_2023/dilanur/summ_stats_polarizasyon/chimp/homo_sapiens_ancestor_GRCh37_e71
```

For each SNP:

- if the ancestral allele matches REF, the variant is kept unchanged;
- if the ancestral allele matches ALT, REF and ALT are swapped and genotype allele codes are flipped;
- variants for which the ancestral state cannot be determined are retained unchanged.

### `vcf2ms.pl`

Modified version of the VCF-to-MS conversion script from:

```text
https://github.com/lstevison/vcf-conversion-tools/blob/master/vcf2MS.pl
```

The script converts the filtered VCF files to MS format and normalizes variant positions relative to the start of each 50-kb BED interval.

### `remove_empty_line.py`

Removes the extra empty line occurring after the `positions:` line in the generated MS files.

### `estimate_summary.py`

Modified version of the `Preprocess.py` script from the BaSe repository:

```text
https://github.com/ulasisik/balancing-selection/blob/master/BaSe/Preprocess.py
```

The script calculates the summary statistics used in the analysis.

The complete set calculated by the script includes:

- mean, median and maximum pairwise distance
- Tajima's D
- Watterson's theta
- observed heterozygosity
- observed/expected heterozygosity ratio
- median r²
- H1, H12, H123 and H2/H1
- haplotype diversity and number of haplotypes
- EHH
- iHS
- nSL
- NCD1
- Kelly's ZnS
- nucleotide diversity (pi)
- Fay and Wu's H
- number of singleton variants
- Fu and Li's D*
- Fu and Li's F*
- Zeng's E
- raggedness index

Downstream analyses focused primarily on site-frequency-spectrum-related statistics: Tajima's D, mean observed heterozygosity, Fu and Li's D*, Fu and Li's F*, and Zeng's E.

## Software

The workflow uses:

- BEDTools 2.25.0
- VCFtools 0.1.16
- Python 3
- Perl
- scikit-allel
- NumPy
- pandas
- pysam

## Example

The main pipeline was executed from:

```bash
cd /mnt/NEOGENE1/projects/selection_2023/dilanur/summ_stats_polarizasyon
sbatch full_pipeline.sh
```

Input/output paths and software paths used for the original analysis are defined directly in `full_pipeline.sh`.

## References

BaSe:

```text
https://github.com/ulasisik/balancing-selection
```

VCF conversion tools:

```text
https://github.com/lstevison/vcf-conversion-tools
```

Ensembl ancestral allele FASTA:

```text
https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
```

## Note

The scripts in this repository are provided in the form used for the analysis. Paths in the scripts are specific to the original computing environment and may need to be changed before running the workflow on another system.
