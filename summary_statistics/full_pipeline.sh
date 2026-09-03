#!/bin/bash

#SBATCH -p macaque
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J est_summ
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%j-%N-%u.err

# Paths
bed_dir="/mnt/NEOGENE1/projects/selection_2023/dilanur/summ_stats/P3/bed_files"
input_vcf="/mnt/NEOGENE1/projects/balancing_selection_2023/tgp_genotyping/TGP.realigned.hs37d5.sorted.rd.q20.wRG.vqsrfilter.bi5.2.onlyPASSautosomes.wbcftools.vcf.gz"
vcftools_path="/usr/local/sw/vcftools-0.1.16/bin/vcftools"

temp_output_dir="$(pwd)/temp_outputs"
flip_output_dir="$(pwd)/FLIP"
biallelic_output_dir="$(pwd)/Biallelic"
ms_output_dir="$(pwd)/ms_files"
summary_output_dir="$(pwd)/summary_statistics"

python_script="estimate_summary.py"
perl_script="vcf2ms.pl"

# Create output directories if they don't exist
mkdir -p "$temp_output_dir" "$flip_output_dir" "$biallelic_output_dir" "$ms_output_dir" "$summary_output_dir"

# Optional: clean previous partial outputs (comment out if you want to keep old results)
# rm -f "$temp_output_dir"/*.vcf "$flip_output_dir"/*.vcf "$biallelic_output_dir"/*.recode.vcf "$ms_output_dir"/*.ms "$summary_output_dir"/*.txt

# Step 1: Intersect BED files and create VCF files (keep header)
echo "Starting intersect operations for BED files..."
for bed_file in "${bed_dir}"/*.bed; do
  bed=$(basename "$bed_file" .bed)
  final_output="${temp_output_dir}/${bed}_P3.vcf"

  /usr/local/sw/bedtools2/bin/intersectBed \
    -header \
    -a "$input_vcf" \
    -b "$bed_file" \
    -wa > "$final_output"

  echo "Processed intersect for: ${bed}"
done

# Sanity check: Step 1 produced VCFs
if ! compgen -G "$temp_output_dir/*.vcf" > /dev/null; then
  echo "ERROR: Step 1 produced no VCF files in $temp_output_dir" >&2
  exit 1
fi

# Step 2: Flip genotypes using chimp/ancestor FASTA (Ensembl e71, GRCh37)
echo "Starting genotype flip operations (ancestral polarization via chimp FASTA)..."
for vcf_file in "$temp_output_dir"/*.vcf; do
  base_name=$(basename "$vcf_file" .vcf)
  output_file="${flip_output_dir}/${base_name}_flipped.vcf"
  python3 flip2.py "$vcf_file" "$output_file"
  echo "Flipped: $vcf_file -> $output_file"
done

# Sanity check: Step 2 produced flipped VCFs
if ! compgen -G "$flip_output_dir/*.vcf" > /dev/null; then
  echo "ERROR: Step 2 produced no flipped VCF files in $flip_output_dir" >&2
  exit 1
fi

# Step 3: Filter VCF files to retain only biallelic sites
echo "Starting biallelic filtering..."
for flip_file in "$flip_output_dir"/*.vcf; do
  filebase=$(basename "$flip_file" .vcf)
  out_prefix="${biallelic_output_dir}/${filebase}_biallelic"

  $vcftools_path --vcf "$flip_file" --min-alleles 2 --max-alleles 2 --recode --out "$out_prefix"

  echo "Biallelic filtered: $flip_file -> ${out_prefix}.recode.vcf"
done

# Sanity check: Step 3 produced recode VCFs
if ! compgen -G "$biallelic_output_dir/*_biallelic.recode.vcf" > /dev/null; then
  echo "ERROR: Step 3 produced no biallelic recode VCF files in $biallelic_output_dir" >&2
  exit 1
fi


# Step 4: Convert biallelic VCF files to MS format using vcf2ms.pl
echo "Starting VCF to MS conversion..."
sample_size=16  # Define sample size

for vcf_file in "$biallelic_output_dir"/*_biallelic.recode.vcf; do
  base_name=$(basename "$vcf_file" "_biallelic.recode.vcf")
  output_file="${ms_output_dir}/${base_name}.ms"

  # base_name example: rs988913_P3_flipped_biallelic
  # Our BED files are named as rsID only: rs988913.bed
  rsid="${base_name%%_*}"
  bed_file="${bed_dir}/${rsid}.bed"

  if [[ -f "$bed_file" ]]; then
    perl "$perl_script" "$vcf_file" "$output_file" "$sample_size" "$bed_file"
    echo "Converted: $vcf_file -> $output_file"
  else
    echo "Warning: BED file for ${rsid} not found at ${bed_file}, skipping..."
  fi
done


# Step 4.1: Remove empty lines after 'positions' in .ms files
echo "Removing empty lines after 'positions' in .ms files..."
python3 remove_empty_line.py "$ms_output_dir"
echo "Empty lines removed from MS files."

# Sanity check: Step 4 produced MS files
if ! compgen -G "$ms_output_dir/*.ms" > /dev/null; then
  echo "ERROR: No .ms files produced in $ms_output_dir" >&2
  exit 1
fi

# Step 5: Calculate summary statistics for MS files using estimate_summary.py
echo "Starting summary statistics calculations for .ms files..."
for ms_file in "$ms_output_dir"/*.ms; do
  echo "Processing $ms_file with RUN=1 and SAMPLING_TIME=1"
  python3 "$python_script" "$ms_file" 1 1 "$summary_output_dir"
done

echo "All operations completed successfully."
