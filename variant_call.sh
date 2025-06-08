#!/bin/bash
#SBATCH --job-name=var_call    
#SBATCH --partition=macaque   
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=1   
#SBATCH --time=30-00:00:00    
#SBATCH --output=output_%j.txt 
#SBATCH --error=error_%j.txt

bcftools=/usr/local/sw/bcftools-1.18/bcftools
ref=/mnt/NEOGENE3/share/ref/genomes/hsa/hg38-v0-Homo_sapiens_assembly38.fasta
bedfile=/mnt/NEOGENE1/projects/selection_2023/dilanur/P3/50kb_hg38.bed

# File containing sample names (one sample per line)
INFILE=/mnt/NEOGENE1/projects/selection_2023/dilanur/P3/yeni_P3/alkan_ornekler.txt

OUTDIR=vcf_files
mkdir -p ${OUTDIR}

echo "Starting variant calling process..."

# Loop through sample names and perform variant calling for each sample
while IFS=$'\n' read -r sample
do
    filein=/mnt/NEOGENE3/share/dna/hsa/comparative_seqs/alkan2014/${sample}
    echo "Processing: ${sample}"
    $bcftools mpileup -R -B -q 30 -Q 30 -f ${ref} -a FORMAT/DP,FORMAT/AD -r ${bedfile} ${filein} | \
    $bcftools call -mv -a GQ,GP -Ov -o ${OUTDIR}/${sample}.vcf
done < "$INFILE"

echo "Adding sample names to VCF files..."

# For each VCF file, append the sample name to every line and create a modified file
for file in ${OUTDIR}/*.vcf; do
    filename=$(basename "$file")
    filename_without_extension="${filename%%.*}"
    awk -v filename="$filename_without_extension" '{print $0 "\t" filename}' "$file" > "${file}.modified"
done

echo "Combining all modified VCF files..."
cat ${OUTDIR}/*.vcf.modified > combined.vcf

echo "All steps completed."
