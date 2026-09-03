import gzip
import re
import argparse
import os
import pysam

ANCESTOR_DIR = "/mnt/NEOGENE1/projects/selection_2023/dilanur/summ_stats_polarizasyon/chimp/homo_sapiens_ancestor_GRCh37_e71"

def read_vcf(file_path):
    opener = gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'rt')
    with opener as f:
        lines = f.readlines()
    meta    = [l.rstrip('\n') for l in lines if l.startswith('##')]
    header  = next(l.rstrip('\n') for l in lines if l.startswith('#') and not l.startswith('##'))
    records = [l.rstrip('\n') for l in lines if not l.startswith('#') and l.strip()]
    return meta, header, records

def flip_genotype(gt):
    sep = '|' if '|' in gt else '/'
    alleles = gt.split(sep)
    flipped = ['1' if a == '0' else '0' if a == '1' else a for a in alleles]
    return sep.join(flipped)

def get_ancestral(fa_cache, chrom, pos):
    chrom_num = chrom.replace('chr', '')
    fa_path   = os.path.join(ANCESTOR_DIR, f"homo_sapiens_ancestor_{chrom_num}.fa")

    if not os.path.exists(fa_path):
        return 'N'

    if fa_path not in fa_cache:
        try:
            fasta = pysam.FastaFile(fa_path)
            # header: >ANCESTOR_for_chromosome:GRCh37:1:1:...:1
            seq_name = fasta.references[0]
            fa_cache[fa_path] = (fasta, seq_name)
        except Exception:
            return 'N'

    fasta, seq_name = fa_cache[fa_path]

    try:
        return fasta.fetch(seq_name, pos - 1, pos).upper()
    except Exception:
        return 'N'

def process(input_vcf, output_vcf):
    meta, header, records = read_vcf(input_vcf)
    fa_cache = {}
    output_records = []
    stats = {'flipped': 0, 'kept': 0, 'skipped_N': 0, 'skipped_mismatch': 0, 'non_snp': 0}

    for record in records:
        fields = record.split('\t')
        chrom  = fields[0]
        pos    = int(fields[1])
        ref    = fields[3]
        alt    = fields[4]

        # sadece biallelic SNP
        if len(ref) != 1 or len(alt) != 1 or ',' in alt:
            stats['non_snp'] += 1
            output_records.append('\t'.join(fields))
            continue

        aa = get_ancestral(fa_cache, chrom, pos)

        if aa == 'N':
            stats['skipped_N'] += 1
            output_records.append('\t'.join(fields))

        elif aa == ref.upper():
            # REF zaten ancestral, flip gerekmez
            stats['kept'] += 1
            output_records.append('\t'.join(fields))

        elif aa == alt.upper():
            # ALT ancestral -> REF/ALT swap + genotipleri flip
            fields[3] = alt
            fields[4] = ref
            for i in range(9, len(fields)):
                gt_parts    = fields[i].split(':')
                gt_parts[0] = flip_genotype(gt_parts[0])
                fields[i]   = ':'.join(gt_parts)
            stats['flipped'] += 1
            output_records.append('\t'.join(fields))

        else:
            # AA ne REF ne ALT -> polarize edilemiyor, oldugu gibi birak
            stats['skipped_mismatch'] += 1
            output_records.append('\t'.join(fields))

    with open(output_vcf, 'w') as f:
        for line in meta:
            f.write(line + '\n')
        f.write(header + '\n')
        for record in output_records:
            f.write(record + '\n')

    total = len(records)
    print(f"  Total SNPs     : {total}")
    print(f"  Kept (REF=AA)  : {stats['kept']}")
    print(f"  Flipped        : {stats['flipped']}")
    print(f"  Skipped (N)    : {stats['skipped_N']}")
    print(f"  Skipped (mismatch): {stats['skipped_mismatch']}")
    print(f"  Non-SNP/multi  : {stats['non_snp']}")

def main():
    parser = argparse.ArgumentParser(
        description="Flip VCF REF/ALT using Ensembl ancestor FASTA so that REF=ancestral"
    )
    parser.add_argument('input',  help="Input VCF (plain or .gz)")
    parser.add_argument('output', help="Output flipped VCF")
    args = parser.parse_args()
    process(args.input, args.output)

if __name__ == "__main__":
    main()

### usage: python flip_chimp.py input.vcf output_flipped.vcf
