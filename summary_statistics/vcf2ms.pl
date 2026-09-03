#! /usr/bin/perl

#program converts vcf file to MS input file format
#Modified to normalize positions based on a starting position from a BED file
#and replace './.' genotypes with '0'

use strict;
use warnings;

# Input files and sample size from command line arguments
my $vcf = $ARGV[0];
my $output = $ARGV[1];
my $sample_size = $ARGV[2];
my $bed_file = $ARGV[3]; # BED file containing the starting position

unless ($#ARGV==3) {
    print STDERR "Please provide name of input vcf file, filename for MS formatted output, sample size, and BED file on command line\n\n";
    die;
} #end unless

# Read the starting position from the BED file
open(BED, $bed_file) or die "Could not open BED file: $!";
my $bed_line = <BED>;  # Read the first line from the BED file
chomp($bed_line);
close(BED);

# Extract the starting position from the BED file line (second column)
my @bed_fields = split(/\t/, $bed_line);
my $start_position = $bed_fields[1]; # Second column is the start position (index 1)

open(VCF, $vcf);

my @positions = ();
my @names = ();
my $loop_size = $sample_size + 8;
my %genotypes = ();

print STDERR "Reading in VCF file...";

while(<VCF>) {
    chomp;
    if ($_=~/\#\#/) {
        next;
    } elsif ($_=~/\#/) {
        my @input_line = split(/\s+/, $_);
        for (my $a=9; $a<=$loop_size; $a++) {
            push(@names, $input_line[$a]);
        } #end for
        next;
    } #end elsif
    
    my @input_line = split(/\s+/, $_);
    push(@positions, $input_line[1]);

    for (my $i=9; $i<=$loop_size; $i++) {
        my $o = $i - 9;
        my $hap1 = $names[$o] . "_1";
        my $hap2 = $names[$o] . "_2";
        my @genotype = split(":", $input_line[$i]);
        
        # Check for './.' and replace with '0'
        if ($genotype[0] eq "./.") {
            push @{$genotypes{$hap1}}, 0;
            push @{$genotypes{$hap2}}, 0;
        } else {
            my @haplotypes = split(/[\|\/]/, $genotype[0]);
            push @{$genotypes{$hap1}}, $haplotypes[0];
            push @{$genotypes{$hap2}}, $haplotypes[1];
        } #end if
    } #end for
} #end while

print STDERR "done.\nNow printing output...";

open(OUTPUT, ">$output");
my $number_loci = $#positions + 1;

# Print header
print OUTPUT "//\n";
print OUTPUT "segsites: $number_loci\n";
print OUTPUT "positions: ";

# Print normalized positions
for (my $b=0; $b<=$#positions; $b++) {
    my $normalized_position = ($positions[$b] - $start_position) / 50000;
    print OUTPUT "$normalized_position ";
} #end for
print OUTPUT "\n";

# Print genotypes
for (my $c=0; $c<=$#names; $c++) {
    print OUTPUT "\n";

    my $hap1 = $names[$c] . "_1";
    my @hap1_geno = @{$genotypes{$hap1}};

    for (my $d=0; $d<=$#hap1_geno; $d++) {
        print OUTPUT "$hap1_geno[$d]";
    } #end for
    print OUTPUT "\n";

    my $hap2 = $names[$c] . "_2";
    my @hap2_geno = @{$genotypes{$hap2}};

    for (my $e=0; $e<=$#hap2_geno; $e++) {
        print OUTPUT "$hap2_geno[$e]";
    } #end for
} #end for

print OUTPUT "\n";

print STDERR "done.\n";
