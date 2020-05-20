#!/usr/bin/env python3

# Input: EHDN and EH output
# Output: Correlation between EHDN sizes and EH genotypes
# NOTE: It is unlikely that other users would want to run this script directly; this script is provided mainly to document the
# algorithm used to calculate the EHDN-EH correlation.

import argparse
import BTlib
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from collections import defaultdict
import os

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("EHDN_filename", type=str) # File containing EHDN motifs and TRF coordinates
parser.add_argument("genotype_table_filename", type=str) # EH genotype_table file
parser.add_argument("good_sample_list_filename", type=str) # file containing samples to use (as filtered by Bank). Any samples not in this file will be discarded.  = /hpf/largeprojects/tcagstor/users/worrawat/Expansion/SampleInfo/Final.17404.samples.MSSNG.SSC.1000G.tsv
parser.add_argument("correlation_function", type=str)
args = parser.parse_args()
#####################################

if args.correlation_function == "pearson":
    correlation_function = pearsonr
elif args.correlation_function == "spearman":
    correlation_function = spearmanr

output_dir = "{}.{}.out".format(args.EHDN_filename, args.correlation_function)

os.makedirs(output_dir, exist_ok=True)

good_samples = BTlib.get_set_from_file(args.good_sample_list_filename, col=0, header=True)
genotype_table = pd.read_csv(args.genotype_table_filename, sep="\t", index_col=0)

raw_output_file = open("{}/EHDN.overlap.correlations_raw.txt".format(output_dir), "w")
nice_output_file = open("{}/EHDN.overlap.correlations_clean.txt".format(output_dir), "w")
merged_region_output_file = open("{}/EHDN.overlap.correlations_merged.txt".format(output_dir), "w")

EHDN_file = open(args.EHDN_filename)

nice_output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("EHDN region", "EH region", "Number of individuals called by EHDN", "Correlation coefficient (largest allele)", "P-value (largest allele)", "Correlation coefficient (sum of alleles)", "P-value (sum of alleles)"))
merged_region_output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("EHDN region", "EH region", "Number of individuals called by EHDN", "Correlation coefficient (largest allele)", "P-value (largest allele)", "Correlation coefficient (sum of alleles)", "P-value (sum of alleles)"))


EHDN_list_merged = defaultdict(lambda: defaultdict(list))
EH_largest_allele_merged = defaultdict(lambda: defaultdict(list))
EH_sum_merged = defaultdict(lambda: defaultdict(list))

merged_EHDN_regions = []

seen = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: False)))

for line in EHDN_file:
    line = line.rstrip()
    fields = line.split("\t")

    merged_EHDN_start = fields[1]
    merged_EHDN_end = fields[2]

    orientation = fields[4]
    TRF_chrom = fields[6]
    TRF_start = fields[7]
    TRF_end = fields[8]
    EHDN_chrom = fields[11]
    EHDN_start = fields[12]
    EHDN_end = fields[13]
    EHDN_motif = fields[14]

    if orientation == "Reverse complement":
        EHDN_motif = BTlib.reverse_complement(EHDN_motif)

    merged_EHDN_region = "{}:{}-{}".format(TRF_chrom, merged_EHDN_start, merged_EHDN_end)
    if merged_EHDN_region not in merged_EHDN_regions:
        merged_EHDN_regions.append(merged_EHDN_region)

    TRF_region = "{}:{}-{}".format(TRF_chrom, TRF_start, TRF_end)

    if seen[merged_EHDN_region][TRF_region][EHDN_motif]:
        continue
    seen[merged_EHDN_region][TRF_region][EHDN_motif] = True

    individual_EHDN_region = "{}:{}-{}:{}".format(EHDN_chrom, EHDN_start, EHDN_end, EHDN_motif)

    genotype_table_key = "{}:{}-{}:{}".format(TRF_chrom, TRF_start, TRF_end, EHDN_motif)

    f = open("{}/{}.{}.list.txt".format(output_dir, individual_EHDN_region, genotype_table_key), "w")

    samples_field = fields[17]

    EHDN_list = []
    EH_largest_allele_list = []
    EH_sum_allele_list = []

    for item in samples_field.split(","):
        sample, EHDN_size = item.split(":")
        if sample not in good_samples:
            continue

        EH_genotype = genotype_table.loc[sample][genotype_table_key]

        if "." in EH_genotype or "-" in EH_genotype:
            continue
        if "/" in EH_genotype:
            allele1, allele2 = EH_genotype.split("/")
            largest_allele = max(int(allele1), int(allele2))
            sum_of_alleles = int(allele1) + int(allele2)
        else:
            largest_allele = int(EH_genotype)
            sum_of_alleles = int(EH_genotype)

        EHDN_list.append(float(EHDN_size))
        EH_largest_allele_list.append(largest_allele)
        EH_sum_allele_list.append(sum_of_alleles)

        EHDN_list_merged[merged_EHDN_region][TRF_region].append(float(EHDN_size))

        f.write("{:.3f}\t{}\n".format(float(EHDN_size), largest_allele))

        EH_largest_allele_merged[merged_EHDN_region][TRF_region].append(largest_allele)
        EH_sum_merged[merged_EHDN_region][TRF_region].append(sum_of_alleles)

    if len(EHDN_list) > 1:
        correlation_largest, correlation_largest_pvalue = correlation_function(EHDN_list, EH_largest_allele_list)
        correlation_sum, correlation_sum_pvalue = correlation_function(EHDN_list, EH_sum_allele_list)
    else:
        correlation_largest = 0
        correlation_largest_pvalue = 1
        correlation_sum = 0
        correlation_sum_pvalue = 1

    raw_output_file.write("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(line, len(EHDN_list), len(EHDN_motif), ":".join([str(x) for x in EHDN_list]), ":".join([str(x) for x in EH_largest_allele_list]), correlation_largest, correlation_sum))
    nice_output_file.write("{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(individual_EHDN_region, genotype_table_key, len(EHDN_list), correlation_largest, correlation_largest_pvalue, correlation_sum, correlation_sum_pvalue))

for merged_EHDN_region in merged_EHDN_regions:
    for TRF_region in EHDN_list_merged[merged_EHDN_region]:
        correlation_largest, correlation_largest_pvalue = correlation_function(EHDN_list_merged[merged_EHDN_region][TRF_region], EH_largest_allele_merged[merged_EHDN_region][TRF_region])
        correlation_sum, correlation_sum_pvalue = correlation_function(EHDN_list_merged[merged_EHDN_region][TRF_region], EH_sum_merged[merged_EHDN_region][TRF_region])
        merged_region_output_file.write("{}\t{}\t{}\t{:.3f}\t{:.1E}\t{:.3f}\t{:.1E}\n".format(merged_EHDN_region, TRF_region, len(EHDN_list_merged[merged_EHDN_region][TRF_region]), correlation_largest, correlation_largest_pvalue, correlation_sum, correlation_sum_pvalue))
