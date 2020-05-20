#!/usr/bin/env python3

# Input: JSON file from ExpansionHunter Denovo
# Output: The number of raw calls (anchored in-repeat reasds) detected

import argparse
import BTlib

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("json_filename", type=str)
parser.add_argument("--ignore-alt-contigs", default=False, action="store_true")
args = parser.parse_args()
#####################################

chromosomes = BTlib.get_chromosomes()

f = open(args.json_filename)
sample_count = 0

for line in f:
    if "RegionsWithIrrAnchors" in line:
        counting = True
    elif "RegionsWithIrrs" in line:
        counting = False
    elif "-" in line and counting:
        chrom = line.split(":")[0].replace('"', "").replace(" ", "")
        if not args.ignore_alt_contigs or chrom.replace("chr", "") in chromosomes:
            sample_count += 1
f.close()

print(sample_count)
