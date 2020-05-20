#!/usr/bin/env bash

# Input: a file containing motifs (repeat units), one per line (and which may be repeated). For example:
# AGT
# ACG
# AGT
# AA

# Output: Three plots giving the size distribution of the motifs, the most common motifs, and the base distribution of the motifs

repeat_units_file=$1

# Analyze repeats
EHDN_STR_summary.py $repeat_units_file

# Make plots from summary data
print_status "Creating repeat unit size distribution file...\n"
repeat_unit_size_distribution.R $repeat_units_file.repeat_unit_size_distribution.txt $repeat_units_file.repeat_unit_size_distribution.pdf

print_status "Creating repeat unit distribution file...\n"
repeat_unit_distribution.R $repeat_units_file.repeat_unit_distribution.txt $repeat_units_file.repeat_unit_distribution.pdf

print_status "Making base composition bar chart...\n"
base_composition_bar_chart.R $repeat_units_file.base_composition.txt $repeat_units_file.base_composition.pdf
