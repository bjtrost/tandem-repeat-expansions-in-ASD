#!/usr/bin/env python3

# HELPER SCRIPT - not to be called directly.
# This script is called by EHDN_STR_summary.sh


import argparse
import BTlib

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("repeat_unit_filename", type=str)
args = parser.parse_args()
#####################################

# Get the distribution of the sizes of the repeat units
def repeat_unit_size_distribution(repeat_units):
    output_file = open(args.repeat_unit_filename + ".repeat_unit_size_distribution.txt", "w")
    output_file.write("Size of motif (bp)\tNumber of tandem repeats\tPercentage of tandem repeats\n")

    repeat_unit_size_counts = BTlib.nested_dict(1, int)
    total_count = 0
    for repeat_unit in repeat_units:
        repeat_unit_size = len(repeat_unit)
        repeat_unit_size_counts[repeat_unit_size] += 1
        total_count += 1

    for size in sorted(repeat_unit_size_counts):
        percentage = "{:.2f}".format(repeat_unit_size_counts[size] / total_count * 100)
        output_file.write(f"{size}\t{repeat_unit_size_counts[size]}\t{percentage}\n")
    output_file.close()

def repeat_unit_distribution(repeat_units):
    output_file = open(args.repeat_unit_filename + ".repeat_unit_distribution.txt", "w")
    output_file.write("Motif\tNumber of tandem repeats\tPercentage of tandem repeats\n")

    repeat_unit_counts = BTlib.nested_dict(1, int)
    total_count = 0
    for repeat_unit in repeat_units:
        repeat_unit_counts[repeat_unit] += 1
        total_count += 1

    total_outputted = 0
    for repeat_unit in sorted(repeat_unit_counts, key=repeat_unit_counts.get, reverse=True):
        if total_outputted > 19:
            break
        percentage = "{:.2f}".format(repeat_unit_counts[repeat_unit] / total_count * 100)
        output_file.write(f"{repeat_unit}\t{repeat_unit_counts[repeat_unit]}\t{percentage}\n")
        total_outputted += 1
    output_file.close()

def base_distribution(repeat_units):
    output_file = open(args.repeat_unit_filename + ".base_composition.txt", "w")
    output_file.write("Motif\tPercentage A/T\n")

    for repeat_unit in repeat_units:
        percent_A_or_T = (repeat_unit.count("A") + repeat_unit.count("T")) / len(repeat_unit) * 100
        output_file.write("{}\t{:.1f}\n".format(repeat_unit, percent_A_or_T))

    output_file.close()

repeat_unit_file = open(args.repeat_unit_filename)
repeat_units = [x.rstrip("\n") for x in repeat_unit_file.readlines()]
repeat_unit_size_distribution(repeat_units)
repeat_unit_distribution(repeat_units)
base_distribution(repeat_units)
