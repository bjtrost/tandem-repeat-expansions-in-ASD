#!/usr/bin/env python3

# Input: File containing overlap between EHDN results, TRF loci, and AsmVar calls
# Output: Table indicating which EHDN calls are confirmed/not confirmed.
# NOTE: It is unlikely that other users would want to run this script directly; this script is provided mainly to document the
# algorithm used to calculate the EHDN PacBio validation

import argparse
import BTlib
from collections import defaultdict

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("TRF_overlap_filename", type=str)
args = parser.parse_args()
#####################################

def get_AsmVar_adjustment(AsmVar_matches, TRF_start, TRF_end):
    AsmVar_adjustment = 0

    if TRF_start == "NA":
        TRF_start_padded = "NA"
        TRF_end_padded = "NA"
    else:
        TRF_start_padded = int(TRF_start) - 100
        TRF_end_padded = int(TRF_end) + 100

    for AsmVar_match in AsmVar_matches:
        if AsmVar_match == "NA" or AsmVar_match[:3] == "tig":
            continue

        _, AsmVar_coords, AsmVar_size, AsmVar_type = AsmVar_match.split(":")
        AsmVar_start, AsmVar_end = [int(x) for x in AsmVar_coords.split("-")]
        AsmVar_size = int(AsmVar_size)

        if AsmVar_type == "INS" or AsmVar_type == "DEL":
            if TRF_start == "NA" or BTlib.overlaps(AsmVar_start, AsmVar_end, TRF_start_padded, TRF_end_padded):
                AsmVar_adjustment += AsmVar_size

    return(AsmVar_adjustment)


###################################################################
###################################################################

TRF_overlap_file = open(args.TRF_overlap_filename)

results = defaultdict(lambda: defaultdict(lambda: defaultdict("N/A")))
print("EHDN chr\tEHDN start\tEHDN end\tEHDN motif\tCanu/AsmVar variants\tTRF match\tTRF size\tCanu/AsmVar adjustment\tTotal size\tConfirmed?")
keys = []
caller_end_index = 5

for line in TRF_overlap_file:
    fields = line.rstrip("\n").split("\t")
    key = "\t".join(fields[:4])
    if key not in keys:
        keys.append(key)
        results[key]["LARGEST_TOTAL_SIZE"] = -999999
    results[key]["LINE"] = "\t".join(line.rstrip("\n").split("\t")[:caller_end_index])
    EHDN_repeat_unit = fields[3]

    TRF_match = defaultdict(lambda: "NA")
    TRF_match["TRF_size"] = 0
    if fields[caller_end_index] != ".":
        TRF_match["chr"], TRF_match["start"], TRF_match["end"], TRF_match["repeat_unit"] = fields[caller_end_index:caller_end_index+4]
        TRF_match["TRF_size"] = int(TRF_match["end"]) - int(TRF_match["start"]) + 1

    AsmVar_matches = fields[4].split("|")
    TRF_match["AsmVar_adjustment"] = get_AsmVar_adjustment(AsmVar_matches, TRF_match["start"], TRF_match["end"])
    TRF_match["total_size"] = TRF_match["TRF_size"] + TRF_match["AsmVar_adjustment"]

    if TRF_match["total_size"] >= 150:
        TRF_match["confirmed"] = "yes"
    else:
        TRF_match["confirmed"] = "no"

    if TRF_match["total_size"] > results[key]["LARGEST_TOTAL_SIZE"]:
        results[key]["LARGEST_TOTAL_SIZE"] = TRF_match["total_size"]
        results[key]["MATCH"] = TRF_match

for key in keys:
    m = results[key]["MATCH"]
    TRF_match = "{}:{}-{}:{}".format(m["chr"], m["start"], m["end"], m["repeat_unit"])

    print("\t".join([results[key]["LINE"], TRF_match, str(m["TRF_size"]), str(m["AsmVar_adjustment"]), str(m["total_size"]), m["confirmed"]]))
