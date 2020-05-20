#!/usr/bin/env Rscript

# HELPER SCRIPT - not to be called directly.
# This script is called by EHDN_STR_summary.sh

suppressMessages(library(BTlib))
suppressMessages(library(extrafont))
suppressMessages(library(ggplot2))

args = commandArgs(TRUE)
file = args[1]
outfile = args[2]

data = read.table(file, header=TRUE, check.names=FALSE, sep="\t")

histogram(
    data=data,
    column="Percentage A/T",
    outfile=outfile,
    ylab="Percentage of motifs",
    num_bins=40,
    x_axis_hjust=0.5,
    center=50,
    text_size=30
)
