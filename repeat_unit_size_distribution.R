#!/usr/bin/env Rscript

# HELPER SCRIPT - not to be called directly.
# This script is called by EHDN_STR_summary.sh

suppressMessages(library(BTlib))
suppressMessages(library(extrafont))

args = commandArgs(TRUE)
file = args[1]
outfile = args[2]
data = read.table(file, header=TRUE, check.names=FALSE, sep="\t")

data[,1] = as.character(data[,1])
data[,1] = factor(data[,1], levels=data[,1])

bar_graph(
    data=data,
    xcol="Size of motif (bp)",
    ycol="Percentage of tandem repeats",
    outfile=outfile,
    x_text_angle=0,
    ylim=c(0,40),
    continuous_y_breaks=c(0,10,20,30),
    x_axis_hjust=0.5,
    text_size=30,
    x_text_size=16
    )
