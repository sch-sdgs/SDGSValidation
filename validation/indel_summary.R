#!/usr/bin/Rscript

args<-commandArgs(TRUE)
dir <- args[1]

final_stats <- paste(dir,"final_stats", sep = "/")
coverage <- paste(dir,"coverage", sep = "/")
depth_tp <- paste(dir,"depth_for_true_pos", sep = "/")
qual_tp <- paste(dir,"quality_for_true_pos", sep = "/")
qual_fp <- paste(dir,"quality_for_false_pos", sep = "/")
out_file <- paste(dir,"indel_summary.html", sep="/")

require('ggplot2')
require('reshape2')
require('rmarkdown')
rmarkdown::render("/home/bioinfo/Natalie/wc/sdgsValidation/validation/indel_comp_summary.Rmd",
                    params=list(final_stats=final_stats,cov=coverage,depth_tp=depth_tp,qual_tp=qual_tp,wual_fp=qual_fp),
                    output_file=out_file)


