#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2)
   {
   print("Usage: $ ./processing-raw-csvs_merging-sets.R input_filename output_filename")
   stop("Wrong arguments.")
   }

input_filename  = args[1]
output_filename = args[2]

library('plyr', lib="/afs/cern.ch/user/o/otoldaie/work/soft/R-packages/")

distr <- read.csv(input_filename, header=T)

distr_per_set <- ddply( distr, .(type, dtag), function(df) {
    c <- colSums(df[4:length(df)])
    return(data.frame( bins=as.numeric(gsub("X", "", gsub("X.", "-", names(c), fixed=T))), value=as.numeric(c) ))
    })

write.csv( distr_per_set, output_filename, row.names=F)




