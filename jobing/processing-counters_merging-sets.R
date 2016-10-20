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

counters <- read.csv(input_filename, header=T)
# type,dtag,job_num,value

if (dim(counters)[1]==0) {
  counters_per_set <- counters[c(1,2,4)] # <---- !!!! all except job_num!! considering job_num is 3rd
} else {
  # sum all jobs
  counters_per_set <- ddply( counters, .(type, dtag), function(df) {
    r <- data.frame(sum(df[4]))
    names(r) <- names(df)[4]
    return(data.frame( r ))
    })
}

write.csv( counters_per_set, output_filename, row.names=F)




