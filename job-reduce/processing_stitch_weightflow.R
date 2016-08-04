#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#length(args)
#dim(args)
#mode(args)
#class(args)
#str(args)
#print(args)

# args contains all filenames of per-set weight-flow values to be merged


output_filename <- args[1]

weight_flow <- read.csv(args[2], header=T)
dtags <- weight_flow$dtag

for (f in args[3:length(args)])
   {
   d <- read.csv(f, header=T)
   # first 3 columns are job definition: (MC/Data, dtag, job_number)
   # the 4th column is the value
   # 
   #n <- names(d)[4]
   #v <- d[4]
   #names(v) <- d$dtag
   weight_flow <- merge(weight_flow, d, all=T)
   }

weight_flow[is.na(weight_flow)] <- 0

write.csv( weight_flow, output_filename, row.names=F)




