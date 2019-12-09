#!/usr/bin/env Rscript
source('myfunctions.R')

df <- read.table('/home/tb37w/project/Research/lonely_mice/master_loops/master_requested_loops', header=TRUE, sep="\t")
df <- clean_qvals(df)
write.table(df, '/home/tb37w/project/Research/lonely_mice/master_loops/clean_master_requested_loops', quote=FALSE, 
            sep="\t", row.names=FALSE)
