library(ggplot2)
library("reshape2")

args <- commandArgs(TRUE)
if (length(args) < 1){
    cat("Enter argument")
}

readVal <- read.csv(args[1],stringsAsFactors=FALSE,na.string=".")
samples <- unique(readVal$sample)

#cat("GQ == 99")
#print(NROW(readval[readVal$GQ == 99,]))

#Count samples
#for (sample in samples){
#   print(sample)
#   print(NROW(readVal[readVal$sample == sample && readVal$GQ == 99,]))
#}




