args<-commandArgs(TRUE)
suppressMessages(library(HMMcopy))

sample.name <- args[1] # sample name
rfile <- args[2] # wig file produced by 03.readCounter.sh

gfile <- "/mnt/data/yanzhiqiang/138/software/HMMcopy/hg38/hg38.gc.1000000.wig" # gc file produced by hmmcopy toolkit
mfile <- "/mnt/data/yanzhiqiang/138/software/HMMcopy/hg38/hg38.map.1000000.bw" # map file produced by hmmcopy toolkit

normal_reads <- wigsToRangedData(rfile, gfile, mfile)
#normal_reads[100:110, ]

normal_copy <- correctReadcount(normal_reads)
#normal_copy[100:110, ]

final_segments <- HMMsegment(normal_copy)
#final_segments[10:110, ]

write.csv(normal_copy,file = paste0(sample.name,".hmmcopy.win1000000.corrected.csv"),row.names = F,quote = F)
write.csv(final_segments$segs,file = paste0(sample.name,".hmmcopy.win1000000.segment.csv"),row.names = F,quote = F)
