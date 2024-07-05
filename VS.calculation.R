rm(list = ls())
sammples <- read.table("../samples.txt",header = F,stringsAsFactors = F)
# samples.txt records the name of each cell

sd.out <- c()
for(s in sammples$V1){
  df <- read.csv(paste0("win1M.hmmcopy/",s,".hmmcopy.win1000000.corrected.csv"),header = T,stringsAsFactors = F) # this is the corrected reads file produced by 04.hmmcopy.R
  df <- df[!is.na(df$copy),]
  aggdata <- aggregate(df[,11],by=list(V1=df$chr),FUN=sd)
  sd.now <- round(median(aggdata$x),4)
  sd.out <- c(sd.out,sd.now)
}

data.vs <- data.frame(sammples$V1,sd.out)
colnames(data.vs) <- c("samples","vs")
write.table(data.vs,file = "1M.sd.xls",quote = F,sep = "\t",col.names = T,row.names = F)
