rm(list = ls())

############## plot cnv plot of representative cells in fig2
#E09_052, euploid
#E19_038, gain
#E15_027(real is E01_051), loss

############## plot cnv plot of complementary cells in fig2
#E03_069,E03_081, whole
#E06_026,E06_029, seg


evalue=0.9999999;

samplename = 'E19_038'

library(HMMcopy)
rfile <- paste0("/home/kongsiming/project/PGT-A/wig/",samplename,".window1000000.readcounts.wig")
gfile <- "/home/yanzhiqiang/software/HMMcopy/hg38/hg38.gc.1000000.wig" #gc file
mfile <- "/home/yanzhiqiang/software/HMMcopy/hg38/hg38.map.1000000.bw" # mappbility file
x1_untrimmed_reads <- wigsToRangedData(rfile,gfile,mfile)
corrected <- correctReadcount(x1_untrimmed_reads)

param <- HMMsegment(corrected, getparam = TRUE)
param$e = evalue

segments <- HMMsegment(corrected, param = param )

chromosomes <- levels(corrected$space) # List of chromosomes...

lengths <- aggregate(end(corrected), by = list(chr = space(corrected)), max) # Obtain length of each chromosome
lengths <- lengths[c(1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 11, 13, 14, 15, 16, 17, 18, 19, 21, 20, 24, 23, 8, 22), ] # Just sorting
lengths$offset <- c(0, cumsum(as.numeric(lengths$x))[-nrow(lengths)]) # Obtain the cumulative sum of chromosomes lengths
corrected$state <- segments$state

# Merging offset information back into original corrected data frame
data <- merge(corrected, lengths[, c(1, 3)], by.x = "space", by.y = "chr") # Creates an offset column, with values assigned according chromosomes


if(1<0){
  data1 <- data.frame(data)
  head(data1)
  data1[data1$space == "chr4", 13] <- 3
}


pdf(paste0(samplename,".dot.plot.pdf"), height = 4, width = 20)
# Reproducing internals of plotSegment for maximum customization
cols <- stateCols() # Predefined colours
cols <- c("green","#00008B","#00008B","#00008B","red","red")

range <- quantile(data$copy, na.rm = TRUE, prob = c(0.01, 0.99)) # Optional range, tweak this as desired

plot(x = start(data$ranges) + data$offset, y = data$copy,xaxt='n',main=samplename, xlab="", ylab="", col = cols[data$state], pch = 20, ylim = c(-1.5,1.5) )
#plot(x = start(data$ranges) + data$offset, y = data$copy,xaxt='n',main=samplename, xlab="", ylab="", col = cols[data1$state], pch = 20, ylim = c(-1.5,1.5) )

abline(v = lengths$offset) # Separate chromosomes

#segments(x0 = min(start(data$ranges) + data$offset), y0 = 0.44, x1 = max(start(data$ranges) + data$offset), y1 = 0.44, col = "red", lwd = 2) #
#segments(x0 = min(start(data$ranges) + data$offset), y0 = -0.57, x1 = max(start(data$ranges) + data$offset), y1 = -0.57, col = "green", lwd = 2) #

#text(lengths$offset, 1, labels = lengths$chr, adj = -0.1) # Label chromosomes
text(lengths$offset, 1.5, labels = gsub("chr", "", lengths$chr), adj = -0.2) # Label chromosomes

# Inserting the segmentation lines if desired...
segs <- merge(segments$segs, lengths[, c(1, 3)], all.x = TRUE)
segments(x0 = segs$start + segs$offset, y0 = segs$median, x1 = segs$end + segs$offset, col = "black", lwd = 3) # Note, segments is a function... not the output, unfortunately naming...
dev.off()

