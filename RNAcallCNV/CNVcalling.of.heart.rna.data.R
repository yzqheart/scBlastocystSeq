## Aneuploid detection
rm(list = ls())
library(readxl)
library(data.table)
library(edgeR)
library(dplyr)
library(scploid)
library(cowplot)
library(ggplot2)
library(ggsci)

here <- here::here 
summarize <- dplyr::summarize  
## Read in metadata sheet
metasheet <- read_xlsx("example/heart.metadata.xlsx", sheet = 1)
head(metasheet)

metasheet <- metasheet %>%
  as.data.table()
metasheet <- metasheet[order(metasheet$Sample),]
metasheet$Group <- paste0(metasheet$Estage, "_", metasheet$lineage)   

metasheet[, qc_passed := (`Percent Mapped` > quantile(`Percent Mapped`, 0.1) & `Mapped Reads` > quantile(`Mapped Reads`, 0.1))]
metasheet[, Embryo := gsub("_", ".", Embryo)]

## read in count matrix
raw_count <- read.table("example/merge.heart.counts.txt",row.names = 1,header = T)
raw_count <- raw_count[, metasheet[metasheet$qc_passed == TRUE,]$Sample]

## read in gene pos file; keep autosomal genes only
genepos <- read.table("example/UCSC.hg38.gene.pos",header = F)
head(genepos)

genepos <- genepos [genepos $V2 %in% paste0("chr",seq(1,22,1)), ] 
genepos  <- genepos %>%
  as.data.table()
genepos[, V2 := gsub("chr", "", V2)]

raw_count <- raw_count[rownames(raw_count) %in% genepos$V1, ]
raw_count_cpm <- edgeR::cpm(raw_count, normalized.lib.sizes = TRUE, log = FALSE)
dim(raw_count_cpm)

## expression >1 in at least three cells
raw_count_cpm <- raw_count_cpm[which(rowSums(raw_count_cpm>1)>=3), ]
dim(raw_count_cpm)
## gene number over 2000
raw_count_cpm <- raw_count_cpm[,which(colSums(raw_count_cpm>0)>=2000) ]
dim(raw_count_cpm)  #

# Apply gene filter used in 2020GR analysis.
gene_filter <- apply(raw_count_cpm, 1, median) > 50
filtered_cpm <- raw_count_cpm[gene_filter, ] 
dim(filtered_cpm) 

filtered_counts <- raw_count[rownames(raw_count) %in% rownames(filtered_cpm), 
                             colnames(raw_count) %in% colnames(filtered_cpm)]
dim(filtered_counts) 
filtered_counts[1:2,1:3]

filtered_counts <- filtered_counts[, colnames(filtered_counts) %in% metasheet$Sample]
filtered_counts <- filtered_counts[order(rownames(filtered_counts)), order(colnames(filtered_counts))]

dim(filtered_counts)

genepos <- genepos[genepos$V1 %in% rownames(filtered_counts),]
genepos <- genepos[order(genepos$V1),] 
dim(genepos)

## create scploid object
ploidytest_dt <- makeAneu(counts = as.matrix(filtered_counts),
                          genes = rownames(filtered_counts),
                          chrs = genepos$V2,
                          cellNames = metasheet[metasheet$Sample %in% colnames(filtered_counts)]$Sample,
                          cellGroups = metasheet[metasheet$Sample %in% colnames(filtered_counts)]$Group) # split data by EStage and cell type
spt <- splitCellsByGroup(ploidytest_dt) 

# run scploid
expression_results <- do.call(rbind, lapply(spt, calcAneu)) %>%
  as.data.table()
expression_results[1:6,]
# add cell type information
metadata <- metasheet[, c("Sample", "Embryo", "Estage", "Stage", "lineage","cell")] %>%
  setnames(., c("cell", "embryo", "EStage", "stage", "lineage","sample"))
expression_results <- merge(expression_results, metadata, "cell")
setnames(expression_results, "z", "scploid_z")
setnames(expression_results, "score", "scploid_score")
setnames(expression_results, "p", "scploid_p")
expression_results[, cell := gsub("_", ".", cell)]
expression_results[, cell := gsub("-", ".", cell)]

## add ASE data  methods from 2020 GR
if (file.exists(here("example/ase_by_chr.txt"))) {
  ase_by_chr <- fread(here("example/ase_by_chr.txt"))
} else {
  file_list <- list.files(here("example/heart.table"), pattern = "*.table", full.names = TRUE)
  read_ase <- function(file_name) {
    id <- sub('\\.table*$', '', basename(file_name))
    dt <- fread(file_name)
    dt[, cell := id]
    return(dt)
  }
  ase <- do.call(rbind, lapply(file_list, function(x) read_ase(x)))
  ase[, snp_id := paste(contig, position, sep = "_")]
  ase[, embryo := sub("^(.*)[.].*", '\\1', cell)]
  ase[, embryo := gsub("_", ".", embryo)]
  ase[, cell := gsub("_", ".", cell)]
  ase[, cell := gsub("-", ".", cell)]
  ase[, embryo_snp_id := paste(embryo, snp_id, sep = "_")]
  ase[, minCount := pmin(refCount, altCount)]
  ase[, maxCount := pmax(refCount, altCount)]
  
  # summarize ASE per cell-chromosome
  ase_by_chr <- group_by(ase, cell, contig) %>%
    summarize(., allelic_ratio = sum(minCount) / sum(totalCount), 
              min_count_sum = sum(minCount), 
              total_reads = sum(totalCount)) %>%
    as.data.table()
  fwrite(ase_by_chr, file = here("example/ase_by_chr.txt"), 
         sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# merge with expression data
ase_by_chr[1:6,]
ase_by_chr[, chrom := paste(cell, contig, sep = ".")]
expression_results[, sample := gsub("-", ".", sample)]
expression_results[, sample := gsub("_", ".", sample)]
expression_results[, chrom := paste(sample, paste0("chr",chr), sep = ".")]
results <- merge(expression_results, 
                 ase_by_chr[, c("chrom", "allelic_ratio", "total_reads")], 
                 "chrom", all.x = TRUE)
head(results)

# get mapped reads metadata
coverage_metadata <- metasheet[, c("cell", "Mapped Reads")] %>%
  setnames(., c("sample", "mapped_reads"))
coverage_metadata[, sample := gsub("_", ".", sample)]
coverage_metadata[, sample := gsub("-", ".", sample)]
head(coverage_metadata)
results <- merge(results, coverage_metadata, "sample")

# if no SNPs discovered for a chromosome (despite sufficient coverage), set allelic ratio to 0
results[is.na(allelic_ratio), allelic_ratio := 0]
results[is.na(allelic_ratio), total_reads := 0]

# correct allelic imbalance based on coverage
results[, resid_allelic_ratio := resid(lm(data = results, formula = allelic_ratio ~ mapped_reads))]

# estimate variance and compute allelic imbalance z-scores
ase_iqr <- quantile(results$allelic_ratio, c(0.25, 0.75))
iqr_indices <- which(results$allelic_ratio > ase_iqr[1] & results$allelic_ratio < ase_iqr[2])
m_iqr <- results$resid_allelic_ratio[iqr_indices]
null_allelic_ratio <- mean(m_iqr)
null_var_allelic_ratio <- {IQR(results$resid_allelic_ratio)/(2 * qnorm(.75))} ^ 2
results[, ase_z := (resid_allelic_ratio - null_allelic_ratio) / sqrt(null_var_allelic_ratio)]
results[, ase_p := pnorm(ase_z)]
results[, scploid_effect_p := scploid_p]

results[,state := 2]
results[results$scploid_z < -3.65 &  results$ase_z < -1.6, state := 1]
results[results$scploid_z > 5.6 &  results$ase_z > 0.05, state := 3]
results[results$total_reads < 100, state := NA]
head(results)



tmp4 <- results[results$state != 2,.N,by=c("embryo","sample","state","EStage")] 
tmp4
dim(tmp4)

tmp5 <- results[,c("embryo","sample")] 
tmp5
index <- duplicated(tmp5$sample)
tmp5 <- tmp5[!index,] 
dim(tmp5)
tmp6 <- tmp5[,.N,by=c("embryo")]
tmp6

### abnormal number of cells
index2 <- duplicated(tmp4$sample)
tmp7 <- tmp4[!index2,] 
dim(tmp7)
tmp7 <- tmp7[,.N,by=c("embryo")]
tmp7

tmp8 <- merge(tmp6,tmp7,by="embryo")
tmp8[,ratio := round(N.y/N.x,3)]
tmp8
colnames(tmp8) <- c("embryo","all.cell","aneu.cell","aneu.ratio")

################ nature color with needed numbers
adaptive_pal <- function(values) {
  force(values)
  function(n = 10) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

pal_npg_adaptive <- function(palette = c("nrc"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"npg"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}

scale_color_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "npg", pal_npg_adaptive(palette, alpha), ...)
}

scale_fill_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "npg", pal_npg_adaptive(palette, alpha), ...)
}

tmp8$euratio <- 1 - tmp8$aneu.ratio
HE7W.1 <- data.frame(embryo="HE7W.1",all.cell=147,aneu.cell=0,aneu.ratio=0,euratio=1)
tmp8 <- rbind(tmp8,HE7W.1)
short.tmp8 <- melt(tmp8[,c("embryo","aneu.ratio","euratio")])
#short.tmp8$value <- short.tmp8$value*100
name.order <- c("HE5W.1","HE6W","HE6W.2","HE7W.1","HE7W.2","HE7W.3","HE9W.1","HE10W","HE13W.1","HE13W.2",
                "HE15W","HE17W.1","HE17W.2","HE20W.1","HE22W.1","HE23W.2","HE24W","HE25W")
short.tmp8$embryo = factor(short.tmp8$embryo, levels=name.order) ## 设置柱条的顺序

pdf("example/aneu.ratio.in.heart.per.embryo.pdf", width=6,height=4)
ggplot(short.tmp8,aes(x=embryo,y=value,fill=variable))+
  geom_bar(stat = 'identity',position = "stack",width = 0.9)+
  theme_bw() +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  ggtitle('Aneu ratio in different embryos of brain')+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size=12,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 1,size=12,colour = "black"),
        panel.grid=element_blank())+
  scale_x_discrete(limits=name.order)+
  scale_fill_npg_adaptive()+
  ylab("Percentage (%)")
dev.off()


#### combine embryos to sample week
aneu.byweek <- data.frame(week=c("5W","6W","7W","9W","10W","13W","15W","17W","20W","22W","23W","24W","25W"),
                          mean=c(mean(tmp8[grep("5W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("6W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("7W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("9W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("10W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("13W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("15W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("17W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("20W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("22W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("23W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("24W",tmp8$embryo),]$aneu.ratio),
                                 mean(tmp8[grep("25W",tmp8$embryo),]$aneu.ratio)),
                          sd=c(sd(tmp8[grep("5W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("6W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("7W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("9W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("10W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("13W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("15W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("17W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("20W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("22W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("23W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("24W",tmp8$embryo),]$aneu.ratio),
                               sd(tmp8[grep("25W",tmp8$embryo),]$aneu.ratio)),
                          se= c(sd(tmp8[grep("5W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("5W",tmp8$embryo),])),
                                sd(tmp8[grep("6W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("6W",tmp8$embryo),])),
                                sd(tmp8[grep("7W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("7W",tmp8$embryo),])),
                                sd(tmp8[grep("9W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("9W",tmp8$embryo),])),
                                sd(tmp8[grep("10W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("10W",tmp8$embryo),])),
                                sd(tmp8[grep("13W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("13W",tmp8$embryo),])),
                                sd(tmp8[grep("15W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("15W",tmp8$embryo),])),
                                sd(tmp8[grep("17W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("17W",tmp8$embryo),])),
                                sd(tmp8[grep("20W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("20W",tmp8$embryo),])),
                                sd(tmp8[grep("22W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("22W",tmp8$embryo),])),
                                sd(tmp8[grep("23W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("23W",tmp8$embryo),])),
                                sd(tmp8[grep("24W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("24W",tmp8$embryo),])),
                                sd(tmp8[grep("25W",tmp8$embryo),]$aneu.ratio)/sqrt(nrow(tmp8[grep("25W",tmp8$embryo),]))),
                          median=c(median(tmp8[grep("5W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("6W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("7W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("9W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("10W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("13W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("15W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("17W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("20W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("22W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("23W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("24W",tmp8$embryo),]$aneu.ratio),
                                   median(tmp8[grep("25W",tmp8$embryo),]$aneu.ratio)))

name.order <- c("5W","6W","7W","9W","10W","13W","15W","17W","20W","22W","23W","24W","25W")
aneu.byweek$week = factor(aneu.byweek$week, levels=c("5W","6W","7W","9W","10W","13W","15W","17W","20W","22W","23W","24W","25W")) 

pdf("example/aneu.ratio.in.heart.combineweek.mean.pdf")
ggplot(aneu.byweek)+
  geom_bar(aes(x=week,y=mean,fill=week),stat = 'identity')+
  geom_errorbar(aes(x = week ,ymin = mean - se ,ymax = mean + se), width=0.2, colour="black",position = position_dodge(0.9))+
  theme_bw() +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  ggtitle('Aneu ratio in different embryos of post implantation combine weeks')+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size=12,colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1,size=12,colour = "black"),
        panel.grid=element_blank())+
  scale_x_discrete(limits=name.order)+
  scale_fill_npg_adaptive()
dev.off()