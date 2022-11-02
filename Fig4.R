## Aneuploid detection
#.libPaths("E:/Program Files (x86)/R/R-4.1.1/library")
rm(list = ls())
library(dplyr)
library(scran)
library(tidyr)
library(readxl)
library(scploid)
library(edgeR)
library(data.table)
library(tidyverse)
library(broom)
library(survcomp)
library(TreeBH)
library(ggrepel)
library(lme4)
library(margins)
library(ggdendro)
library(scales)
library(cowplot)
library(here)

here <- here::here
summarize <- dplyr::summarize
## Read in metadata sheet
#metasheet <- read_xlsx("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/renal.metadata.xlsx", sheet = 1)
metasheet <- read_xlsx("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/brain.metadata.all.regions.xlsx", sheet = 1)
head(metasheet)

metasheet <- metasheet %>%
  as.data.table()
dim(metasheet)
metasheet <- metasheet[order(metasheet$Sample),]
metasheet$Group <- paste0(metasheet$Estage, "_", metasheet$lineage)
dim(metasheet)
head(metasheet)
#View(metasheet)
metasheet[, qc_passed := (`Percent Mapped` > quantile(`Percent Mapped`, 0.1) & `Mapped Reads` > quantile(`Mapped Reads`, 0.1))]
metasheet[, Embryo := gsub("_", ".", Embryo)]

raw_count <- read.table("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/merge.brain.19embryos.counts.txt",row.names = 1,header = T)
raw_count[1:5,1:2]
dim(raw_count)  ## 26485 genes and 3600 cells renal   26485 genes and 2607 cells heart    26485 genes and 4124 cells brain 6embryos 26485 genes and 6044 cells brain in 10 embryos
raw_count <- raw_count[, metasheet[metasheet$qc_passed == TRUE,]$Sample]
dim(raw_count)  ## 26485 genes and 3508 cells brain

## read in gene pos file; keep autosomal genes only
genepos <- read.table("D:/QiaoLab/multi-omics-zhaifan/clinical.type/UCSC.hg38.gene.pos",header = F)
dim(genepos)
head(genepos)
genepos <- genepos [genepos $V2 %in% paste0("chr",seq(1,22,1)), ] # 25225 genes
dim(genepos)  # 25225 genes
genepos  <- genepos %>%
  as.data.table()

genepos[, V2 := gsub("chr", "", V2)]
raw_count <- raw_count[rownames(raw_count) %in% genepos$V1, ]
dim(raw_count)   
## kong: calculate cpm using edgeR package
## raw counts ???��� ??edgeR ת????normalize ֮?󣬵??ǲ?ȡlog??֮???Լ?ȡlog2(emtab3929_cpm+1)
raw_count_cpm <- edgeR::cpm(raw_count, normalized.lib.sizes = TRUE, log = FALSE)
dim(raw_count_cpm)
## expression >1 in at least three cells
raw_count_cpm <- raw_count_cpm[which(rowSums(raw_count_cpm>1)>=3), ]
dim(raw_count_cpm)
## gene number over 2000
raw_count_cpm <- raw_count_cpm[,which(colSums(raw_count_cpm>0)>=2000) ]
dim(raw_count_cpm)  # 21618 genes 3471 cells
## log2 transform
raw_count_log2cpm <- log2(raw_count_cpm + 1)
dim(raw_count_log2cpm)
# Apply gene filter used in 2020GR analysis.
gene_filter <- apply(raw_count_cpm, 1, median) > 50
filtered_cpm <- raw_count_cpm[gene_filter, ] 
dim(filtered_cpm) # 815 genes and 3048 cells renal 794genes and 2002 cells heart 1393genes and 3471 cells brain in 6 embryos and 5017 cells in 10 embryos in brain and 10074 cells in 19 embryos
filtered_log2cpm <- log2(filtered_cpm + 1)
dim(filtered_log2cpm)
filtered_counts <- raw_count[rownames(raw_count) %in% rownames(filtered_cpm), 
                                    colnames(raw_count) %in% colnames(filtered_cpm)]
dim(filtered_counts) #1393 genes and 3471 cells
filtered_counts[1:2,1:3]

filtered_counts <- filtered_counts[, colnames(filtered_counts) %in% metasheet$Sample]
filtered_counts <- filtered_counts[order(rownames(filtered_counts)), order(colnames(filtered_counts))]
#View(filtered_counts)
dim(filtered_counts)  ## 799 3257 renal 794 2002 heart  1393 genes and 4208 cells brain
filtered_counts[1:6,1:6]
genepos <- genepos[genepos$V1 %in% rownames(filtered_counts),]
dim(genepos)
genepos <- genepos[order(genepos$V1),] 
 
head(genepos)
dim(genepos)

## create scploid object
## chrs ??Ҫ???ݵ?һ???Լ?????
ploidytest_dt <- makeAneu(counts = as.matrix(filtered_counts),
                          genes = rownames(filtered_counts),
                          chrs = genepos$V2,
                          cellNames = metasheet[metasheet$Sample %in% colnames(filtered_counts)]$Sample,
                          cellGroups = metasheet[metasheet$Sample %in% colnames(filtered_counts)]$Group) # split data by EStage and cell type

# data split into 13 subsets, one for each EStage_celltype combination
spt <- splitCellsByGroup(ploidytest_dt) 

# run scploid
expression_results <- do.call(rbind, lapply(spt, calcAneu)) %>%
  as.data.table()

# add cell type information from Stirparo et al.
metasheet <- metasheet %>% 
  as.data.table()

metadata <- metasheet[, c("Sample", "Embryo", "Estage", "Stage", "lineage","cell")] %>%
  setnames(., c("cell", "embryo", "EStage", "stage", "lineage","sample"))
expression_results <- merge(expression_results, metadata, "cell")
setnames(expression_results, "z", "scploid_z")
setnames(expression_results, "score", "scploid_score")
setnames(expression_results, "p", "scploid_p")
expression_results[, cell := gsub("_", ".", cell)]
expression_results[, cell := gsub("-", ".", cell)]
#expression_results <- expression_results[order(expression_results[,8],decreasing=F),] 
## add ASE data
if (file.exists(here("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/ase_by_chr.txt"))) {
  ase_by_chr <- fread(here("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/ase_by_chr.txt"))
} else {
  file_list <- list.files(here("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/brain.table"), pattern = "*.table", full.names = TRUE)
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
  fwrite(ase_by_chr, file = here("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/ase_by_chr.8embryo.txt"), 
         sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# merge with expression data
ase_by_chr[, chrom := paste(cell, contig, sep = ".")]
expression_results[, sample := gsub("-", ".", sample)]
expression_results[, sample := gsub("_", ".", sample)]
expression_results[, chrom := paste(sample, paste0("chr",chr), sep = ".")]
results <- merge(expression_results, 
                 ase_by_chr[, c("chrom", "allelic_ratio", "total_reads")], 
                 "chrom", all.x = TRUE)
dim(results)
head(results)
##########################################test here 2021.8.6 ####################################################
# get mapped reads metadata
coverage_metadata <- metasheet[, c("cell", "Mapped Reads")] %>%
  setnames(., c("sample", "mapped_reads"))
coverage_metadata[, sample := gsub("_", ".", sample)]
coverage_metadata[, sample := gsub("-", ".", sample)]
head(coverage_metadata)
results <- merge(results, coverage_metadata, "sample")
dim(results)
results[1:6,]
# View(results)
# if no SNPs discovered for a chromosome (despite sufficient coverage), set allelic ratio to 0
results[is.na(allelic_ratio), allelic_ratio := 0]
results[is.na(allelic_ratio), total_reads := 0]

# correct allelic imbalance based on converage
results[, resid_allelic_ratio := resid(lm(data = results, formula = allelic_ratio ~ mapped_reads))]

# estimate variance and compute allelic imbalance z-scores
ase_iqr <- quantile(results$allelic_ratio, c(0.25, 0.75))
iqr_indices <- which(results$allelic_ratio > ase_iqr[1] & results$allelic_ratio < ase_iqr[2])
m_iqr <- results$resid_allelic_ratio[iqr_indices]
null_allelic_ratio <- mean(m_iqr)
null_var_allelic_ratio <- {IQR(results$resid_allelic_ratio)/(2 * qnorm(.75))} ^ 2
results[, ase_z := (resid_allelic_ratio - null_allelic_ratio) / sqrt(null_var_allelic_ratio)]
results[, ase_p := pnorm(ase_z)]
#write.table(results,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/brain.results.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)


cor.test(results$scploid_p, results$ase_p, method = "kendall")

fisher_wrapper <- function(pval_1, pval_2, wt_1 = 1, wt_2 = 1) {
  return(tryCatch(combine.test(p = c(pval_1, pval_2), weight = c(wt_1, wt_2), method = "fisher"), error = function(e) NA))
}



# impose effect size threshold, consistent with Griffiths et al.; set p-values to 1
results[, scploid_effect_p := scploid_p]
results[(scploid_score > 0.8 & scploid_score < 1.2), scploid_effect_p := 1]
results[, fisher_p := mapply(fisher_wrapper, results$ase_p, results$scploid_effect_p)]


# set FDR = 1%
fdr <- 0.01
sc_groups <- as.matrix(results[, c("embryo", "sample", "chrom")])
calls <- suppressWarnings(get_TreeBH_selections(results$fisher_p,
                                                sc_groups,
                                                q = c(fdr, fdr, fdr)))
head(calls)

results[, sig_embryo := calls[, 1]]
results[, sig_cell := calls[, 2]]
results[, sig_chrom := calls[, 3]]

results[cell %in% unique(results[sig_cell == 1]$cell), sig_cell := 1]
results[embryo %in% unique(results[sig_embryo == 1]$embryo), sig_embryo := 1]

length(unique(results[sig_chrom == 1]$chrom))
length(unique(results[sig_chrom == 1]$cell))
length(unique(results[sig_chrom == 1]$embryo))

cell_fraction <- group_by(results[!duplicated(cell)], EStage, embryo) %>%
  summarize(., prop_aneuploid = mean(sig_cell), n = n()) %>%
  as.data.table()

cell_fraction_plot <- ggplot(data = cell_fraction) +
  geom_histogram(aes(x = prop_aneuploid), bins = 30) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Prop. aneuploid cells") +
  ylab("Number of embryos")

#plot_grid(fdr_plot, cell_fraction_plot, labels = c('A', 'B'))

chr_fraction <- group_by(results, EStage, embryo, chr) %>%
  summarize(., prop_aneuploid = mean(sig_chrom), n = n()) %>%
  as.data.table()


# calculate aneuploidies per chromosome
aneuploid_by_chr <- group_by(results, chr) %>%
  summarize(aneuploid_cells = sum(sig_chrom == 1), euploid_cells = sum(sig_chrom == 0)) %>%
  as.data.table()
aneuploid_by_chr[, chr:= paste0("chr",chr)]
#aneuploid_by_chr[, chrom:= NULL]
#aneuploid_by_chr$chr <- factor(aneuploid_by_chr$chr, paste0("chr",c(1:22)))
aneuploid_by_chr[, chr_numeric := gsub("chr", "", chr)]
aneuploid_by_chr$chr_numeric <- factor(aneuploid_by_chr$chr_numeric, c(1:22))

gencode <- fread("F:/GR/gencode.v32.annotation.gtf.gz")
genes_per_chr <- data.table(table(gencode[V3 == "gene" & grepl("protein_coding", V9)][, 1])) %>%
  setnames(., c("chr", "n_genes"))
chrom_lengths <- fread("F:/GR/hg38.chrom.sizes.txt") %>%
  setnames(., c("chr", "len"))

aneuploid_by_chr <- merge(chrom_lengths, merge(aneuploid_by_chr, genes_per_chr, "chr"), "chr")

cor.test(aneuploid_by_chr$aneuploid_cells, aneuploid_by_chr$len)
cor.test(aneuploid_by_chr$aneuploid_cells, aneuploid_by_chr$n_genes)

by_chrom_plot <- ggplot(data = aneuploid_by_chr, aes(x = n_genes , y = aneuploid_cells, label = chr)) +
  theme_classic() +
  xlab("Number of protein-coding genes") +
  ylab("Number of aneuploid cells") +
  geom_point() +
  geom_label_repel(size = 4) +
  ylim(0, 125) +
  xlim(0, 2250)

m1 <- glmer(data = results, formula = (sig_chrom == 1) ~ (1 | embryo / cell) + (1 | lineage) + chr, family = binomial, nAGQ = 0)
m0 <- glmer(data = results, formula = (sig_chrom == 1) ~ (1 | embryo / cell) + (1 | lineage), family = binomial, nAGQ = 0)
anova(m1, m0, test = "Chisq") 

# embryo-specific models; to average over levels of random effect, need to extract average marginal effects as below
mx <- margins(m1, type = "response", variables = "chr")
b <- summary(mx)
cov_mat <- attr(mx, "vcov")
k <- diag(nrow = ncol(cov_mat))
kvar <- t(k) %*% cov_mat %*% k
kb <- k %*% b$AME
my_chi <- t(kb) %*% solve(kvar) %*% kb
pchisq(my_chi[1, 1], df = ncol(cov_mat), lower.tail = F)

# cluster with k-means, assign clusters to monosomy and trisomy
km <- results[sig_chrom == 1, c("scploid_z", "ase_z")] %>%
  kmeans(centers = 2)
km_clusters <- results[sig_chrom == 1, c("chrom", "scploid_z", "ase_z")]
km_clusters[, cluster := km$cluster]

results <- merge(results, km_clusters[, c("chrom", "cluster")], "chrom", all.x = TRUE)

results[, ploidy := 2]
if (mean(results[cluster == 1]$ase_z) < mean(results[cluster == 2]$ase_z)) {
  results[cluster == 1, ploidy := 1]
  results[cluster == 2, ploidy := 3] 
} else {
  results[cluster == 2, ploidy := 1]
  results[cluster == 1, ploidy := 3] 
}

########################
## find threshold using log_fisher_p
results[monosomy == TRUE, log_fisher_p := -1 * log_fisher_p]
tmp2 <-results


tmp_triploid <- results[results$log_fisher_p >2,]
tmp_triploid
hist(tmp_triploid$log_fisher_p)
#View(tmp_triploid)
dim(tmp_triploid)
tmp_triploid[tmp_triploid$log_fisher_p > 5,]$log_fisher_p <- 5


hist(tmp_triploid$log_fisher_p,breaks = 100,col = "gray",main = "threshold for triploid in brain",cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
rug(jitter(tmp_triploid$log_fisher_p),side = 1,col = 2)
plot(density(tmp_triploid$log_fisher_p),type = "l",main="threshold for triploid in brain",xlab="range",ylab="density")
rug(jitter(tmp_triploid$log_fisher_p),side = 1,col = 2)
## find local minimum 
d <- density(tmp_triploid$log_fisher_p)
v <- optimize(approxfun(d$x,d$y),interval=c(2,5))$minimum
v
abline(v=v, col="blue")  ## 4.32  10embryos 4.40 19embryos 4.45


tmp_haploid <- results[results$log_fisher_p < -2,]
tmp_haploid

tmp_haploid[tmp_haploid$log_fisher_p < -5,]$log_fisher_p <- -5
#View(tmp_haploid)

hist(tmp_haploid$log_fisher_p,breaks = 100,col = "gray",main = "threshold for haploid in brain",cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
rug(jitter(tmp_haploid$log_fisher_p),side = 1,col = 2)
plot(density(tmp_haploid$log_fisher_p),type = "l",main="threshold for haploid in brain",xlab="range",ylab="density")
rug(jitter(tmp_haploid$log_fisher_p),side = 1,col = 2)
## find local minimum 
d <- density(tmp_haploid$log_fisher_p)
v <- optimize(approxfun(d$x,d$y),interval=c(-5,-3))$minimum
v
abline(v=v, col="blue")  ## -4.24  10brain embryos: -4.46  19brain embryos:-4.63
results[,state := 2]
results[results$log_fisher_p < -4.49, state := 1]
results[results$log_fisher_p > 4.42, state := 3]
write.table(results,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/brain.results.08.20.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

### kong: this is exactly what i want

plot_mca <- function(embryo_id, combine = FALSE, legend = TRUE, cluster_method = "ward.D2", dist_method = "euclidean") {
  heatmap_data <- pivot_wider(results[(embryo == embryo_id), c("chr_num", "cell", "ploidy")], names_from = chr_num, values_from = ploidy)
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$cell
  
  distance.row <- dist(heatmap_matrix, method = dist_method)
  cluster.row <- hclust(distance.row, method = cluster_method)
  
  dendrogram <- ggplot(segment(dendro_data(cluster.row))) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0.2, 0)) + 
    theme_dendro() + 
    scale_x_reverse()
  
  sample_order <- rev(cluster.row$labels[cluster.row$order])
  
  dt_to_plot <- results[embryo == embryo_id]
  dt_to_plot$cell <- factor(dt_to_plot$cell, levels = sample_order)
  
  exp_heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = scploid_z)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradientn(name = "Z-score", colors = c("blue", rep("white", 3), "red"), limits = c(-5, 5), na.value = 1, oob = squish) +
    xlab("Chromosome") +
    ylab("Cell") +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  
  ase_heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = ase_z)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradientn(name = "Z-score", colors = c("blue", rep("white", 3), "red"), limits = c(-3, 3), na.value = 1, oob = squish) +
    xlab("Chromosome") +
    ylab("Cell") +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  
  #dt_to_plot[monosomy == TRUE, log_fisher_p := -1 * log_fisher_p]
  
  # heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell, fill = log_fisher_p)) +
  # geom_tile() +
  # theme_bw() +
  # scale_fill_gradientn(name = "-log10(p)", colors = c("blue", "white", "white", "red"), limits = c(-5, 5), oob = squish) +
  # xlab("Chromosome") +
  # ylab("") +
  # theme(axis.text.y = element_blank(), plot.margin = unit(c(5.5, 5.5, 5.5, -3), "pt"), panel.grid = element_blank())
  dt_to_plot$state <- as.factor(dt_to_plot$state)
  heatmap <- ggplot(data = dt_to_plot, aes(x = chr_num, y = cell)) +
    geom_tile(aes(fill = state)) +
    theme_bw() +
   # scale_fill_manual(name = "Copy number", values = c( "white","#BB0021FF"))+ ## 2,3
    scale_fill_manual(name = "Copy number", values = c("#008280FF", "white","#BB0021FF")) +  ## green red
    #scale_fill_manual(name = "Copy number", values = c("#0342E5","white","#E40000")) +   
    ## blue red
    xlab("Chromosome") +
    ylab("") +
    theme(axis.text.y = element_blank(), plot.margin = unit(c(5.5, 5.5, 5.5, -3), "pt"), panel.grid = element_blank())
  
  
  if (combine == TRUE) {
    plot_grid(dendrogram, heatmap, align = "h", axis = "b", rel_widths = c(0.3, 1), scale = c(1, 0.95))
  } else {
    plot_grid(exp_heatmap, ase_heatmap, align = "h", axis = "b", rel_widths = c(1, 1))
  }
}


for (embryo_id in unique(results$embryo)) {
  pdf(file = paste0(here("D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain"), embryo_id, ".20211026.pdf"), height = 4, width = 6)
  try(print(plot_mca(embryo_id, combine = TRUE, cluster_method = "ward.D2", dist_method = "euclidean")))
  dev.off()
}
pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/原版.all.embryos.pdf")
plot_mca("10W3", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("12W2", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("14W2", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("16W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("7W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("9W2", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("16W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("21W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
plot_mca("26W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()

#plot_grid(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'))
pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE10W3.pdf")
plot_mca("10W3", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE12W2.pdf")
plot_mca("12W2", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE14W2.pdf")
plot_mca("14W2", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE16W.pdf")
plot_mca("16W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE7W.pdf")
plot_mca("7W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE9W2.pdf")
plot_mca("9W2", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()


pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE26W.pdf")
plot_mca("26W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()


## ONLY chromosome gain
pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/原版.HE21W.pdf")
plot_mca("21W", combine = TRUE, cluster_method = "average", dist_method = "euclidean")
dev.off()


tmp4 <- results[results$state != 2,.N,by=c("embryo","sample","state","EStage")] 
tmp4
dim(tmp4)
#View(tmp4)
#write.table(tmp4,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/CNV.rna.code/ASE_counts/brain/brain.state.num.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)


tmp5 <- results[,c("embryo","sample")] 
tmp5
index <- duplicated(tmp5$sample)
tmp5 <- tmp5[!index,] 
dim(tmp5)
tmp6 <- tmp5[,.N,by=c("embryo")]
tmp6
#write.table(tmp6,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.11.19/器官表格/brain.cell.num.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

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
write.table(tmp8,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/gastrulation/gas.aneu.ratio.txo.fig5b.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/brain/aneu.ratio.in.brain.pdf")
name.order <- c("7W","9W2","10W3","12W2","14W2","16W","21W","26W")
tmp8$embryo = factor(tmp8$embryo, levels=c("7W","9W2","10W3","12W2","14W2","16W","21W","26W")) ## 设置柱条的顺序
ggplot(tmp8,aes(x=embryo,y=aneu.ratio,fill=embryo))+
  geom_bar(stat = 'identity',fill =c("#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF","#BC3C29FF","#0072B5FF"))+
  theme_bw() +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  ggtitle('Aneu ratio in different embryos of brain')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_discrete(limits=name.order)
dev.off()

#### combine embryos to sample week
tmp5 <- results[,c("EStage","sample")] 
tmp5
index <- duplicated(tmp5$sample)
tmp5 <- tmp5[!index,] 
dim(tmp5)
tmp6 <- tmp5[,.N,by=c("EStage")]
tmp6
#write.table(tmp6,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.11.19/器官表格/brain.cell.num.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

### abnormal number of cells
index2 <- duplicated(tmp4$sample)
tmp7 <- tmp4[!index2,] 
dim(tmp7)
tmp7 <- tmp7[,.N,by=c("EStage")]
tmp7

tmp8 <- merge(tmp6,tmp7,by="EStage")
tmp8[,ratio := round(N.y/N.x,3)]
tmp8
colnames(tmp8) <- c("embryo","all.cell","aneu.cell","aneu.ratio")
write.table(tmp8,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/brain.aneu.ratio.to.EStage.fig5b.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/aneu.ratio.in.brain.combineweek.pdf")
name.order <- c("7W","8W","9W","10W","11W","12W","14W","16W","17W","20W","21W","26W")
ggplot(tmp8,aes(x=embryo,y=aneu.ratio))+
  geom_bar(stat = 'identity',fill ="#DC00007F")+
  theme_bw() +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  ggtitle('Aneu ratio in different embryos of brain combine weeks')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_discrete(limits=name.order)
dev.off()



##################################################
## for all embryos pheatmap
library(pheatmap)
results <- read.table("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/brain.results.12.28.txt",header = T)
head(results)
dim(results)
tmp10 <- results[,c("sample","chr","state")]
head(tmp10)
tmp11 <- dcast(tmp10,sample~tmp10$chr,value.var = "state")
dim(tmp11)
head(tmp11)

HE7W <- tmp11[substr(tmp11$sample, 1, 4) == "HE7W", ]
head(HE7W)
HE8W <- tmp11[substr(tmp11$sample, 1, 4) == "HE8W", ]
head(HE8W)
HE9W <- tmp11[substr(tmp11$sample, 1, 4) == "HE9W", ]
head(HE9W)
HE10 <- tmp11[substr(tmp11$sample, 1, 4) == "HE10", ]
head(HE10)
HE11 <- tmp11[substr(tmp11$sample, 1, 4) == "HE11", ]
head(HE11)
HE12 <- tmp11[substr(tmp11$sample, 1, 4) == "HE12", ]
head(HE12)
HE14 <- tmp11[substr(tmp11$sample, 1, 4) == "HE14", ]
head(HE14)
HE16 <- tmp11[substr(tmp11$sample, 1, 4) == "HE16", ]
head(HE16)
HE17 <- tmp11[substr(tmp11$sample, 1, 4) == "HE17", ]
head(HE17)
HE20 <- tmp11[substr(tmp11$sample, 1, 4) == "HE20", ]
head(HE20)
HE21 <- tmp11[substr(tmp11$sample, 1, 4) == "HE21", ]
head(HE21)
HE26 <- tmp11[substr(tmp11$sample, 1, 4) == "HE26", ]
head(HE26)


newdata <- rbind(HE7W,HE8W,HE9W, HE10,HE11,HE12,HE14,HE16,HE17,HE20,HE21,HE26)
newdata <- data.frame(newdata)
rownames(newdata) <- newdata$sample
newdata <- newdata[,-1]
head(newdata)
tail(newdata)
#newdata <- as.numeric(newdata)
summary(newdata)


library(pheatmap)
annotation_row = data.frame(
  Group = factor(substr(rownames(newdata), 1, 4), levels = c("HE7W", "HE8W","HE9W", "HE10","HE11","HE12","HE14","HE16","HE17","HE20",
                                                             "HE21","HE26")
                                                      
  )
)
table(annotation_row)

#ann_colors = list(
  Embryo = c(HE7W="#BC3C29FF",HE9W="#0072B5FF", HE10="#E18727FF",HE12="#20854EFF",HE14="#7876B1FF",HE16="#6F99ADFF")
)


rownames(annotation_row) <- rownames(newdata)

pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/19brain.pheatmap.green.pdf",width = 6,height = 6)
pheatmap(newdata,scale = "none",annotation_row = annotation_row,annotation_legend = T,legend = T, 
         #color = c("#008280FF",  "white", "#BB0021FF"), 
         #color = c("#0342E5","white","#E40000"), #blue
         color = c("#008280FF", "white", "#BB0021FF"),
        # ann_colors = ann_colors,
         show_rownames=F,show_colnames = T, border_color="black",cluster_rows = F,cluster_cols = F)


dev.off()
########################### this part is for statistics ###############################################
## hist-density
rm(list = ls())
state_num <- read.table("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/brain/brain.results.08.20.txt",header = T)
dim(state_num)
data <- state_num[,c("sample","embryo","chrom","state")]
data <- data[data$state != 2,]
tmp <- data
dim(tmp)
head(tmp)
head(data)
tmp <- tmp %>% as.data.table()


tmp[,chr:= do.call(rbind,strsplit(as.character(tmp$chrom),split = "chr"))[,2]] 

tmp_state1 <- tmp[tmp$state ==1,c("sample","chr")]
tmp_state1$value <- 1
tmp_state3 <- tmp[tmp$state ==3,c("sample","chr")]
tmp_state3$value <- 1
state1 <- reshape(tmp_state1, idvar = "sample", timevar = "chr", direction = "wide")
head(state1)
state3 <- reshape(tmp_state3, idvar = "sample", timevar = "chr", direction = "wide")
head(state3)


### abnormal chromosome in each embryo

chr_freq <- tmp  %>%
  .[,.N,by=c("embryo","chr")]
chr_freq
write.table(chr_freq,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/brain/brain.chr.ab.freq.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

#scale_x_discrete(limits=name.order)
## Density plot of error number in different embryos of brain
abnormal <- tmp  %>%
  .[,.N,by=c("embryo","sample")]
abnormal
head(abnormal)
abnormal$embryo = factor(abnormal$embryo, levels=c("7W","9W2","10W3","12W2","14W2","16W","21W","26W"))
pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/brain.statistic.pdf")
ggplot(data=abnormal,mapping = aes(x=N,color=embryo,..scaled..,fill = embryo))+
  geom_density(alpha = 0.2)+
  scale_x_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme_bw() +
  ggtitle('Density plot of error number in different embryos of brain')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
#######################################
## hist density from 0 on
rm(list = ls())
library(ggridges)
state_num <- read.table("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/brain/brain.results.08.20.txt",header = T)
data <- state_num[,c("sample","embryo","chrom","state")]
tmp1 <- data[data$state == 2,]

tmp1 <- tmp1 %>% as.data.table()
dim(tmp1)
index <- duplicated(tmp1$sample)
tmp1 <- tmp1[!index,] 
dim(tmp1)
head(tmp1)
tmp2 <- tmp1[,.N,by=c("embryo","sample")]
'%!in%' <- function(x,y)!('%in%'(x,y))
tmp3 <- tmp2[tmp2$sample %!in% abnormal$sample,]
tmp3$N <- 0
tmp3

newdata <- rbind(abnormal,tmp3)

name.order <- c("26W","21W","16W","14W2","12W2","10W3","9W2","7W")
newdata$EStage = factor(newdata$EStage, levels=rev(name.order))
write.table(newdata,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/brain.burden.ratio.0.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

class(newdata)

burden.embryo <- newdata  %>%
  .[,.N,by=c("embryo","N")]
write.table(burden.embryo,"D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/brain.burden.embryo.ratio.0.txt",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

burden <- newdata  %>%
  .[,.N,by=c("N")]


pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/20220820/brain/new.brain.statistic.from0.on2.pdf",height = 6,width = 8)
ggplot(data=newdata,mapping = aes(x=N,y=embryo,color=embryo,fill = embryo))+
  geom_density_ridges(alpha = 0.8) +
  scale_x_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme_bw() +
  scale_y_discrete(limits=name.order)+
  ggtitle('Density plot of error number in different embryos of brain from 0 on')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF"))


ggplot(data=newdata,mapping = aes(x=N,color=embryo,..scaled..,fill = embryo))+
  geom_density(alpha = 0.2)+
  scale_x_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme_bw() +
  ggtitle('Density plot of error number in different embryos of brain from 0 on')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggplot(data=newdata,mapping = aes(x=N,color=embryo,..scaled..))+
  geom_density(alpha = 0.9,size=0.5) +
  scale_x_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme_bw() +
  ggtitle('Density plot of error number in different embryos of brain from 0 on')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF"))

ggplot(data=newdata,mapping = aes(x=N,color=embryo,..scaled..))+
  geom_density(alpha = 0.9,size=0.5) +
  scale_x_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme_bw() +
  ggtitle('Density plot of error number in different embryos of brain from 0 on')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_manual(values=c("#BC3C29FF"))
dev.off()



### gain/loss ratio in different embryos of brain
tmp$state[which(tmp$state=="3" )] <- "gain"
tmp$state[which(tmp$state=="1" )] <- "loss"

gloss <- tmp %>%
  .[,.N,by=c("embryo","state")]

long <- dcast(gloss, embryo ~ state)  
long[is.na(long)] <- 0
long <- long %>%
  as.data.table() %>%
  .[,gainratio := gain/(gain+loss)] %>%
  .[,lossratio := loss/(gain+loss)]

short <- melt(long[,c(1,4,5)])
name.order <- c('7W','8W','9W2','10W','10W3','11W','12W','12W2','12W3','12W4','14W2','16W',
                '17W2','20W1','20W2','21W','21W2','26W','26W3')
pdf("D:/QiaoLab/multi-omics-zhaifan/clinical.type/文章撰写/2021.12.8 器官更新/brain.19.embryos/19embryos.gain.loss.pdf",height = 6,width = 8)

ggplot(data = short, mapping = aes(x = embryo, y = value, fill = variable)) + 
  geom_bar(stat= 'identity', position = 'stack')+
  scale_fill_manual(values=c("#BB0021FF","#008280FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  ggtitle('gain/loss in brain')+
  scale_x_discrete(limits=name.order)+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

dev.off()
