####### cell type identification in kidney ###########
rm(list=ls())
.libPaths(c("/mnt/data/kongsiming/software/Rlib" , .libPaths()))
palette <- c("#4DBBD5FF","#E64B35FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF",
                  "#A23B2CFF","#815D45FF","#B09C85FF","#B1716CFF","#2FB0B7FF","#7C98A3FF","#197F87FF","#BB9599FF")


show_col(palette)


library(Seurat)
library(dplyr)
library(patchwork)

count <- read.table("postimplantaion.raw.counts.txt",header = T,sep = "\t") # each row is a gene and each column is a single cell
dim(count)
head(count[1:6,1:6])
temp <- sweep(count,2,colSums(count),"/")
tpm <- temp*1000000
mytpm <- log2(tpm +1)
head(mytpm[1:6,1:6])

colnames(mytpm)<-gsub("X.mnt.data.kongsiming.clinical_ZF.mosaic.post.implantation.05.star.map.map2.","",colnames(mytpm))
colnames(mytpm)<-gsub("Aligned.sortedByCoord.out.bam","",colnames(mytpm))

#rownames(count) <- count[,1]
#count <- count[,-1]
# Initialize the Seurat object with the raw (non-normalized data).

data <- CreateSeuratObject(counts =mytpm,project = "post.implantation",min.cells =3, min.features = 2000)
str(data)
data@assays$RNA@counts[c("A2MP1","A4GALT"),1:30]

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["percent.ercc"]] <- PercentageFeatureSet(data, pattern = "^ERCC-")
# Show QC metrics for the first 5 cells
head(data@meta.data, 5)
head(PercentageFeatureSet(data, pattern = "^MT-"))

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
## calculate the pearson correlation
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)
top100 <- head(VariableFeatures(data), 100)

totalHVG <- VariableFeatures(data)
totalHVG
#write.table(as.data.frame(totalHVG),file="RNA/20201201_results/totalHVG.txt",sep="\t",quote = F,row.names = F)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top100, repel = TRUE)
plot2
CombinePlots(plots = list(plot1, plot2),legend="bottom")

# Scale the data (needed for PCA etc)
data <- ScaleData(object = data)
all.genes <- rownames(data)
###################### Perform linear dimensional reduction ######################################
data <- RunPCA(data, features = VariableFeatures(object = data)) ## using HVG
# Examine and visualize PCA results a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
data <- JackStraw(data, num.replicate = 1000)
data <- ScoreJackStraw(data, dims = 1:20)

JackStrawPlot(data, dims = 1:15)
ElbowPlot(data)
###################### Cluster the cells ####################################################
set.seed(123)
data <- FindNeighbors(data, dims = 1:12)
##my tpm
data <- FindClusters(data, resolution = 0.6)
# Look at cluster IDs of the first 5 cells
head(Idents(data), 5)
###################### Run non-linear dimensional reduction (UMAP/tSNE) #####################
data <- RunUMAP(data, dims = 1:12)
data <- RunTSNE(data, dims = 1:12)

# individual clusters
DimPlot(data, reduction = "umap",label = TRUE, cols = palette)
DimPlot(data, reduction = "tsne",label = TRUE, cols = palette)

#mytpm
genes_EPI <- c("NANOG","POU5F1","SOX2","KLF4","PRMD14","FOXD3","GDF3","NR5A2",
               "UTF1","CLDN19","JMJD2B","JMJD1A","IFITM2")
genes_PE <- c("GATA4","FGFR4","TDGF1","KRT8","KRT19","CLDN3","IFITM1",
              "DPPA5","JMJD4","NODAL")
genes_TE <- c("TFAP2A","ACKR2","CDX2","GATA3","GATA2","KLF5","KRT18",
              "CLDN10","CLDN4","CD46","CD164","ATP1A1","ATP1B3","TET2")
genes_ysTE <- c("CD44")




genes_all <- c(genes_EPI,genes_PE,genes_TE,genes_ysTE)
#DotPlot(data,features = genes_all)
pdf("/mnt/data/kongsiming/clinical_ZF/mosaic/post.implantation/plot/figures/cell.type.pdf",width = 20,height = 40)
VlnPlot(data,features = genes_all)
FeaturePlot(data,genes_all,pt.size = 1,reduction = "tsne")
FeaturePlot(data,genes_all,pt.size = 1,reduction = "umap")
dev.off()

pdf("/mnt/data/kongsiming/clinical_ZF/mosaic/post.implantation/plot/figures/cell.type.dotplot.pdf",width = 20,height = 10)
DotPlot(object = data, features = genes_all)
dev.off()



## annotation 

data@meta.data$free.annotation <- plyr::mapvalues(from = c(9,13,0,1,2,3,4,5,6,7,8,10,11,12,14,15),
                                                  to = c("EPI","PE","TE","TE","TE","TE","TE","TE","TE","TE","TE","TE","TE","TE","TE","ysTE"),
                                                  x= data@meta.data$seurat_clusters)
pdf("/mnt/data/kongsiming/clinical_ZF/mosaic/post.implantation/plot/figures/post.implantation.cell.type.cluster.2.15.pdf")
TSNEPlot(object = data,group.by= "free.annotation",label=T, cols = palette)    
UMAPPlot(object = data,group.by= "free.annotation",label=T, cols = palette)   
dev.off()
# saveRDS(data, file = "/mnt/data/kongsiming/clinical_ZF/mosaic/post.implantation/plot/figures/data.post.implantation.2.15.rds")
###################### extract cluster ###############################
table(data@active.ident) 
table(Idents(data))
table(data$RNA_snn_res.0.6)

# extract cell type
tmp2 <- as.data.frame(data@active.ident)
write.table(data@meta.data,file = "/mnt/data/kongsiming/clinical_ZF/mosaic/post.implantation/plot/figures/cell.cluster.final.2.15.txt",quote = F,sep = '\t',row.names = T,col.names = T)


