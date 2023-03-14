library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
rm(list = ls())
#### input data
data_LN <- read.csv("GSE164522_CRLM_LN_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_LN <- CreateSeuratObject(
  data_LN,
  project = "GSE164522", 
  #min.cells = 10,
  #min.features = 200
 )
save(sce_LN,file=paste0("./","sce_LN.RData"))
rm(data_LN)

data_MN <- read.csv("GSE164522_CRLM_MN_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_MN <- CreateSeuratObject(
  data_MN,
  project = "GSE164522", 
  #min.cells = 10,
  #min.features = 200
)
save(sce_MN,file=paste0("./","sce_MN.RData"))
rm(data_MN)

data_MT <- read.csv("GSE164522_CRLM_MT_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_MT <- CreateSeuratObject(
  data_MT,
  project = "GSE164522", 
  #min.cells = 10,
  #min.features = 200
)
save(sce_MT,file=paste0("./","sce_MT.RData"))
rm(data_MT)

data_PBMC <- read.csv("GSE164522_CRLM_PBMC_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_PBMC <- CreateSeuratObject(
  data_PBMC,
  project = "GSE164522", 
  #min.cells = 10,
  #min.features = 200
)
save(sce_PBMC,file=paste0("./","sce_PBMC.RData"))
rm(data_PBMC)

data_PN<- read.csv("GSE164522_CRLM_PN_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_PN <- CreateSeuratObject(
  data_PN,
  project = "GSE164522", 
  #min.cells = 10,
  #min.features = 200
)
save(sce_PN,file=paste0("./","sce_PN.RData"))
rm(data_PN)

data_PT <- read.csv("GSE164522_CRLM_PT_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_PT <- CreateSeuratObject(
  data_PT,
  project = "GSE164522", 
  #min.cells = 10,
  #min.features = 200
)
save(sce_PT,file=paste0("./","sce_PT.RData"))
rm(data_PT)

# merge 合并seurat除LN外的5种组织
sce<- merge(sce_MN, y=c(sce_MT, sce_PBMC, sce_PN, sce_PT))
save(sce,file=paste0("./","sce0_of_5_tissue.RData"))
rm(sce_MN,sce_MT, sce_PBMC, sce_PN, sce_PT)
load(file = 'sce0_of_5_tissue.RData')

##output data file
sam.name <- "test1"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
### QC
#############################################################################################################
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#画图的代码需要选中从pdf一直到dev.off()的所有代码，一起运行
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

#### 5. 筛选细胞 ####
cat("Before filter :",nrow(sce@meta.data),"cells\n")
#sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
#cat("After filter :",nrow(sce@meta.data),"cells\n")

# ## Normalization
# ############################################################################################################
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)  # Checkpoint
# ## Feature selection
# ############################################################################################################
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

#均一化
#这是2000高变基因的归一化
sce <- ScaleData(sce,
                 vars.to.regress = c("percent.mt"))
#############################################################################################################
# ## Reduction
#############################################################################################################
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(sce, reduction = "pca",ndims = 50)
dev.off()

### Cluster
dim.use <- 1:50
# ############################################################################################################
sce <- RunUMAP(sce, dims = dim.use) # umap tsne
sce <- FindNeighbors(sce, dims = 1:50)  # louvain cluster, graph based
sce <- FindClusters(sce, resolution = 0.8)

##seurat_clusters_umap
pdf(paste0("./",sam.name,"/CellCluster-DimPlot_umap_",max(dim.use),".pdf"),width = 10,height = 9)
DimPlot(sce, reduction = "umap",group.by = 'seurat_clusters',label = T)
dev.off()

##major_marker_gene_umap
pdf(paste0("./",sam.name,"/umap.",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","TNFRSF17","COL1A1","VWF","COL1A2","PLVAP","LILRA4","FCGR3B"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

## 气泡图
features <- c("TFF3","AGR2","EPCAM",#Epithelial
              "CD3E","CD3D","CCL5",#T
              "PLVAP","VWF","PECAM1",#
              "JCHAIN","MZB1","IGHA1",#plasma
              "COL1A1","COL1A2","DCN",#Fibroblast
              "LYZ","CD68","C1QC",#Myeloid
              "MS4A1","CD79A","CD19",#B
              "TPSB2","TPSAB1","CPA3"#Mast
)
features= c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","FCGR3B","VWF",
            "LILRA4","CD68","COL1A2","FCGR3B","PLVAP","CD14")

## major_marker_gene 气泡图
pdf(paste0("./",sam.name,"/DotPlot2_umap.",max(dim.use),"PC.pdf"),width = 10,height = 10)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()

#细胞类群相似树
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/CellCluster-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()

##在meta.data中重命名seurat_clusters为celltype
#使用unname函数配合向量：
cluster2celltype <- c( "0"="T",
                       "1"="T", 
                       "2"="T", 
                       "3"= "NK", 
                       "4"= "T", 
                       "5"= "B",
                       "6"= "Plasma", 
                       "7"= "T", 
                       "8"= "T",
                       "9"= "NK",
                       "10"="Myeloid",
                       "11"="NK", 
                       "12"="T", 
                       "13"= "Myeloid", 
                       "14"= "T", 
                       "15"= "T",
                       "16"= "Plasma", 
                       "17"= "T", 
                       "18"= "B",
                       "19"= "Myeloid",
                       "20"="T",
                       "21"="Myeloid", 
                       "22"="Myeloid", 
                       "23"= "T",
                       "24"= "T", 
                       "25"= "Plasma",
                       "26"= "T", 
                       "27"= "B", 
                       "28"= "Epithelial",
                       "29"= "NK",
                       "30"="T",
                       "31"="B", 
                       "32"="NK", 
                       "33"= "Plasma", 
                       "34"= "B", 
                       "35"= "Mast",
                       "36"= "Myeloid", 
                       "37"= "T", 
                       "38"= "T",
                       "39"= "T",
                       "40"="T",
                       "41"="Stromal", 
                       "42"="T", 
                       "43"= "Myeloid",
                       "44"= "Stromal", 
                       "45"= "B",
                       "46"= "T", 
                       "47"= "Stromal",
                       "48"= "NK"
)
sce[['celltype']] = unname(cluster2celltype[sce@meta.data$seurat_clusters])
View(sce@meta.data) 

##celltype UMAP图
pdf(paste0("./",sam.name,"/celltype-DimPlot_umap_",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "umap",group.by = 'celltype',label = T)
dev.off()

#### 2.1 计算major_marker gene
sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
write.table(sce.markers,file=paste0("./",sam.name,"/","seurat_clusters",max(dim.use),"PC.txt"),sep="\t",quote = F)

#补充meta.data数据
ldf=read.table("patient.txt",header = T)[,1]
sce@meta.data$patient <- ldf
ldf=read.table("tissue.txt",header = T)[,1]
sce@meta.data$tissue <- ldf
ldf=read.table("age.txt",header = T)[,1]
sce@meta.data$age <- ldf
ldf=read.table("sex.txt",header = T)[,1]
sce@meta.data$sex<- ldf
ldf=read.table("sample.txt",header = T)[,1]
sce@meta.data$sample <- ldf
ldf=read.table("TNM_T.txt",header = T)[,1]
sce@meta.data$TNM_T <- ldf
ldf=read.table("TNM_N.txt",header = T)[,1]
sce@meta.data$TNM_N <- ldf
ldf=read.table("TNM_M.txt",header = T)[,1]
sce@meta.data$TNM_M <- ldf
ldf=read.table("stage.txt",header = T)[,1]
sce@meta.data$stage <- ldf
ldf=read.table("lymph_node_metastasis.txt",header = T)[,1]
sce@meta.data$lymph_node_metastasis <- ldf

sce@meta.data$sourse<-"GSE164522"
rm(ldf)

#4.存储数据大类sce1
save(sce,file=paste0("./",sam.name,"/","sce1.RData"))

## 2.重新读取数据
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'sce1.RData')
View(sce@meta.data)
#确定用于细胞分群的PC,放置数据的位置"test2"
dim.use <- 1:50
sam.name <- "test2"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#按照数据来源分组展示细胞异同-sample
#sce= SetIdent(sce,value="subname")
pdf(paste0("./",sam.name,"/umap_patient.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="patient",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_tissue.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="tissue",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_TNM_T.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_T",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_TNM_N.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_N",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_TNM_M.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_M",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_sex.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="sex",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/lymph_node_metastasis.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="lymph_node_metastasis",reduction='umap')
dev.off()

features= c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","FCGR3B","VWF",
            "LILRA4","CD68","COL1A2","FCGR3B","PLVAP","CD14")

##major_marker gene 气泡图
pdf(paste0("./",sam.name,"/DotPlot_umap_celltype.",max(dim.use),"PC.pdf"),width = 10,height = 4)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()

##major_marker gene 气泡图更好的可视化
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
features <- c("TFF3","AGR2","EPCAM",#Epithelial
              "CD3E","CD3D","CCL5",#T
              "PLVAP","VWF","PECAM1",#
              "JCHAIN","MZB1","IGHA1",#plasma
              "COL1A1","COL1A2","DCN",#Fibroblast
              "LYZ","CD68","C1QC",#Myeloid
              "MS4A1","CD79A","CD19",#B
              "TPSB2","TPSAB1","CPA3"#Mast
)
pdf(paste0("./",sam.name,"/niche_marker_DotPlot_CRC1",".pdf"),width =4.5,height = 5)
DotPlot(sce, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()

load(file = 'sce1.RData')

#4.提取Mast亚群
sce_Mast<-subset(sce,celltype %in% c("Mast"))
table(sce_Mast@meta.data$celltype)

#4.存储sce_Mast数据
save(sce_Mast,file=paste0("./",sam.name,"/","sce_Mast.RData"))

#计算Mast_NormalvsTumor基因的差异
sce_Mast= SetIdent(sce_Mast,value="tissue") 
markers <- FindMarkers(sce_Mast, ident.1="NC", ident.2="CRC",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","Mast_NormalvsTumor",max(dim.use),"PC.txt"),sep="\t",quote = F)
table(sce_Mast@meta.data$tissue)

## Mast亚群重新聚类分组
sce_Mast <- NormalizeData(sce_Mast, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_Mast <- FindVariableFeatures(sce_Mast, selection.method = 'vst', nfeatures = 2000)
sce_Mast<- ScaleData(sce_Mast, vars.to.regress = "percent.mt")
sce_Mast<- RunPCA(sce_Mast, features = VariableFeatures(object = sce)) 
ElbowPlot(sce_Mast, reduction = "pca",ndims = 50)

dim.use<-1:50
sce_Mast <- FindNeighbors(sce_Mast, dims = 1:50)
sce_Mast <- FindClusters(sce_Mast, resolution = 0.8 )
sce_Mast <- RunUMAP(sce_Mast, dims = 1:50)

##sce_Mast umap图
pdf(paste0("./",sam.name,"/sce_Mast_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Mast, reduction = 'umap',label=T)
dev.off()

##sce_Mast umap图split.by = "tissue"
pdf(paste0("./",sam.name,"/sce_Mast_umap_split_tissue",max(dim.use),"PC.pdf"),width = 15,height = 4)
DimPlot(sce_Mast, split.by = "tissue",reduction = 'umap',label=T)
dev.off()

#### 计算mast_marker基因
all.markers <- FindAllMarkers(sce_Mast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

#### 计算main_marker基因
sce= SetIdent(sce,value="celltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","main_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

