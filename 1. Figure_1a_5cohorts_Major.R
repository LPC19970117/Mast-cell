rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)

#output data file
sam.name <- "test1"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#Epithelial_Count data input
experiment.data<-Read10X_h5('Epithelial_Count_matrix.h5')

#创建Seurat分析对象
sce_Epithelial <- CreateSeuratObject(
  experiment.data,
  project = "Epithelial_Count", 
  #min.cells = 10,
  #min.features = 200
  )

#NonEpithelial_Count data input
experiment.data<-Read10X_h5('NonEpithelial_Count_matrix.h5')

#创建Seurat分析对象
sce_NonEpithelial <- CreateSeuratObject(
  experiment.data,
  project = "NonEpithelial_Count", 
 #min.cells = 10,
 #min.features = 200
 )

#merge Epithelial和NonEpithelial数据
sce <- merge(sce_Epithelial, y=c(sce_NonEpithelial))
rm(sce_Epithelial,sce_NonEpithelial)
#4.存储数据sce0
save(sce,file=paste0("./",sam.name,"/","sce0.RData"))

## 2.重新读取数据
#rm(list = ls())
#library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'sce0.RData')

#### 5. QC ####原作者qc过，我就不进一步qc了
### QC——线粒体
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#画图的代码需要选中从pdf一直到dev.off()的所有代码，一起运???
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

#cat("Before filter :",nrow(sce@meta.data),"cells\n")
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
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "pca")
dev.off()
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(sce, reduction = "pca",ndims = 50)
dev.off()

### Cluster
dim.use <- 1:50
# ############################################################################################################
sce <- RunUMAP(sce, dims = dim.use) # umap tsne
sce <- FindNeighbors(sce, dims = 1:50)  # louvain cluster, graph based
sce <- FindClusters(sce, resolution = 0.8)

##seurat_clusters_DimPlot
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
features= c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","FCGR3B","VWF",
            "LILRA4","CD68","COL1A2","FCGR3B","PLVAP","CD14")
##major_marker_gene_DimPlot
pdf(paste0("./",sam.name,"/DotPlot_umap.",max(dim.use),"PC.pdf"),width = 10,height = 10)
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

#以下重命名能在meta.data中显示
#方法二：使用unname函数配合向量：
celltype <- c(         "0"="T",
                       "1"="T", 
                       "2"="B", 
                       "3"= "T", 
                       "4"= "Myeloid", 
                       "5"= "T",
                       "6"= "Fibroblast", 
                       "7"= "Epithelial", 
                       "8"= "Myeloid",
                       "9"= "T",
                       "10"="Plasma",
                       "11"="Endothelial", 
                       "12"="Fibroblast", 
                       "13"= "Fibroblast", 
                       "14"= "Fibroblast", 
                       "15"= "Plasma",
                       "16"= "Plasma", 
                       "17"= "Epithelial", 
                       "18"= "Epithelial",
                       "19"= "Epithelial",
                       "20"="Myeloid",
                       "21"="Plasma", 
                       "22"="Plasma", 
                       "23"= "Plasma",
                       "24"= "Fibroblast", 
                       "25"= "Epithelial",
                       "26"= "Plasma", 
                       "27"= "Plasma", 
                       "28"= "Mast",
                       "29"= "T",
                       "30"="Epithelial",
                       "31"="T", 
                       "32"="Epithelial", 
                       "33"= "Plasma", 
                       "34"= "Fibroblast", 
                       "35"= "Epithelial",
                       "36"= "T", 
                       "37"= "Myeloid", 
                       "38"= "Epithelial",
                       "39"= "Fibroblast",
                       "40"="T",
                       "41"="Epithelial", 
                       "42"="Epithelial", 
                       "43"= "Plasma", 
                       "44"= "Epithelial", 
                       "45"= "Plasma",
                       "46"= "Endothelial", 
                       "47"= "Fibroblast", 
                       "48"= "T",
                       "49"= "Myeloid",
                       "50"="Epithelial",
                       "51"="Epithelial", 
                       "52"="T", 
                       "53"= "Plasma", 
                       "54"= "Plasma", 
                       "55"= "Plasma",
                       "56"= "Endothelial", 
                       "57"= "Fibroblast", 
                       "58"= "Plasma",
                       "59"= "T",
                       "60"="Plasma",
                       "61"="T"
)
sce[['celltype']] = unname(celltype[sce@meta.data$seurat_clusters])
View(sce@meta.data) 

##celltype UMAP图
pdf(paste0("./",sam.name,"/celltype_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce, reduction = "umap",group.by = 'celltype',label = T)
dev.off()

#### 2.1 计算major_marker gene
sce= SetIdent(sce,value="celltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","main_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

#输出meta.data数据，在excel中修改后补充
write.table(sce@meta.data,file=paste0("./",sam.name,"/","sce@meta.data",".txt"),sep="\t",quote = F)
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
ldf=read.table("site.txt",header = T)[,1]
sce@meta.data$site <- ldf
ldf=read.table("MSI.txt",header = T)[,1]
sce@meta.data$MSI<- ldf
ldf=read.table("source.txt",header = T)[,1]
sce@meta.data$source<- ldf
rm(ldf)

##临床参数的UMAP图
pdf(paste0("./",sam.name,"/umap_patient.",max(dim.use),"PC.pdf"),width = 8,height = 4)
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
pdf(paste0("./",sam.name,"/umap_source.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="source",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_site.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="site",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_stage.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="stage",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_MSI.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="MSI",reduction='umap')
dev.off()

#4.存储数据sce1
save(sce,file=paste0("./",sam.name,"/","sce1.RData"))

##输出样本临床和细胞数，用于后续比较分析
table(sce@meta.data$sample,sce@meta.data$celltype)
table(sce@meta.data$celltype,sce@meta.data$tissue)
table(sce@meta.data$sample,sce@meta.data$tissue)
write.csv(table(sce@meta.data$celltype,sce@meta.data$tissue),file="celltype_tissue.csv")
write.csv(table(sce@meta.data$sample,sce@meta.data$celltype),file="sample_celltype.csv")
write.csv(table(sce@meta.data$sample,sce@meta.data$tissue),file="sample_tissue.csv")

##major_marker_gene点状热图
features= c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","FCGR3B","VWF",
            "LILRA4","CD68","COL1A2","FCGR3B","PLVAP","CD14")
pdf(paste0("./",sam.name,"/DotPlot_umap_celltype.",max(dim.use),"PC.pdf"),width = 10,height = 4)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()

##major_marker_gene点状热图
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

#4.提取mast_cell单细胞亚群
sce_Mast<-subset(sce,celltype %in% c("Mast"))

#dim(sce)
#rm(sce_Macrophage)
#4.存储数据sce1
save(sce_Mast,file=paste0("./",sam.name,"/","sce_Mast.RData"))

#Mast_NormalvsTumor基因的差异
sce_Mast= SetIdent(sce_Mast,value="tissue") 
markers <- FindMarkers(sce_Mast, ident.1="Normal", ident.2="Tumor",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","Mast_NormalvsTumor",max(dim.use),"PC.txt"),sep="\t",quote = F)
table(sce_Mast@meta.data$tissue)

## Mast_cell重新聚类分组
sce_Mast <- NormalizeData(sce_Mast, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_Mast <- FindVariableFeatures(sce_Mast, selection.method = 'vst', nfeatures = 2000)
sce_Mast<- ScaleData(sce_Mast, vars.to.regress = "percent.mt")
sce_Mast<- RunPCA(sce_Mast, features = VariableFeatures(object = sce)) 
ElbowPlot(sce_Mast, reduction = "pca",ndims = 50)

dim.use<-1:50
sce_Mast <- FindNeighbors(sce_Mast, dims = 1:50)
sce_Mast <- FindClusters(sce_Mast, resolution = 0.8 )
sce_Mast <- RunUMAP(sce_Mast, dims = 1:50)

##umap
pdf(paste0("./",sam.name,"/sce_Mast_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Mast, reduction = 'umap',label=T)
dev.off()

#umap split.by = "tissue"
pdf(paste0("./",sam.name,"/sce_Mast_umap_split_tissue",max(dim.use),"PC.pdf"),width = 15,height = 4)
DimPlot(sce_Mast, split.by = "tissue",reduction = 'umap',label=T)
dev.off()
#### 2.1 计算mast_marker基因
all.markers <- FindAllMarkers(sce_Mast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

