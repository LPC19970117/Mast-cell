
##载入SpaceRanger数据
rm(list = ls())
#载入R包；
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(hdf5r)
# Error in library(hdf5r) : 不存在叫‘hdf5r’这个名字的程辑包
# 缺什么就安装什么
# BiocManager::install("hdf5r")
#library(hdf5r)

################################################################################
CRC <- Load10X_Spatial(data.dir ="./", #本地文件目录下运行
                       filename = "filtered_feature_bc_matrix.h5",#读取h5表达谱文件
                       slice ="CRC")

#这里文件夹的名字可以修改，但最好只用英文字母
sam.name <- "Mast1"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
# An object of class Seurat 
# 17943 features across 4371 samples within 1 assay 
# Active assay: Spatial (17943 features, 0 variable features)

p1 <- VlnPlot(CRC, 
              features ="nCount_Spatial",
              pt.size = 0,
              cols ="red")
pdf(paste0("./",sam.name,"/nCount.pdf"),width = 8,height = 4.5)
SpatialFeaturePlot(CRC, features ="nCount_Spatial")
dev.off()
p1 | p2

##################################################
## spatial天生就是Seurat对象
# 数据标准化，使用SCTransform方法进行标准化

CRC <- SCTransform(CRC, assay = "Spatial", verbose = FALSE)
CRC <- RunPCA(CRC, assay = "SCT", verbose = FALSE) 
plot1 <- DimPlot(CRC, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(CRC, ndims=20, reduction="pca") 
plot1+plot2
#选取主成分
pc.num=1:15
## 细胞聚类
# 图聚类分群
CRC <- FindNeighbors(CRC, reduction = "pca", dims = pc.num)
CRC <- FindClusters(CRC, verbose = FALSE)
# UMAP降维可视化
CRC <- RunUMAP(CRC, reduction = "pca", dims =  pc.num)
p1 <- DimPlot(CRC, reduction = "umap", label = TRUE)
# 使用SpatialDimPlot函数进行可视化
p2 <- SpatialDimPlot(CRC, label = TRUE, label.size = 3)
p1 + p2
pdf(paste0("./",sam.name,"/cluster.pdf"),width = 8,height = 4.5)
p1 + p2
dev.off()
save(CRC,file = "CRC.Rds")
View(CRC@meta.data)
################################################################################
##Figure_6A
pdf(paste0("./",sam.name,"/MC_3markers.pdf"),width = 15,height = 5)
SpatialFeaturePlot(CRC, features =c("COL1A2","PECAM1","TPSB2"),
)
dev.off()
pdf(paste0("./",sam.name,"/MC_5markers_Dot.pdf"),width = 5,height = 5)
DotPlot(CRC, features =c("TPSAB1","TPSB2","MS4A2","CPA3","HPGDS","CMA1","KIT","KITLG","IL33","IL1RL1","PECAM1","CLDN5","IGHE"))+ RotatedAxis()
dev.off()
pdf(paste0("./",sam.name,"/MC_5markers_Vln.pdf"),width = 15,height = 15)
VlnPlot(CRC, features =c("TPSAB1","TPSB2","MS4A2","CPA3","HPGDS","CMA1","KIT","KITLG","IL33","IL1RL1","PECAM1","CLDN5"))+ RotatedAxis()
dev.off()

##Figure_6C
#MC_activation_signature
#4.计算MC_activation_signature
library(tidyverse)
library(Matrix)
library(cowplot)
MC_activation_signature <- readxl::read_xlsx("MC_activation_signature.xlsx")
View(MC_activation_signature)
#转换成list
gene <- as.list(MC_activation_signature)
CRC <- AddModuleScore(
  object = CRC,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'MC_activation_signature'
)
#colnames(CRC@meta.data)["26"] <- 'MC_activation_signature' 
pdf(paste0("./",sam.name,"/MC_activation_signature_Vln.pdf"),width = 5,height = 4)
VlnPlot(CRC,features = 'MC_activation_signature1',pt.size=0)
dev.off()

pdf(paste0("./",sam.name,"/MC_activation_signature.pdf"),width = 8,height = 10)
SpatialFeaturePlot(CRC, features = 'MC_activation_signature1')
dev.off()

#以下重命名髓系亚群
#方法二：使用unname函数配合向量：
Region_type <-c("0"="Tumor",
                      "1"="Tumor",
                      "2"="Stromal",
                      "3"="Tumor",
                      "4"="Normal",
                      "5"="Normal",
                      "6"="Normal",
                      "7"="Tumor",
                      "8"="Tumor",
                      "9"="Normal")
CRC[['Region']] = unname(Region_type[CRC@meta.data$seurat_clusters])
View(CRC@meta.data) 

pdf(paste0("./",sam.name,"/Region_type_",".pdf"),width = 8,height = 7)
DimPlot(CRC, reduction = "umap",group.by = 'Region',label = T)
dev.off()

##Figure_6A(Left2)
pdf(paste0("./",sam.name,"/niche_region_type_ST",".pdf"),width = 8,height = 7)
SpatialDimPlot(CRC, group.by = 'Region',label = TRUE, label.size = 3)
dev.off()

CRC= SetIdent(CRC,value="Region") 
pdf(paste0("./",sam.name,"/region_MC_5markers_Dot.pdf"),width = 6,height = 3)
DotPlot(CRC, features =c("EPCAM","UQCRB","CEACAM5","MYC","ZG16","MUC2","FCGBP","DCN","COL1A2","PECAM1","PECAM1","TPSAB1","TPSB2","CPA3","MS4A2","HPGDS","KIT","KITLG"))+ RotatedAxis()
dev.off()

features =c("EPCAM","UQCRB","CEACAM5","MYC","ZG16","MUC2","FCGBP","DCN","COL1A2",
            "PECAM1","TPSAB1","TPSB2","CPA3","MS4A2","HPGDS","KIT","KITLG")

#### 2.1 计算marker基因
CRC= SetIdent(CRC,value="Region") 
all.markers <- FindAllMarkers(CRC, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/Region","PC.txt"),sep="\t",quote = F)
###Figure_6B
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_character_gene_umap2",".pdf"),width =3.5,height = 5)
DotPlot(CRC, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()
