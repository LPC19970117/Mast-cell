rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)
#make a file for output data
sam.name <- "肥大细胞特征基因热图3"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

##read mast cell data in GSE178341
load(file = 'sce_Mast.RData')

## mast cell data in GSE178341重新聚类分组
sce_Mast <- NormalizeData(sce_Mast, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_Mast <- FindVariableFeatures(sce_Mast, selection.method = 'vst', nfeatures = 2000)
sce_Mast<- ScaleData(sce_Mast, vars.to.regress = "percent.mt")
sce_Mast<- RunPCA(sce_Mast, features = VariableFeatures(object = sce_Mast)) 
ElbowPlot(sce_Mast, reduction = "pca",ndims = 50)
#dev.off()
dim.use<-1:50
sce_Mast <- FindNeighbors(sce_Mast, dims = 1:50)
sce_Mast <- FindClusters(sce_Mast, resolution = 0.8 )
sce_Mast <- RunUMAP(sce_Mast, dims = 1:50)
pdf(paste0("./",sam.name,"/sce_Mast_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce_Mast, reduction = 'umap',label=T)
dev.off()
pdf(paste0("./",sam.name,"/sce_Mast_umap_split_tissue",max(dim.use),"PC.pdf"),width = 8,height = 4)
DimPlot(sce_Mast, split.by = "tissue",reduction = 'umap',label=T)
dev.off()

## Figure_S1A
features <- c("FTH1","LTC4S","PTGS2","PIGR",
              "IL1RL1","CMA1","MZB1","CD3D",#n= 15cytokines and GF
              "TPSAB1","IFIT1","MKI67","MT-CO2")
pdf(paste0("./",sam.name,"/mast_character_gene_umap.",max(dim.use),"PC.pdf"),width = 10,height = 8)
DotPlot(sce_Mast, features = unique(features)) + RotatedAxis()
dev.off()
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_character_gene_umap2",".pdf"),width =5,height =3)
DotPlot(sce_Mast, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()

######
pdf(paste0("./",sam.name,"/umap.",max(dim.use),"_mast_gene_umap.pdf"),width = 25,height = 24)
FeaturePlot(sce_Mast, features = features,
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()
#细胞类群相似树
sce <- BuildClusterTree(
  sce_Mast,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/Mast-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()
pdf(paste0("./",sam.name,"/umap_tissue_mast_.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce_Mast, group.by="tissue",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_patient.",max(dim.use),"PC.pdf"),width = 8,height = 4)
DimPlot(object = sce, group.by="patient",reduction='umap')
dev.off()
######################热图
sce_Mast= SetIdent(sce_Mast,value="tissue") 
heatmap_AveE <- AverageExpression(sce_Mast, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA
gene_num <- c(11,8,6,8)
gaps_row <- cumsum(gene_num)

cluster_num <- c(0)
gaps_col <- cumsum(cluster_num)
######################热图
sce_Mast= SetIdent(sce_Mast,value="sample") 
heatmap_AveE <- AverageExpression(sce_Mast, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA

gene_num <- c(11,8,6,8)
gaps_row <- cumsum(gene_num)
cluster_num <- c(0)
gaps_col <- cumsum(cluster_num)

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
#annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
#                             `subcelltype2` = rep(factor(subcelltype2,levels = subcelltype2),gene_num))
#annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
#                             `subcelltype2` = rep(factor(subcelltype2,levels = subcelltype2),cluster_num))
#annotation_colors = list(`subcelltype2` = subcelltype2_colors)
#names(annotation_colors$`subcelltype2`) = subcelltype2
pdf(paste0("./",sam.name,"/heatmap20220508_sample.pdf"),width = 30,height = 15)
pheatmap(heatmap_AveE,cluster_cols = F,cluster_rows = F,show_colnames=T,show_rownames=T,
         border=F,#border_color = "white",
         color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                   colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
         breaks=bk,scale="row",legend_breaks=seq(-2,2,2),
         gaps_row = gaps_row,gaps_col = gaps_col,
         #annotation_row = annotation_row,annotation_col = annotation_col,
         #annotation_colors = annotation_colors,
         annotation_names_row = F,annotation_names_col = T)
dev.off()

重命名
#方法二：使用unname函数配合向量：
cluster5celltype <-c("0"="MC06_FTH",
                     "1"="MC04_LTC4S",
                     "2"="MC03_PTGS2",
                     "3"="MC07_PIGR",
                     "4"="MC01_IL1RL1",
                     "5"="MC05_CMA1",
                     "6"="MC11_MZB1",
                     "7"="MC12_CD3D",
                     "8"="MC02_TPSAB1",
                     "9"="MC09_IFITs",
                     "10"="MC08_MKI67",
                     "11"="MC10_MTs"
                     
)
sce_Mast[['subcelltype']] = unname(cluster5celltype[sce_Mast@meta.data$seurat_clusters])
View(sce_Mast@meta.data) 

pdf(paste0("./",sam.name,"/Mast_subcelltype_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce_Mast, reduction = "umap",group.by = 'subcelltype',label = T)
dev.off()

cluster5celltype <-c("0"="Resting_MC",
                     "1"="Activated_MC",
                     "2"="Activated_MC",
                     "3"="Resting_MC",
                     "4"="Activated_MC",
                     "5"="Resting_MC",
                     "6"="Other_MC",
                     "7"="Other_MC",
                     "8"="Activated_MC",
                     "9"="Other_MC",
                     "10"="Proliferating_MC",
                     "11"="Other_MC"
                     
)
sce_Mast[['subcelltype2']] = unname(cluster5celltype[sce_Mast@meta.data$seurat_clusters])

##Figure_3A_left
pdf(paste0("./",sam.name,"/Mast_subcelltype2_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce_Mast, reduction = "umap",group.by = 'subcelltype2',label = T)
dev.off()


sce_Mast= SetIdent(sce_Mast,value="subcelltype2") 
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_subcelltype2_",".pdf"),width =6,height =10)
DotPlot(sce_Mast, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()
#### 2.1 计算mast_marker基因
#先去除Mt等干扰基因
GeneMT <- read.table("MT.txt",stringsAsFactors = F)
GeneHSP <- read.table("HSP.txt",stringsAsFactors = F)
GeneRP <- read.table("RP.txt",stringsAsFactors = F)
GeneSC <- read.table("sc_dissociation.txt",stringsAsFactors = F)

GenesRM <- rbind(GeneMT, GeneHSP, GeneRP,GeneSC)
GenesRM <- GenesRM$V1
# gene filter 
sce_Mast2 <- sce_Mast[!rownames(sce_Mast) %in% GenesRM,]
all.markers <- FindAllMarkers(sce_Mast, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes_subcelltype2_","PC.txt"),sep="\t",quote = F)

sce_Mast= SetIdent(sce_Mast,value="subcelltype") 
library(RColorBrewer)
pdf(paste0("./",sam.name,"/mast_character_subcelltype_umap.",max(dim.use),"PC.pdf"),width = 10,height = 8)
DotPlot(sce_Mast, features = unique(features)) + RotatedAxis()
dev.off()
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_subcelltype_",".pdf"),width =20,height =10)
DotPlot(sce_Mast, features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()

#MC_resting_vs MC_activated基因的差异
#Table_S6
sce_Mast2= SetIdent(sce_Mast2,value="subcelltype2") 
markers <- FindMarkers(sce_Mast2, ident.1="Resting_MC", ident.2="Activated_MC",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","MC_resting_vs MC_activated",max(dim.use),"PC.txt"),sep="\t",quote = F)

#方法二：使用unname函数配合向量：
cluster5celltype <-c("0"="MC06",
                     "1"="MC04",
                     "2"="MC03",
                     "3"="MC07",
                     "4"="MC01",
                     "5"="MC05",
                     "6"="MC11",
                     "7"="MC12",
                     "8"="MC02",
                     "9"="MC09",
                     "10"="MC08",
                     "11"="MC10")
sce_Mast[['subcelltype3']] = unname(cluster5celltype[sce_Mast@meta.data$seurat_clusters])
sce_Mast= SetIdent(sce_Mast,value="subcelltype3") 

##Figure_3C
features <- c("IL4","IL5","IL9","IL13","IL18",
              "CCL1",
              #"TNF","CSF2","IFNG","KITLG",
              "LIF","CSF1","AREG","VEGFA","TGFB1",#n= 15cytokines and GF
              "TPSAB1","TPSB2","CPA3","CTSG","CTSD","CTSW","CMA1","HDC",#n=8 proteases 和组胺
              "LTC4S","ALOX5AP","ALOX5","HPGDS","PTGS2","PTGS1",#n=6 lipid mediator
              "KIT","FCER1G","FCER1A","MS4A2","MRGPRX2","CSF2RB","IL1RL1","AHR"#n=8 receptor
)
library(RColorBrewer)
pdf(paste0("./",sam.name,"/mast_character_subcelltype3_umap.",max(dim.use),"PC.pdf"),width = 10,height = 8)
DotPlot(sce_Mast,group.by = 'subcelltype3', features = unique(features)) + RotatedAxis()
dev.off()

color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_subcelltype3_",".pdf"),width =5,height =7)
DotPlot(sce_Mast, group.by = 'subcelltype3',features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()

####计算每个样本MC的subcelltype差异
table(sce_Mast@meta.data$subcelltype,sce_Mast@meta.data$tissue)
table(sce_Mast@meta.data$subcelltype,sce_Mast@meta.data$sample)
write.csv(table(sce_Mast@meta.data$subcelltype,sce_Mast@meta.data$tissue),file="MC_subcelltype_tissue.csv")
write.csv(table(sce_Mast@meta.data$subcelltype,sce_Mast@meta.data$sample),file="MC_subcelltype_sample.csv")

#输出MC meta.data数据
write.table(sce_Mast@meta.data,file=paste0("./",sam.name,"/","sce_Mast@meta.data",".txt"),sep="\t",quote = F)
#4.存储数据
save(sce_Mast,file=paste0("./sce_Mast2.RData"))

######subcelltype特征gene图
cluster5celltype <-c("0"="MC06_FTH1",
                     "1"="MC04_LTC4S",
                     "2"="MC03_PTGS2",
                     "3"="MC07_PIGR",
                     "4"="MC01_IL1RL1",
                     "5"="MC05_CMA1",
                     "6"="MC11_MZB1",
                     "7"="MC12_CD3D",
                     "8"="MC02_TPSAB1",
                     "9"="MC09_IFITs",
                     "10"="MC08_MKI67",
                     "11"="MC10_MTs"
                     
)
features <- c("IL1RL1","TPSAB1","PTGS2","LTC4S",
              "CMA1","FTH1","PIGR","MKI67",#n= 15cytokines and GF
              "IFIT1","MT-CO2","MZB1","CD3D"
)
#Figure_S1A
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./mast_subcelltype",".pdf"),width =5,height =3.5) 
DotPlot(sce_Mast, group.by = 'subcelltype',features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+ RotatedAxis()+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()
