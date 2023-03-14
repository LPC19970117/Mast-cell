rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)
sam.name <- "肥大细胞特征基因热图2"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

load(file = 'sce_Mast.RData')

sce_Mast<-subset(sce_Mast,tissue %in% c("Normal","Tumor"))




## 重新聚类分组
sce_Mast <- NormalizeData(sce_Mast, normalization.method = "LogNormalize", scale.factor = 10000) 
sce_Mast <- FindVariableFeatures(sce_Mast, selection.method = 'vst', nfeatures = 2000)
#all.genes <- rownames(sce) 
#sce <- ScaleData(sce, features = all.genes, vars.to.regress = c("percent.mt"))
sce_Mast<- ScaleData(sce_Mast, vars.to.regress = "percent.mt")
sce_Mast<- RunPCA(sce_Mast, features = VariableFeatures(object = sce_Mast)) 
#pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
#DimPlot(sce, reduction = "pca")
#dev.off()
#pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
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
#### 2.1 计算mast_marker基因
sce_Mast= SetIdent(sce_Mast,value="seurat_clusters") 
all.markers <- FindAllMarkers(sce_Mast, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes","PC.txt"),sep="\t",quote = F)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes_noMIT","PC.txt"),sep="\t",quote = F)

#Mast_NormalvsTumor基因的差异
sce_Mast= SetIdent(sce_Mast,value="tissue") 
markers <- FindMarkers(sce_Mast, ident.1="N", ident.2="T",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","Mast_NormalvsTumor_count",max(dim.use),"PC.txt"),sep="\t",quote = F)


markers <- FindMarkers(sce_Mast, ident.1="Normal", ident.2="Tumor",
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","Mast_NormalvsTumor-2",max(dim.use),"PC.txt"),sep="\t",quote = F)

## 气泡图
features <- c("IL4","IL5","IL9","IL13","IL18",
              "CCL1",
              #"TNF","CSF2","IFNG","KITLG",
              "LIF","CSF1","AREG","VEGFA","TGFB1",#n= 15cytokines and GF
              "TPSAB1","TPSB2","CPA3","CTSG","CTSD","CTSW","CMA1","HDC",#n=8 proteases 和组胺
              "LTC4S","ALOX5AP","ALOX5","HPGDS","PTGS2","PTGS1",#n=6 lipid mediator
              "KIT","FCER1G","FCER1A","MS4A2","MRGPRX2","CSF2RB","IL1RL1","AHR"#n=8 receptor
              
)
pdf(paste0("./",sam.name,"/mast_character_gene_umap.",max(dim.use),"PC.pdf"),width = 10,height = 2)
DotPlot(sce_Mast, features = unique(features)) + RotatedAxis()
dev.off()
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_character_gene_umap2",".pdf"),width =8,height = 12)
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

#输出meta.data数据，在excel中修改后补充
sce_Mast= SetIdent(sce_Mast,value="sample") 
heatmap_AveE <- AverageExpression(sce_Mast, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA

write.table(sce_Mast@meta.data,file=paste0("./",sam.name,"/","sce_mast@meta.data",".txt"),sep="\t",quote = F)
#补充meta.data数据
ldf=read.table("sample2.txt",header = T)[,1]
sce_Mast@meta.data$sample2 <- ldf

#输出heatmap_AveE数据，在excel中排序后输入
write.table(heatmap_AveE,file=paste0("./",sam.name,"/","heatmap_AveE",".txt"),sep="\t",quote = F)

heatmap_AveE2 <-read.table("Mast_heatmap_input.txt",header = TRUE)
#将数据从dataframe转为matrix
heatmap_AveE2 <-data.matrix(heatmap_AveE2)
gene_num <- c(11,8,6,8)
gaps_row <- cumsum(gene_num)

cluster_num <- c(31)
gaps_col <- cumsum(cluster_num)

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
#annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
#                             `subcelltype2` = rep(factor(subcelltype2,levels = subcelltype2),gene_num))
#annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
#                             `subcelltype2` = rep(factor(subcelltype2,levels = subcelltype2),cluster_num))
#annotation_colors = list(`subcelltype2` = subcelltype2_colors)
#names(annotation_colors$`subcelltype2`) = subcelltype2
pdf(paste0("./",sam.name,"/heatmap20220508_sample2.pdf"),width = 4,height = 8)
pheatmap(heatmap_AveE2,cluster_cols = F,cluster_rows = F,show_colnames=T,show_rownames=T,
         border=F,#border_color = "white",
         color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                   colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
         breaks=bk,scale="row",legend_breaks=seq(-2,2,2),
         gaps_row = gaps_row,gaps_col = gaps_col,
         #annotation_row = annotation_row,annotation_col = annotation_col,
         #annotation_colors = annotation_colors,
         annotation_names_row = F,annotation_names_col = T)
dev.off()

