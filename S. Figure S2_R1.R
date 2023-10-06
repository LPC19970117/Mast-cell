rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)
#make a file for output data
sam.name <- "Mast_umap2"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
##read mast cell data in 5-cohort
load(file = 'sce_Mast.RData')


## mast cell data in GSE178341 Recluster grouping
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


######
pdf(paste0("./",sam.name,"/umap.",max(dim.use),"_mast_gene_umap.pdf"),width = 25,height = 24)
FeaturePlot(sce_Mast, features = features,
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()
#ClusterTree
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
pdf(paste0("./",sam.name,"/umap_source.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="source",reduction='umap')
dev.off()
######
pdf(paste0("./",sam.name,"/umap.",max(dim.use),"_mast_gene_umap.pdf"),width = 25,height = 24)
FeaturePlot(sce, features = features,
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()


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
pdf(paste0("./",sam.name,"/mast_character_seurat_clusters_umap.",max(dim.use),"PC.pdf"),width = 10,height = 8)
DotPlot(sce_Mast,group.by = 'seurat_clusters', features = unique(features)) + RotatedAxis()
dev.off()

color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_seurat_clusters.pdf"),width =6,height =7)
DotPlot(sce_Mast, group.by = 'seurat_clusters',features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)
dev.off()

all.markers <- FindAllMarkers(sce_Mast, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes_subcelltype2_","PC.txt"),sep="\t",quote = F)


####计算每个样本MC的seurat_clusters差异#### Calculate the seurat_clusters difference of each sample MC
write.csv(table(sce_Mast@meta.data$seurat_clusters,sce_Mast@meta.data$tissue),file="MC_seurat_clusters_tissue.csv")
write.csv(table(sce_Mast@meta.data$seurat_clusters,sce_Mast@meta.data$sample),file="MC_seurat_clusters_sample.csv")

cluster5celltype <-c("0"="MC10",
                     "1"="MC06",
                     "2"="MC08",
                     "3"="MC03",
                     "4"="MC09",
                     "5"="MC12",
                     "6"="MC11",
                     "7"="MC04",
                     "8"="MC01",
                     "9"="MC02",
                     "10"="MC05",
                     "11"="MC15",
                     "12"="MC14",
                     "13"="MC13",
                     "14"="MC07"
)
sce_Mast[['subcelltype_new']] = unname(cluster5celltype[sce_Mast@meta.data$seurat_clusters])
sce_Mast= SetIdent(sce_Mast,value="subcelltype_new") 

cluster5celltype <-c("0"="Resting_MC",
                     "1"="Activated_MC",
                     "2"="Activated_MC",
                     "3"="Activated_MC",
                     "4"="Activated_MC",
                     "5"="Resting_MC",
                     "6"="Resting_MC",
                     "7"="Activated_MC",
                     "8"="Activated_MC",
                     "9"="Activated_MC",
                     "10"="Activated_MC",
                     "11"="Other_MC",
                     "12"="Other_MC",
                     "13"="Proliferating_MC",
                     "14"="Activated_MC"
)
sce_Mast[['subcelltype2']] = unname(cluster5celltype[sce_Mast@meta.data$seurat_clusters])
sce_Mast= SetIdent(sce_Mast,value="subcelltype2") 

##Figure_3C
features <- c("IL4","IL5","IL9","IL13","IL18",
              "CCL1",
              #"TNF","CSF2","IFNG","KITLG",
              "LIF","CSF1","AREG","VEGFA","TGFB1",#n= 15cytokines and GF
              "TPSAB1","TPSB2","CPA3","CTSG","CTSD","CTSW","CMA1","HDC",#n=8 proteases 和组胺
              "LTC4S","ALOX5AP","ALOX5","HPGDS","PTGS2","PTGS1",#n=6 lipid mediator
              "KIT","FCER1G","FCER1A","MS4A2","MRGPRX2","CSF2RB","IL1RL1","AHR"#n=8 receptor
)
features <- c("CD81","IGHA2","CPA3","WASHC1","MS4A1",
              "TPSD1",
              "CLC","TPSB2","CREM","FTH1","KRT1",#n= 15cytokines and GF
              "CMA1","MKI67","MZB1","CD3D"
)
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_subcelltype_new.pdf"),width =6,height =7)
DotPlot(sce_Mast, group.by = 'subcelltype_new',features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+ RotatedAxis()
dev.off()

library(RColorBrewer)
pdf(paste0("./",sam.name,"/mast_character_subcelltype2.",max(dim.use),"PC.pdf"),width = 10,height = 8)
DotPlot(sce_Mast,group.by = 'subcelltype2', features = unique(features)) + RotatedAxis()
dev.off()

color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_subcelltype2.pdf"),width =6,height =7)
DotPlot(sce_Mast, group.by = 'subcelltype2',features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+ RotatedAxis()
dev.off()

pdf(paste0("./",sam.name,"/subcelltype2.pdf"),width = 5,height = 4)
DimPlot(object = sce_Mast, group.by="subcelltype2",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/subcelltype_new.pdf"),width = 5,height = 4)
DimPlot(object = sce_Mast, group.by="subcelltype_new",label=T,reduction='umap')
dev.off()
#MC_resting_vs MC_activated基因的差异
#
sce_Mast= SetIdent(sce_Mast,value="subcelltype2") 
markers <- FindMarkers(sce_Mast, ident.1="Resting_MC", ident.2="Activated_MC",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","MC_resting_vs MC_activated",max(dim.use),"PC.txt"),sep="\t",quote = F)

sce_Mast= SetIdent(sce_Mast,value="subcelltype_new") 
all.markers <- FindAllMarkers(sce_Mast, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(all.markers,file=paste0("./",sam.name,"/","sce_Mast_marker_genes_subcelltype_new.txt"),sep="\t",quote = F)


features <- c("CD81","CTSD","CPA3","WASHC1","ENO2",
              "TPSD1",
              "CLC","TPSB2","TPSAB1","FAU","KRT1",#n= 15cytokines and GF
              "CMA1","MKI67","MZB1","CD3D"
)
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/mast_subcelltype_new_marker2.pdf"),width =5.5,height =4)
DotPlot(sce_Mast, group.by = 'subcelltype_new',features = features)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+ RotatedAxis()
dev.off()

####组织中不同subcelltype2####Different subcelltype2 in the organization
Ratio <- sce_Mast@meta.data %>%group_by(tissue,subcelltype2) %>%
  count() %>%
  group_by(tissue) %>%
  mutate(Freq = n/sum(n)*100)

pdf(paste0("./",sam.name,"/Ratio_3.",max(dim.use),"PC.pdf"),width = 4,height = 8)
ggplot(Ratio, aes(x = tissue, y = Freq, fill = subcelltype2))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")), 
            position = position_stack(vjust = 0.5))+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                               "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                               "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
dev.off()


