rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)
sam.name <- "肥大细胞特征基因热图_主要细胞2（重制版）"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
#read major cell data (5_cohorts)
load(file = 'sce1.RData')

##
features <- c("IL4","IL5","IL9","IL13","IL18",
              "CCL1",
              #"TNF","CSF2","IFNG","KITLG",
              "LIF","CSF1","AREG","VEGFA","TGFB1",#n= 15cytokines and GF
              "TPSAB1","TPSB2","CPA3","CTSG","CTSD","CTSW","CMA1","HDC",#n=8 proteases 和组胺
              "LTC4S","ALOX5AP","ALOX5","HPGDS","PTGS2","PTGS1",#n=6 lipid mediator
              "KIT","FCER1G","FCER1A","MS4A2","MRGPRX2","CSF2RB","IL1RL1","AHR"#n=8 receptor
              
)


######################热图Fiugre2d
sce= SetIdent(sce,value="celltype") 
heatmap_AveE <- AverageExpression(sce, assays = "RNA", features = features,verbose = TRUE) %>% .$RNA

gene_num <- c(15,8,6,8)
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
pdf(paste0("./",sam.name,"/heatmap20220508_2.pdf"),width = 3,height = 7)
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

#Figure5c_left
library(RColorBrewer)
color1 = colorRampPalette(rev(brewer.pal(n= 7, name = "RdYlBu")))(100)
pdf(paste0("./",sam.name,"/4gene_umap_",".pdf"),width =4,height = 1.8)
DotPlot(sce, features = c("KITLG","KIT","IL33","IL1RL1"))+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = color1)+ RotatedAxis()
dev.off()

#Figure5c_Right
pdf(paste0("./",sam.name,"/双基因图KIT.pdf"),width = 10,height = 3)
FeaturePlot(sce,features=c("KITLG","KIT"),cols = c("#2166ac","#b2182b"),blend=TRUE)
dev.off()
pdf(paste0("./",sam.name,"/双基因图IL33.pdf"),width = 10,height = 3)
FeaturePlot(sce,features=c("IL33","IL1RL1"),blend=TRUE)
dev.off()
pdf(paste0("./",sam.name,"/双基因图IL33_0.5.pdf"),width = 10,height = 3)
FeaturePlot(sce,features=c("IL33","IL1RL1"),blend.threshold = 0,blend=TRUE)
dev.off()


#Table.S7A
sce_endo<-subset(sce,celltype %in% c("Endothelial"))
#Endothelial_NormalvsTumor基因的差异
sce_endo= SetIdent(sce_endo,value="tissue") 

markers <- FindMarkers(sce_endo,  ident.1="Normal", ident.2="Tumor",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","sce_endo_NormalvsTumor","PC.txt"),sep="\t",quote = F)
table(sce_endo@meta.data$tissue)

pdf(paste0("./",sam.name,"/Endo_DotPlot_KITLG.","PC.pdf"),width = 7,height =8)
DotPlot(sce_endo, features = c("KITLG","KIT","IL33","IL1RL1")) + RotatedAxis()
dev.off()

#Table.S7B
sce_fibro<-subset(sce,celltype %in% c("Fibroblast"))
sce_fibro= SetIdent(sce_fibro,value="tissue") 

markers <- FindMarkers(sce_fibro, ident.1="Normal", ident.2="Tumor",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","sce_fibro_NormalvsTumor","PC.txt"),sep="\t",quote = F)
table(sce_endo@meta.data$tissue)
pdf(paste0("./",sam.name,"/Fibro_DotPlot_KITLG.","PC.pdf"),width = 7,height =8)
DotPlot(sce_fibro, features = c("KITLG","KIT")) + RotatedAxis()
dev.off()
