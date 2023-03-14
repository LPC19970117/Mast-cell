rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)
sam.name <- "Mast_angiogenesis"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

load(file = 'sce_Mast2.RData')

Angiogenesis_score <- readxl::read_xlsx("Angiogenesis_score.xlsx")
#转换成list
sce_Mast= SetIdent(sce_Mast,value="subcelltype2") 
#4.提取指定单细胞亚群
sce_Mast<-subset(sce_Mast,subcelltype2 %in% c("Activated_MC","Resting_MC"))

gene <- as.list(Angiogenesis_score)
sce_Mast <- AddModuleScore(
  object = sce_Mast,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'Angiogenesis_score'
)
pdf(paste0("./",sam.name,"/Angiogenesis_score",".pdf"),width = 4,height = 4)
VlnPlot(sce_Mast,features = 'Angiogenesis_score1',pt.size = 0
        )+geom_boxplot(width=0.2,col="black",fill="white")
dev.off()
write.table(sce_Mast@meta.data,"sce_Mast@meta.data_Angiogenesis",sep="\t",quote = F)
