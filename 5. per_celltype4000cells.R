rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(hdf5r)
library(pheatmap)
sam.name <- "GSE178341抽样"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

load(file = 'sce1.RData')

# 每个细胞亚群抽50 # Take 50 samples from each cell subpopulation
sce= SetIdent(sce,value="celltype")
allCells=names(Idents(sce))
allType = levels(Idents(sce))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce)== x ]
  cg=sample(cgCells,4000)
  cg
}))
cg_sce = sce[, allCells %in% choose_Cells]
cg_sce
as.data.frame(table(Idents(cg_sce)))
#输出meta.data数据，在excel中修改后补充#Output meta.data data, modify it in excel and add it
write.table(cg_sce@meta.data,file=paste0("./",sam.name,"/","cg_sce@meta.data",".txt"),sep="\t",quote = F)

sce<-cg_sce
#4.存储数据大类sce1#4. Storage data category sce1
save(sce,file=paste0("./",sam.name,"/","sce抽样.RData"))
