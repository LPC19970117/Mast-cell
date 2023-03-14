
#4.计算signature
library(tidyverse)
library(Matrix)
library(cowplot)
MC_activation_signature <- readxl::read_xlsx("MC_activation_signature.xlsx")
View(MC_activation_signature)
#转换成list
gene <- as.list(MC_activation_signature)
sce_Mast_Mast <- AddModuleScore(
  object = sce_Mast_Mast,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'MC_activation_signature'
)

colnames(sce_Mast_Mast@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"   
# [4] "percent.mt"      "RNA_snn_res.0.5" "seurat_clusters"
# [7] "cell_type"       "CD_Features1" 
#colnames(sce_Mast@meta.data)["26"] <- 'MC_activation_signature' 
sce_Mast_Mast= SetIdent(sce_Mast_Mast,value="seurat_clusters") 
pdf(paste0("./",sam.name,"MC_activation_signature",".pdf"),width = 5,height = 4)
VlnPlot(sce_Mast_Mast,features = 'MC_activation_signature1')
dev.off()
pdf(paste0("./",sam.name,"/MC_activation_signature_VlnPlot.",max(dim.use),"PC.pdf"),width = 4,height = 5)
VlnPlot(sce_Mast_Mast, features ='MC_activation_signature1') + RotatedAxis()
dev.off()
sce_Mast_Mast= SetIdent(sce_Mast_Mast,value="tissue") 
pdf(paste0("./",sam.name,"/MC_activation_signature_tissue_VlnPlot_.",max(dim.use),"PC.pdf"),width = 4,height = 5)
VlnPlot(sce_Mast_Mast, features ='MC_activation_signature1') + RotatedAxis()
dev.off()

pdf(paste0("./",sam.name,"/MC_activation_signature_VlnPlot2.",max(dim.use),"PC.pdf"),width = 4,height = 5)
VlnPlot(sce_Mast_Mast, features ='MC_activation_signature1') 
dev.off()

pdf(paste0("./",sam.name,"/MC_activation_signature_FeaturePlot.",max(dim.use),"PC.pdf"),width = 5,height = 5)
FeaturePlot(sce_Mast_Mast, features = 'MC_activation_signature1',
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()



########
sce_Mast_Mast= SetIdent(sce_Mast_Mast,value="subcelltype2") 
pdf(paste0("./",sam.name,"/MC_activation_signature_subcelltype2",".pdf"),width = 5,height = 4)
VlnPlot(sce_Mast_Mast,features = 'MC_activation_signature1')
dev.off()

##加入angiogenesis和proliferation
proliferation_score <- readxl::read_xlsx("proliferation_score.xlsx")
#转换成list
gene <- as.list(proliferation_score)
sce_Mast <- AddModuleScore(
  object = sce_Mast,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'proliferation_score'
)
pdf(paste0("./",sam.name,"/proliferation_score",".pdf"),width = 5,height = 4)
VlnPlot(sce_Mast,features = 'proliferation_score1')
dev.off()

Angiogenesis_score <- readxl::read_xlsx("Angiogenesis_score.xlsx")
#转换成list
gene <- as.list(Angiogenesis_score)
sce_Mast <- AddModuleScore(
  object = sce_Mast,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'Angiogenesis_score'
)
pdf(paste0("./",sam.name,"/Angiogenesis_score",".pdf"),width = 5,height = 4)
VlnPlot(sce_Mast,features = 'Angiogenesis_score1')
dev.off()

