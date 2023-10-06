
BiocManager::install("monocle")
# monocle2
rm(list=ls())
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)


load(file="sce_Mast2.RData")
sce<-sce_Mast
#4.提取Mast_cell亚组#4. Extract Mast_cell subgroup
sce<-subset(sce,subcelltype2 %in% c("Activated_MC","Proliferating_MC","Resting_MC"))
table(sce_Mast@meta.data$celltype)

exp.matrix<-as(as.matrix(sce@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann)<-rownames(exp.matrix)
exp_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-sce@meta.data
rownames(sample_ann)<-colnames(exp.matrix)
exp_pd<-new("AnnotatedDataFrame", data =sample_ann)

#生成monocle对象#Generate monocle object
exp.monocle<-newCellDataSet(exp.matrix,phenoData =exp_pd,featureData =exp_fd,expressionFamily=negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#calculate sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)


#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因#Calculate differentially expressed genes based on seurat cluster and select genes used to construct pseudo-time series trajectories
table(sce$subcelltype)
if(F){
  ex1 <- c(cc.genes$s.genes, cc.genes$g2m.genes)
  ex2 <- grep('^MT-', rownames(sce), value = T)
  ex3 <- grep('^RP[LS]', rownames(sce), value = T)
  ex <- unique(c(ex1, ex2, ex3))
  #确定需要加入的基因#Determine the genes that need to be added
  include <- c("TPSAB1","TPSB2","CPA3","HPGDS","MS4A2","KIT","CMA1")
}
#order.genes <-  VariableFeatures()
#order.genes <- order.genes[!order.genes %in% ex]
#exp.monocle <- setOrderingFilter(exp.monocle, order.genes)
#plot_ordering_genes(exp.monocle)

if(F){
  disp_table <- dispersionTable(exp.monocle)
  order.genes <- subset(disp_table, mean_expression >= 0.01 & 
                          dispersion_empirical >= 1 * dispersion_fit) %>% 
    pull(gene_id) %>% as.character()
  #order.genes <- order.genes[!order.genes %in% ex]
  exp.monocle <- setOrderingFilter(exp.monocle, order.genes)
  plot_ordering_genes(exp.monocle)
}
#DDRTree方法降维并构建拟时序#DDRTree method reduces dimensionality and constructs pseudo-time series
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
save(exp.monocle, file = "exp.monocle.Rdata")
#画图#making plot
plot_cell_trajectory(exp.monocle,color_by = "subcelltype")
ggsave("subcelltype.pdf",device = "pdf",width = 18,height = 19,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "subcelltype2")
ggsave("subcelltype2.pdf",device = "pdf",width = 18,height = 19,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "tissue")
ggsave("tissue.pdf",device = "pdf",width = 18,height = 19,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "State")
ggsave("State.pdf",device = "pdf",width = 18,height = 19,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "Pseudotime")
ggsave("Pseudotime.pdf",device = "pdf",width = 18,height = 19,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "subcelltype")+facet_wrap(~subcelltype,nrow=2)
ggsave("subcelltype.pdf",device = "pdf",width = 21,height = 19,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "subcelltype2")+facet_wrap(~subcelltype2,nrow=2)
ggsave("subcelltype2.pdf",device = "pdf",width = 21,height = 19,units = c("cm"))

s.genes <- c("TPSAB1","TPSB2","CPA3","HPGDS","MS4A2","KIT","CMA1")
plot_genes_violin(exp.monocle[s.genes,], grouping = "State", color_by = "State")
plot5 <-plot_genes_in_pseudotime(exp.monocle[s.genes,], color_by = "Pseudotime")
pdf("plot5_Pseudotime.pdf",width = 10,height = 5)
print(plot5)
dev.off()

plot5 <-plot_genes_in_pseudotime(exp.monocle[s.genes,], color_by = "tissue")
pdf("plot5_tissue.pdf",width = 5,height = 5)
print(plot5)
dev.off()

s.genes <- c("KIT","TPSAB1","CPA3","CMA1")
plot5 <-plot_genes_in_pseudotime(exp.monocle[s.genes,], color_by = "subcelltype2")
pdf("plot5_subcelltype2_2.pdf",width = 5,height = 3)
print(plot5)
dev.off()
