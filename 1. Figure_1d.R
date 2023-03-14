# ???ߣ????㺣??
# ??λ???й?ҩ?ƴ?ѧ??????Ȼҩ???ص?ʵ???ң?????ͳ????????ҩѧ?о?????

# ???󣺷???????????????????????????֯?еĲ?????????????????ƴ?ӵĵ?ͼ????״ͼ

# ???ù???·??
workdir <- "/Users/yimang/Desktop/20221218mast_cell数据分析/泛癌表达/Mast"; setwd(workdir)

# ????R??
library(ggplot2)
library(data.table)
library(randomcoloR)
library(ggpubr)
library(GSVA)
library(clusterProfiler)
library(impute)
library(ComplexHeatmap)
source("twoclasslimma.R")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
# ??????ɫ
blue <- "#4577FF"
red <- "#C2151A"
orange <- "#E45737"
green <- "#6F8B35"
darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"

# ????ͬʱ????????????????????????
tumors <- c("BLCA","BRCA","CESC","CHOL","COAD",
            "ESCA","GBM","HNSC","KICH","KIRC",
            "KIRP","LIHC","LUAD","LUSC","PAAD",
            "PRAD","READ","STAD","THCA","UCEC")

tumors2 <- c(#"ACC",
             "BLCA","BRCA","CESC","CHOL","COAD",
             #"DLBC",
             "ESCA","GBM","HNSC","KICH","KIRC",
            "KIRP",#"LGG",
            "LIHC","LUAD","LUSC",#"MESO",
            #"OV",
            "PAAD",
            #"PCRG",
            "PRAD","READ","SARC","SKCM","STAD",
            #"TGCT",
            "THCA","THYM","UCEC"
            #"UCS",
            #"UVM"
            )

# ???ø???Ȥ?Ļ?????(TTC35??EMC2??ͬ??)
frg2<- c("TPSAB1","TPSB2","CPA3","HPGDS","MS4A2")
# ????TCGA????
# https://gdc.cancer.gov/about-data/publications/pancanatlas
rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # ????��??PanCanAtlas
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
write.table(samAnno,"simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)

#-----------#
# Figure 1B #

# ???ٶ?ȡ??????
# https://gdc.cancer.gov/about-data/publications/pancanatlas
expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1) # ????????
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),] # ?Ƴ??ظ?????
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]

comgene <- intersect(rownames(expr),frg2) # ȡ???ֱ????ף?????Ȥ?Ļ??򼯣?
expr_sub <- expr[comgene,]
colnames(expr_sub) <- substr(colnames(expr_sub),1,15)
expr_sub <- expr_sub[,!duplicated(colnames(expr_sub))]

exprTab <- ndegs <- NULL
log2fc.cutoff <- log2(1.5) # ???ò???????????ֵ (FoldChange = 1.5)
fdr.cutoff <- 0.05 # ???ò???????????ֵ (FDR = 0.05)
for (i in tumors2) {
  message("--",i,"...")
  sam <- samAnno[which(samAnno$`cancer type` == i),"simple_barcode"]
  comsam <- intersect(colnames(expr_sub), sam)
  
  tumsam <- comsam[substr(comsam,14,14) == "0"] # ????????????
  norsam <- comsam[substr(comsam,14,14) == "1"] # ????????????
  
  expr_subset <- expr_sub[,c(tumsam,norsam)]
  expr_subset[expr_subset < 0] <- 0 # ?????????????ڸ?ֵ?????㸺ֵ?Ƚ?С????ҲҪ??????????ʹ???????????????׸???????????
  expr_subset <- as.data.frame(impute.knn(as.matrix(expr_subset))$data)
  write.table(expr_subset, paste0("TCGA_",i,"_expr_subset.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  
  subt <- data.frame(condition = rep(c("tumor","normal"),c(length(tumsam),length(norsam))),
                     row.names = colnames(expr_subset),
                     stringsAsFactors = F)
  twoclasslimma(subtype  = subt, # ??????Ϣ (???뺬??һ?н?'condition')
                featmat  = expr_subset, # ?????? (???Զ??ж?????��??)
                treatVar = "tumor", # ???????顱??????????Ҫ?Ƚϵ??飩
                ctrlVar  = "normal", # ???????顱?????????Ǳ??Ƚϵ??飩
                prefix   = paste0("TCGA_",i), # ???????????ļ???ǰ׺
                overwt   = T, # ?Ƿ񸲸??Ѿ????ڵĲ????????ļ?
                sort.p   = F, # ?Ƿ?????pֵ
                verbose  = TRUE, # ?Ƿ?????????
                res.path = workdir) # ????????
  
  # ???ز????????ļ?
  res <- read.table(paste0("TCGA_",i,"_limma_test_result.tumor_vs_normal.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  upgene <- res[which(res$log2fc > log2fc.cutoff & res$padj < fdr.cutoff),] # ??ȡ?ϵ?????
  dngene <- res[which(res$log2fc < -log2fc.cutoff & res$padj < fdr.cutoff),] # ??ȡ?µ?????
  
  # ????????????????Ŀ?Բ???ͼƬ???ϲ?ע??
  if(nrow(upgene) > 0) {
    nup <- nrow(upgene)
  } else {nup <- 0}
  
  if(nrow(dngene) > 0) {
    ndn <- nrow(dngene)
  } else {ndn <- 0}
  
  exprTab <- rbind.data.frame(exprTab,
                              data.frame(gene = rownames(res),
                                         log2fc = res$log2fc,
                                         FDR = res$padj,
                                         tumor = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
  ndegs <- rbind.data.frame(ndegs,
                            data.frame(tumor = i,
                                       Group = c("UP","DOWN"),
                                       Number = c(nup,ndn),
                                       stringsAsFactors = F),
                            stringsAsFactors = F)
}
ndegs$Group <- factor(ndegs$Group, levels = c("UP","DOWN"))

# ?????ϲ?ע??
p_top <- ggplot(data = ndegs) +
  geom_bar(mapping = aes(x = tumor, y = Number, fill = Group), 
           stat = 'identity',position = 'stack') + 
  scale_fill_manual(values = c(orange,green)) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,0,0,1), "lines"))

# ?????²?????ͼ
exprTab$gene <- factor(exprTab$gene,
                       levels = rev(c("CD3D","CD8A","CD19","MS4A1","CD79A","IGHG1","CD14","CD68","SPP1","PECAM1","EPCAM","DCN","COL1A2")))
my_palette <- colorRampPalette(c(green,"white",orange), alpha=TRUE)(n=128)
p_center <- ggplot(exprTab, aes(x=tumor,y=gene)) +
  geom_point(aes(size=-log10(FDR),color=log2fc)) +
  scale_color_gradientn('log2(FC)', 
                        colors=my_palette) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, color = rep(c(red,blue),c(14,10))),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0,0,1,1), "lines")) 

# ?Ų?ͼƬ
ggarrange(p_top,
          p_center, 
          nrow = 2, ncol = 1,
          align = "v",
          heights = c(2,6),
          common.legend = F)
ggsave("Figure 1B differential expression of interested genes in pancancer_33.pdf", width = 8,height =13)
rm(expr); gc()

