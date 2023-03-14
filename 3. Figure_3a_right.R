library(Seurat)
library(ggplot2)
library(dplyr)
# setwd("D:/KS项目/公众号文章/堆叠柱状图显示比例")

####不同组织tissue中不同subcelltype2比例
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
