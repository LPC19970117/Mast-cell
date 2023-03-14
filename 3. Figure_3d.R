#cell 火山图
setwd("E:/生物信息学/转录组火山图")
df <- read.csv("DEGs_trans.csv",row.names = 1) #导入数据，第一列作为行名
###################图1

plot(df$avg_log2FC, -log10(df$p_val), 
     col="#00000033", pch=19,cex=2.5,
     xlab=paste("log2 (fold change)"),
     ylab="-log10 (p_value)")+title("DEGs")

#筛选上下调
up <- subset(df, df$p_val < 0.01 & df$avg_log2FC > 0.05)
down <- subset(df, df$p_val < 0.01 & df$avg_log2FC < as.numeric(-0.05))
#绘制上下调
points(up$avg_log2FC, -log10(up$p_val), col=1, bg = brewer.pal(9, "YlOrRd")[6], pch=21, cex=2.5)
points(down$avg_log2FC, -log10(down$p_val), col = 1, bg = brewer.pal(11,"RdBu")[9], pch = 21,cex=2.5)
#加上阈值线
abline(h=-log10(0.01),v=c(-5,5),lty=2,lwd=1)
text(-0.4,-log10(3.260000e-64),"OGFOD2")

##################################################################图2
df <- read.csv("DEGs_trans.csv",row.names = 1) #导入数据，第一列作为行名
library(ggpubr)
library(ggthemes)
df$group = 'ns'#添加一组
df$group[which((df$p_val <0.05)&(df$avg_log2FC>0.05))]='Up'#上调基因
df$group[which((df$p_val <0.05)&(df$avg_log2FC< -0.05))]='Down'#下调基因
df$label ='' #标记，用于后续显示基因名称
df <- df[order(df$p_val),]#排序
up <- head(rownames(df)[which(df$group =='Up')],20)#挑选up10
down <- head(rownames(df)[which(df$group =='Down')],20)#挑选down10
top10 <- c(as.character(up), as.character(down))
df$label[match(top10, rownames(df))] <- top10 #将up，down基因名添加在label这一列

df$p_val <-  -log10(df$p_val)#将P转化为-logp
pdf(paste0("./",sam.name,"/火山图2.3.pdf"),width = 8,height = 7)
ggscatter(df,
          x= 'avg_log2FC',
          y= 'p_val',
          color = 'group',
          palette = c("#2f5688", "#BBBBBB","#CC0000"),
          size = 2,#点大小
          label = df$label,
          font.label = 12,#label字体大小
          repel = T,
          xlab="log2 (fold change)",
          ylab="-log10 (p_value)")+ theme_base()+
  geom_hline(yintercept = 1.30, linetype="dashed")+  #添加y轴线
  geom_vline(xintercept= c(-0.05,0.05), linetype="dashed")#添加x轴线
dev.off()
##############################################################图3
df <- read.csv("DEGs_trans.csv",row.names = 1) #导入数据，第一列作为行名
library(ggplot2)

df$group = 'ns'#添加一组
df$group[which((df$p_val <0.05)&(df$avg_log2FC>0))]='Up'#上调基因
df$group[which((df$p_val <0.05)&(df$avg_log2FC<-0))]='Down'#下调基因
df$label ='' #标记，用于后续显示基因名称
df <- df[order(df$p_val),]#排序
up <- head(rownames(df)[which(df$group =='Up')],10)#挑选up10
down <- head(rownames(df)[which(df$group =='Down')],10)#挑选down10
top10 <- c(as.character(up), as.character(down))
df$label[match(top10, rownames(df))] <- top10 #将up，down基因名添加在label这一列

df$p_val <-  -log10(df$p_val)#将P转化为-logp

ggplot(df,aes(avg_log2FC, p_val))+
  geom_hline(yintercept = 1.30, linetype="dashed", color = "#999999")+  #添加y轴线
  geom_vline(xintercept= c(-0.6,0.6), linetype="dashed", color = "#999999")+
  geom_point(aes(size=p_val, color= p_val))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(1,3))+
  theme_bw()+
  theme(panel.grid = element_blank())

library(ggrepel)
pdf(paste0("./","火山图3.pdf"),width = 7,height = 8)
ggplot(df,aes(avg_log2FC, p_val))+
  geom_hline(yintercept = 1.30, linetype="dashed", color = "#999999")+  #添加y轴线
  geom_vline(xintercept= c(-0.05,0.05), linetype="dashed", color = "#999999")+
  geom_point(aes(size=p_val, color= p_val))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(2,4))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.4),
        legend.justification = c(0,1))+
  guides(col = guide_colourbar(title = "-Log10_P-value"),
         size = "none")+
  geom_text_repel(aes(label=label, color = p_val), size = 3, vjust = 1.5, hjust=1)+
  xlab("log2 (fold change)")+
  ylab("-log10 (p_value)")
dev.off()
