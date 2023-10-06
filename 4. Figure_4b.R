options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

library(ggplot2)
library(data.table)
library(survival)
library(ComplexHeatmap)
library(forestplot)
library(survminer)
library(circlize)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

#读取肿瘤注释文件#Read tumor annotation file
rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
write.table(samAnno,"output_simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)

# 读取生存数据# Read survival data
surv <- read.delim("Survival_SupplementalTable_S1_20171025_xena_sp", sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 快速读取表达谱数据并做数据预处理# Quickly read expression profile data and do data preprocessing
#expr <- fread("222.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),]
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
expr[expr < 0] <- 0 # 对于这份泛癌数据，将略小于0的数值拉到0，否则不能取log（其他途径下载的泛癌数据可能不需要此操作）# For this pan-cancer data, pull the value slightly less than 0 to 0, otherwise the log cannot be taken (pan-cancer data downloaded through other channels may not require this operation)
colnames(expr) <- substr(colnames(expr),1,15)
gc()

# 设置感兴趣的基因# Set the genes of interest
geneOfInterest <- "TPSAB1"#or other four MC signature genes
if(!is.element(geneOfInterest, rownames(expr))) {
  warning("The gene ", geneOfInterest," cannot be found!")
} else {
  message("The gene ", geneOfInterest," can be matched!")
}

# 确定肿瘤样本以及对应肿瘤类型# Determine the tumor sample and corresponding tumor type
sam <- samAnno[which(samAnno$`cancer type` != "LAML"),"simple_barcode"] # 去掉白血病样本# Remove leukemia samples
comsam <- intersect(intersect(colnames(expr), sam), rownames(surv)) # 得到与表达谱以及生存的共有样本# Get shared samples with expression profiles and survival
tumsam <- comsam[substr(comsam,14,14) == "0"] # 仅提取肿瘤样本#Extract tumor samples only
tumAnno <- samAnno[which(samAnno$simple_barcode %in% tumsam),] # 获取这些肿瘤样本的注释信息# Get the annotation information of these tumor samples
tumAnno <- tumAnno[order(tumAnno$`cancer type`),] # 根据肿瘤类型排序# Sort by tumor type
tumors <- unique(tumAnno$`cancer type`) # 得到32个肿瘤# Get 32 tumors

# 合并表达和生存数据# Merge expression and survival data
exprSurv <- cbind.data.frame(expr = log2(as.numeric(expr[geneOfInterest,comsam]) + 1),
                             surv[comsam,c("OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")])
write.table(exprSurv, "output_combined dataframe with both gene expression and survival information.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


outTab.cox <- NULL
for(i in tumors) {
  sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"]
  exprSurvSub <- exprSurv[sam,]
  
  ## OS
  coxres <- summary(coxph(Surv(OS.time, OS) ~ expr, data = exprSurvSub))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(tumor = i, 
                                            event = "OS", 
                                            beta = coxres$coefficients[1,1], 
                                            hr = coxres$coefficients[1,2], 
                                            lower = coxres$conf.int[1,3], 
                                            upper = coxres$conf.int[1,4], 
                                            p = coxres$coefficients[1,5], 
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
  
  ## DSS
  coxres <- summary(coxph(Surv(DSS.time, DSS) ~ expr, data = exprSurvSub))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(tumor = i,
                                            event = "DSS",
                                            beta = coxres$coefficients[1,1],
                                            hr = coxres$coefficients[1,2],
                                            lower = coxres$conf.int[1,3],
                                            upper = coxres$conf.int[1,4],
                                            p = coxres$coefficients[1,5],
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
  
  ## DFI
  if(i %in% c("SKCM","THYM","UVM","GBM")) { # 前三个肿瘤没有对应的DFI数据，GBM只有3个样本有对应的DFI，也去除（和原文不符）# The first three tumors do not have corresponding DFI data. Only 3 samples of GBM have corresponding DFI data and are also removed (inconsistent with the original text)
    outTab.cox <- rbind.data.frame(outTab.cox,
                                   data.frame(tumor = i,
                                              event = "DFI",
                                              beta = NA,
                                              hr = NA,
                                              lower = NA,
                                              upper = NA,
                                              p = NA,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  } else {
    coxres <- summary(coxph(Surv(DFI.time, DFI) ~ expr, data = exprSurvSub))
    outTab.cox <- rbind.data.frame(outTab.cox,
                                   data.frame(tumor = i,
                                              event = "DFI",
                                              beta = coxres$coefficients[1,1],
                                              hr = coxres$coefficients[1,2],
                                              lower = coxres$conf.int[1,3],
                                              upper = coxres$conf.int[1,4],
                                              p = coxres$coefficients[1,5],
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  }
  
  ## PFI
  coxres <- summary(coxph(Surv(PFI.time, PFI) ~ expr, data = exprSurvSub))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(tumor = i,
                                            event = "PFI",
                                            beta = coxres$coefficients[1,1],
                                            hr = coxres$coefficients[1,2],
                                            lower = coxres$conf.int[1,4],
                                            upper = coxres$conf.int[1,4],
                                            p = coxres$coefficients[1,5],
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
}
write.table(outTab.cox, file = "output_summary of cox result in pancancer.txt",sep = "\r",row.names = F,col.names = T,quote = F)


minprop <- 0.2
outTab.km <- NULL
for(i in tumors) {
  sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"]
  exprSurvSub <- exprSurv[sam,]
  
  ## OS
  bestcut <- surv_cutpoint(exprSurvSub, 
                           time = "OS.time", 
                           event = "OS", 
                           variables = "expr", 
                           minprop = minprop) #默认组内sample不能低于20%#The sample in the default group cannot be less than 20%
  cutoff <- bestcut$cutpoint[1,1]
  exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("Low","High"))
  fitd <- survdiff(Surv(OS.time, OS) ~ group, data=exprSurvSub, na.action=na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  HR <- (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
  outTab.km <- rbind.data.frame(outTab.km,
                                data.frame(tumor = i, 
                                           event = "OS", 
                                           hr = HR, # Hazard ratio
                                           lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])),
                                           upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])),
                                           p = p.val,
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
  
  ## DSS
  bestcut <- surv_cutpoint(exprSurvSub, 
                           time = "DSS.time", 
                           event = "DSS", 
                           variables = "expr", 
                           minprop = minprop) 
  cutoff <- bestcut$cutpoint[1,1]
  exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("High","Low"))
  fitd <- survdiff(Surv(DSS.time, DSS) ~ group, data=exprSurvSub, na.action=na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  outTab.km <- rbind.data.frame(outTab.km,
                                data.frame(tumor = i, 
                                           event = "DSS",
                                           hr = (fitd$obs[1]/fitd$exp[1])/(fitd$obs[2]/fitd$exp[2]), # Hazard ratio 我把1和2对调后结果才看起来对# Hazard ratio I swapped 1 and 2 and the result looked right.
                                           lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), 
                                           upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])),
                                           p = p.val, # p值
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
  
  ## DFI
  if(i %in% c("SKCM","THYM","UVM","GBM")) { # 前三个肿瘤没有对应的DFI数据，GBM只有3个样本有对应的DFI，也去除（和原文不符）# The first three tumors do not have corresponding DFI data. Only 3 samples of GBM have corresponding DFI data and are also removed (inconsistent with the original text)
    outTab.km <- rbind.data.frame(outTab.km,
                                  data.frame(tumor = i,
                                             event = "DFI",
                                             hr = NA,
                                             lower = NA,
                                             upper = NA,
                                             p = NA,
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
  } else {
    bestcut <- surv_cutpoint(exprSurvSub, 
                             time = "DFI.time", 
                             event = "DFI", 
                             variables = "expr", 
                             minprop = minprop) 
    cutoff <- bestcut$cutpoint[1,1]
    exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("High","Low"))
    fitd <- survdiff(Surv(DFI.time, DFI) ~ group, data=exprSurvSub, na.action=na.exclude)
    p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
    outTab.km <- rbind.data.frame(outTab.km,
                                  data.frame(tumor = i, 
                                             event = "DFI",
                                             hr = (fitd$obs[1]/fitd$exp[1])/(fitd$obs[2]/fitd$exp[2]),# Hazard ratio
                                             lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), 
                                             upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), 
                                             p = p.val, 
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
  }
  
  ## PFI
  bestcut <- surv_cutpoint(exprSurvSub, 
                           time = "PFI.time", 
                           event = "PFI", 
                           variables = "expr", 
                           minprop = minprop) 
  cutoff <- bestcut$cutpoint[1,1]
  exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("High","Low"))
  fitd <- survdiff(Surv(PFI.time, PFI) ~ group, data=exprSurvSub, na.action=na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  outTab.km <- rbind.data.frame(outTab.km,
                                data.frame(tumor = i, 
                                           event = "PFI",
                                           hr = (fitd$obs[1]/fitd$exp[1])/(fitd$obs[2]/fitd$exp[2]),# Hazard ratio
                                           lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), #
                                           upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 
                                           p = p.val, #
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
write.table(outTab.km, file = "output_summary of km result in pancancer.txt",sep = "\r",row.names = F,col.names = T,quote = F)


# 设置颜色# Set color
blue   <- "#A0CEE3"
yellow <- "#EFEFBE"
sea    <- "#37BCDF"
green  <- "#71BC5D"
cherry <- "#E5588C"
red    <- "#E46A6B"
purple <- "#8959A3"

# 制作热图绘制数据# Make heat map drawing data
hmInput <- NULL
for (i in tumors) {
  
  cox.res <- outTab.cox[which(outTab.cox$tumor == i),]
  km.res <- outTab.km[which(outTab.km$tumor == i),]
  
  cox.res$dirct <- ifelse(cox.res$hr > 1 & cox.res$p < 0.05, "Risky",
                          ifelse(cox.res$hr < 1 & cox.res$p < 0.05, "Protective","Nonsense"))
  km.res$dirct <- ifelse(km.res$hr > 1 & km.res$p < 0.05, "Risky",
                         ifelse(km.res$hr < 1 & km.res$p < 0.05, "Protective","Nonsense"))
  hmInput <- rbind.data.frame(hmInput,
                              data.frame(OS.cox = cox.res[1,"dirct"],
                                         OS.km = km.res[1,"dirct"],
                                         DSS.cox = cox.res[2,"dirct"],
                                         DSS.km = km.res[2,"dirct"],
                                         DFI.cox = cox.res[3,"dirct"],
                                         DFI.km = km.res[3,"dirct"],
                                         PFI.cox = cox.res[4,"dirct"],
                                         PFI.km = km.res[4,"dirct"],
                                         row.names = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
}
hmInput[is.na(hmInput)] <- "N/A"
write.table(hmInput, file = "output_summary of gene prognositicationin pancancer.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 注意：结果在某些肿瘤上有较大出入，可能的原因有：
# 1. 泛癌表达谱使用不一致
# 2. 最佳生存曲线cutoff的minprop选择不一致
# 3. 对KM分析中HR的计算不一致，原作者可能依然采用Cox回归估计二分类下的HR
# Note: The results vary greatly for some tumors. Possible reasons are:
# 1. Inconsistent use of pan-cancer expression profiles
# 2. The minprop selection of the optimal survival curve cutoff is inconsistent
# 3. The calculation of HR in KM analysis is inconsistent. The original author may still use Cox regression to estimate HR under two categories.

# 修改热图数据#Modify heat map data
indata <- hmInput
indata[indata == "Risky"] <- 0 # 为了使得risky出现在图例的顶部，所以设置一个最小值0To make risky appear at the top of the legend, set a minimum value of 0
indata[indata == "Protective"] <- 1 # 根据原文，依次设置为1# According to the original text, set to 1 in turn
indata[indata == "N/A"] <- 2 # 根据原文，依次设置为2
indata[indata == "Nonsense"] <- 3 # 根据原文，依次设置为3

# 构建列注释#Build column annotations
annCol <- data.frame("Method" = rep(c("Cox","Log-rank"),4),
                     "Survival Type" = rep(c("OS","DSS","DFI","PFI"), each = 2),
                     check.names = F, 
                     row.names = colnames(indata))
annColors <- list("Method" = c("Cox" = blue, "Log-rank" = yellow),
                  "Survival Type" = c("OS" = sea, "DSS" = green, "DFI" = cherry, "PFI" = purple))
col_fun = circlize::colorRamp2(c(0, 1, 2, 3), c(red,green,"grey50","white")) 

# 绘制热图#Draw a heat map
hm <- Heatmap(matrix = as.matrix(indata),
              border = "black", 
              rect_gp = gpar(col = "black"), 
              name = "Prognostic role",
              cluster_rows = F, 
              cluster_columns = F,
              col = c(red,green,"grey50","white"), 
              show_row_names = T, 
              show_column_names = F, 
              row_names_side = "left",
              top_annotation = HeatmapAnnotation(df = annCol,
                                                 col = annColors,
                                                 gp = gpar(col = "black"),
                                                 border = TRUE),
              width = grid::unit(8, "cm"),
              height = grid::unit(15, "cm"),
              heatmap_legend_param = list(at = c(0, 1, 2, 3), 
                                          legend_gp = grid::gpar(fill = col_fun(c(0,1,2,3))),
                                          labels = c("Risky", "Protective","N/A","Nonsense")))
pdf("prognostic heatmap.pdf", width = 10,height = 10)
draw(hm)
invisible(dev.off())



##################################################################
fpInput <- outTab.cox[which(outTab.cox$event == "OS"),]
hrtable <- fpInput[,c("tumor","event","beta","hr","lower","upper","p")]

tabletext <- cbind(c("Cancers",hrtable$tumor),
                   c("p value",ifelse(round(as.numeric(hrtable$p),3) < 0.001,"<0.001",round(as.numeric(hrtable$p),3))),
                   c("HR (95L-95H)",paste0(round(as.numeric(hrtable$hr),3), " (",
                                           round(as.numeric(hrtable$lower),3),"-",
                                           round(as.numeric(hrtable$upper),3),")")))

pdf("forestplot of os risk table in pancancer.pdf", width = 8, height = 8)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$hr)),# HR
           lower=c(NA,as.numeric(hrtable$lower)), #
           upper=c(NA,as.numeric(hrtable$upper)),# 
           graph.pos = 4,
           graphwidth = unit(.3,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="grey50", lines="grey50", zero = "black"),
           boxsize=c(NA,ifelse(as.numeric(hrtable$p) < 0.05,0.8,0.4)),
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=1,
           lwd.zero=2,
           xticks = c(0,1,2,4,8,13),
           lwd.xaxis=2,
           xlab="Hazard ratio",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "34" = gpar(lwd=2, col="black")),#最后一行底部加黑线???""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.2)),
           lineheight = unit(.55,"cm"),
           colgap = unit(0.4,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())

sessionInfo()

