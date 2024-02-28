# TCGA-LUAD数据下载
setwd('TCGA-LUAD')
setwd('TCGAdata')
library(tidyverse)
# install.packages('BiocManager')
library(BiocManager)
# 安装TCGAbiolinks
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("remotes")
# BiocManager::install("ExperimentHub")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

library(TCGAbiolinks)
cancer_type <- "TCGA-LUAD"#肿瘤类型
# TCGA肿瘤缩写：https://www.jianshu.com/p/3c0f74e85825
# 准备下载
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts")
GDCdownload(expquery,directory = "GDCdata")
#准备为标准需要的数据格式
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
save(expquery2,file="luad.gdc_2022.rda")#保存为rda格式

load("luad.gdc_2022.rda")
# 下载humo_sapiens-https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
if(!require("rtracklayer")) BiocManager::install("rtracklayer")
homo_anno <- rtracklayer::import('Homo_sapiens.GRCh38.111.chr.gtf.gz')
#转换为数据框
homo_anno <- as.data.frame(homo_anno)
# dim(gtf)
# table(gtf$gene_biotype)
gtf <- homo_anno[,c("gene_id",'gene_name','gene_biotype')]
colnames(gtf) <- c('ENSEMBL','symbol','type')
save(gtf,file = "gene_annotation_2022.rda")
duplicated(gtf$ENSEMBL)
gtf2 <- gtf[!duplicated(gtf$ENSEMBL),]
save(gtf2,file = "gene_annotation_duplicated.rda")

load("gene_annotation_2022.rda")#导入gene注释文件
table(gtf$type)#table分组计数
load("gene_annotation_duplicated.rda")
table(gtf2$type)#table分组计数

# 基因名称ENSEMBL symbol
# 提取counts只能进行差异分析；tpms基因表达谱(只针对TCGA数据库)
counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES
 
# counts<-separate(counts,gene_id,into = c("gene_id"),sep="\\.") 
# ID转换

# counts <- counts %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% inner_join("ENSEMBL") %>% .[!duplicated(.$symbol),]
counts <- as.data.frame(counts)
counts <- rownames_to_column(counts,var='ENSEMBL')
counts<-separate(counts,ENSEMBL,into = c("ENSEMBL"),sep="\\.") 
counts <- inner_join(counts,gtf2,'ENSEMBL')
counts <- counts[!duplicated(counts$symbol),]
# 判断两个df是否相同
# identical(a,b)
homo_gtf <- gtf2
save(homo_gtf,file = "homo_gtf.rda")
save(counts,file = "final_counts.rda")


# # mRNA的counts矩阵
# expr_counts_mrna <- assay(se_mrna,"unstranded")

# mRNA的tpm矩阵
# expr_tpm_mrna <- assay(se_mrna,"tpm_unstrand")
expr_tpm_mrna <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(expr_tpm_mrna) <- expquery2@colData@rownames
rownames(expr_tpm_mrna) <- expquery2@rowRanges@ranges@NAMES
tpms <- as.data.frame(expr_tpm_mrna)
tpms <- rownames_to_column(tpms,var='ENSEMBL')
tpms<-separate(tpms,ENSEMBL,into = c("ENSEMBL"),sep="\\.") 
tpms <- inner_join(tpms,gtf2,'ENSEMBL')
tpms <- tpms[!duplicated(tpms$symbol),]
save(tpms,file = "final_tpms.rda")
# # mRNA的fpkm矩阵
# expr_fpkm_mrna <- assay(se_mrna,"fpkm_unstrand")
# 
# # lncRNA的counts矩阵
# expr_counts_lnc <- assay(se_lnc,"unstranded")
# 
# # lncRNA的tpm矩阵
# expr_tpm_lnc <- assay(se_lnc,"tpm_unstrand")
# 
# # lncRNA的fpkm矩阵
# expr_fpkm_lnc <- assay(se_lnc,"fpkm_unstrand")

# 计算患者免疫评分与肿瘤纯度--肿瘤包括基质组分，免疫组分，肿瘤组分
















