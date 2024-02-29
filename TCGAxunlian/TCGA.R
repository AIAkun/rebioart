# 1+1#1+1代表普通的运算
#install.packages("tidyverse")
#也可以通过右下角的操作栏通过点击来安装R包
#chooseBioCmirror()
#a<-6
#a
#以上是试验赋值的用法，还可以用=来表示，不过==代表一种逻辑判断
#setwd("DAY1")#此代码用来设置工作目录（临时工作目录，退出后回到初始目录），以免在使用过程中将储存文件弄混
#dir.create()可以用来直接创建目录，不过麻烦
#getwd()#查看当前目录
#setwd("DAY1.1")
#ctrl+F可以用来查找代码框中的代码，ctrl+Z是撤销一次的含义
#read.table()#读取文本文档
#write.CSV()#是输出为表格
write.table(counts01A,"counts01A.txt",sep="\t",row.names=T,col.names=NA,quote=F)#是输出文本文档#对文件的读取或输出可以直接点击右上角的读取来导入#通常txt格式会更不容易出错
#t()是用于行列反转的函数
#当一个可以判断为数据框的格式但R语言并没有将其判断为数据框时，可以通过函数as.data.frame()来改变

#class()#用来判断数据框类型，若是数据框内的某一列数据通过输入$符号来确定
#[,]代表对数据框进行提取。其中：表示多少到多少。并且，前是行后是列
#%>%表示传导符号

#install.packages("tidyverse")
#library(tidyverse)
#a<-c("a","b","a","b","c")
#duplicated(a)#用来判断数据里的元素是否重复
#a<-!duplicated(a)#感叹号是指将数据集内判断结果反过来，[]中括号可以把判断中是TRUE的元素提取出来
#a[a]#以上操作中将判断结果赋予了a故导致原本数据的丢失
#a<-c("a","b","a","b","c")
#duplicated(a)
#a<-a[!duplicated(a)]
#a
#这才是正确操作结果

#inner_join()#是将两数据框中数据相同的进行合并的函数，包括（表一，表二，by='合并标准'）
#left_join()#则是代表以函数内左的数据框为标准进行合并


###接下来进行数据的下载实战#####
##建立好文件##
setwd("TCGA-LUAD")
setwd('TCGAdata')
library(tidyverse)
# install.packages('BiocManager')
# library(BiocManager)
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("remotes")
# BiocManager::install("ExperimentHub")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
cancer_type<-"TCGA-LUAD"
expquery<-GDCquery(project=cancer_type,
                   data.category="Transcriptome Profiling",
                   data.type="Gene Expression Quantification",
                   workflow.type="STAR - Counts")
GDCdownload(expquery,directory="GDCdata")
expquery2<-GDCprepare(expquery,directory="GDCdata",summarizedExperiment = T)
save(expquery2,file="luad.gdc_2023.rad")#对数据进行保存，rad格式是属于R的格式

setwd('TCGA-LUAD')
setwd("TCGAdata")
load("luad.gdc_2022.rda")

library(tidyverse)
counts<-expquery2@assays@data@listData[["unstranded"]]
colnames(counts)<-expquery2@colData@rownames
rownames(counts)<-expquery2@rowRanges@ranges@NAMES
counts<-as.data.frame(counts)
counts<-rownames_to_column(counts,var = "ENSEMBL")
counts<-separate(counts,ENSEMBL,into = c("ENSEMBL"),sep="\\.") 
counts<-inner_join(counts,homo_gtf,"ENSEMBL")
counts<-counts[!duplicated(counts$symbol),]


rownames(counts)<-NULL#去除行名
counts <- counts %>% drop_na(symbol)
counts<-column_to_rownames(counts,var = "symbol")
table(counts$type)#数一下有多少的type
counts<-counts[counts$type == "protein_coding",]
counts<-counts[,-c(1,ncol(counts))]
colnames(counts)<-substring(colnames(counts),1,16)#将列名中1~16字符提取出来 
counts<-counts[,!duplicated(colnames(counts))]

table(substring(colnames(counts),14,16))#通常01A代表肿瘤样本，11A代表正常样本
counts01A<-counts[,substring(colnames(counts),14,16)=="01A"]
counts11A<-counts[,substring(colnames(counts),14,16)=="11A"]#将样本中01A和11A的样本提出来做分析

write.table(counts01A,"counts01A.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(counts11A,"counts11A.txt",sep="\t",row.names=T,col.names=NA,quote=F)



###接下来提取tpms###前面提的counts是用来做差异分析的，而tpms是用来做后面分析的

tpms1<-expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms1)<-expquery2@colData@rownames
rownames(tpms1)<-expquery2@rowRanges@ranges@NAMES
tpms1<-tpms1%>%
  as.data.frame()%>%
  rownames_to_column("ENSEMBL")%>%separate(ENSEMBL,into = c("ENSEMBL"),sep="\\.") %>% 
  inner_join(homo_gtf,"ENSEMBL")%>%
  .[!duplicated(.$symbol),]

rownames(tpms1)<-NULL
tpms1 <- tpms1 %>% drop_na(symbol)
tpms1<-tpms1 %>%column_to_rownames("symbol")
table(tpms1$type)#查看分类计数

tpms1<-tpms1[tpms1$type == "protein_coding",]
tpms1<-tpms1[,-c(1,ncol(tpms1))]
colnames(tpms1)<-substring(colnames(tpms1),1,16)#将列名中1~16字符提取出来 
tpms1<-tpms1[,!duplicated(colnames(tpms1))]
table(substring(colnames(tpms1),14,16))

tpms01A<-tpms1[,substring(colnames(tpms1),14,16)=="01A"]
tpms11A<-tpms1[,substring(colnames(tpms1),14,16)=="11A"]


identical(rownames(counts01A),rownames(tpms01A))
identical(rownames(counts11A),rownames(tpms11A))
identical(colnames(counts01A),colnames(tpms01A))
identical(colnames(counts11A),colnames(tpms11A))#验证几个样本是否一致

write.table(tpms01A,"tpms01A.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(tpms11A,"tpms11A.txt",sep="\t",row.names=T,col.names=NA,quote=F)


counts<-cbind(counts01A,counts11A)#cbind是以列colnames来合并，rbind是以行rownames来合并
tpms<-cbind(tpms01A,tpms11A)
write.table(tpms,"tpms.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(counts,"counts.txt",sep="\t",row.names=T,col.names=NA,quote=F)


range(tpms)#查看数据范围
range(tpms01A)
range(tpms11A)
tpms_log2<-log2(tpms+1)
range(tpms_log2)
tpms01A_log2<-log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2<-log2(tpms11A+1)
range(tpms11A_log2)
write.table(tpms_log2,"tpms_log2.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(tpms01A_log2,"tpms01A_log2.txt",sep="\t",row.names=T,col.names=NA,quote=F)
write.table(tpms11A_log2,"tpms11A_log2.txt",sep="\t",row.names=T,col.names=NA,quote=F)



###计算肿瘤患者的免疫评分###（基质，免疫，肿瘤）肿瘤纯度=肿瘤所占比例
setwd("TCGA-LUAD")
setwd("ESTIMATE")
#rfoge<-"http://r-forge.r-project.org"#RFORGE就是R语言自己的R包平台
#install.packages("estimate",repos=rfoge,dependencies=TRUE)
library(estimate)#estimate是用来计算免疫评分的
library(tidyverse)
exp<-read.table("tpms01A_log2.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F
                ,header=T)

filterCommonGenes(input.f="tpms01A_log2.txt",
                  output.f="tpms01A_log2.gct",
                  id="GeneSymbol")
estimateScore("tpms01A_log2.gct",
              "tpms01A_log2_estimate_score.txt",
              platform="affymetrix")#affymetrix是一个基因芯片平台，相比于其它平台结果更好

ESTIMATE_result<-read.table("tpms01A_log2_estimate_score.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F
                            ,header=T)
ESTIMATE_result<-ESTIMATE_result[,-1]
# column_to_rownames(ESTIMATE_result,var="NAME")
colnames(ESTIMATE_result)<-ESTIMATE_result[1,]
ESTIMATE_result<-t(ESTIMATE_result[-1,])#通常行列转换后会变得不是数据框形式
ESTIMATE_result2<-as.data.frame(ESTIMATE_result)#可以直接用ESTIMATE_result<-as.data.frame(t(ESTIMATE_result[-1,]))

rownames(ESTIMATE_result)<-colnames(exp)
write.table(ESTIMATE_result,"ESTIMATE_result.txt",sep="\t",row.names=T,col.names=NA,quote=F)


##进行生存分析##
# xena官网：https://xenabrowser.net
setwd("Survival_data")
library(tidyverse)
survival_data<-read.table("OS.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F
                          ,header=T)
survival_data<-survival_data[,2:3]
survival<-rownames_to_column(survival_data,"sample")
survival$name<-paste0(survival$sample,"A")#当粘贴到的数据框中没有name这列时，R会自动生成新的一列
#paste0表示粘贴连接一字符“”
table(substring(survival$name,14,16))
rownames(survival)<-survival$name
survival<-survival[,2:3]

#接下将生存信息与基因表达谱合并起来#
tpms01A_log2<-read.table("tpms01A_log2.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F
                         ,header=T)
a<-intersect(colnames(tpms01A_log2),rownames(survival))#对生存数据和基因表达谱
#~的sample名取交集，主要是为了得到一个额外的sample顺序，用来接下来将两者sample对齐
table(substring(a,14,16))
exp_01A<-tpms01A_log2[,a]
surv_01A<-survival[a,]#以顺序或以集合a作为模板来将两个数据中以行或列进行数据提取
exp_01A<-exp_01A%>%t()%>%as.data.frame()
identical(rownames(exp_01A),rownames(surv_01A))
exp_surv_01A<-cbind(surv_01A,exp_01A)
#将合并后的生存信息和基因表达谱合并的文件保存
write.table(exp_surv_01A,"exp_surv_01A.txt",sep="\t",row.names=T,col.names=NA,quote=F)


#合并生存信息与免疫评分数据（ESTIMATE）
ESTIMATE_result<-read.table("ESTIMATE_result.txt",sep="\t",row.names=1,check.names=F,stringsAsFactors = F
                            ,header=T)
identical(rownames(ESTIMATE_result),rownames(surv_01A))
ESTIMATE_result_surv_01A<-cbind(surv_01A,ESTIMATE_result)
write.table(ESTIMATE_result_surv_01A,"ESTIMATE_result_surv_01A.txt",sep="\t",row.names=T,col.names=NA,quote=F)