library(tidyverse)#加载Error
# 数据读取/保存
content <- read.table('content.txt',sep = '\t',row.names = 1,check.names = F)
write.table(content,'content_save.txt',sep = '\t',row.names = T,col.names = NA)
write.csv(content,file='')
content_csv <- read.csv('',row.names = 1)
# 数据库行列反转
t_content <- t(content)
# as.data.frame
df <- as.data.frame(content)
# class判断数据类型：数字，字符，数据框
class(content)
#$后跟列名，提取列数据
# class(x$)
# substr提取
substr('xunlianying',1,4)
# c()创建集合
jihe <- c('a','b','c')
# 提取/删除数据库行列
q <- content[,1:3]
p <- content[1:3,]
r <- content[1:3,5:10]
s <- content[-1,]#删除第一行
t <- content[,-1]#删除第一列
u <- content[,-(2:4)]#连续删除
v <- content[,-c(1,3)]#非连续删除
# %>%传导符
x <- content %>% t() %>% as.data.frame()
# duplicated函数
dup_data <- c('a','b','a','b','c')
duplicated(dup_data)
r_dup_data <- dup_data[!duplicated(dup_data)]#取不重复（唯一）数据
r_dup_data
# inner_join合并
# tribble创建简易数据框
class1 <- tribble(~'名次',~'姓名',
                  '1','xiaoming',
                  '2','xiaoming2',
                  '3','xiaoming3',
                  )

class2 <- tribble(~'名次',~'姓名',
                  '1','xiaoli',
                  '2','xiaoli2',
                  '4','xiaoli3',
)

class3 <- tribble(~'名次',~'姓名',~'haha',
                  '1','xiaoli',
                  '2','xiaoli2',
                  '3','xiaoli3',
)
class1
class2
inner_join(class1,class2,by='名次')
left_join(class1,class2,by='名次')
right_join(class1,class2,by='名次')

# TCGA-LUAD数据下载
setwd('TCGA-LUAD')
setwd('TCGAdata')




