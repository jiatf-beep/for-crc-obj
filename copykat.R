rm(list=ls())
gc()
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
Epi <- readRDS('F:/rds-liu/8-26/Epi.rds')
Bcell<- readRDS('F:/rds-liu/8-26/Bcell.rds')
Tcell <- readRDS('F:/rds-liu/8-26/allt.rds')

epiMat=as.data.frame(GetAssayData( Epi , slot='counts',assay='RNA'))
 # 这里面仅仅是B淋巴 细胞哦
B_cellMat=as.data.frame(GetAassayData(Bcell, slot='counts',assay='RNA'))
# 这里面仅仅是T淋巴 细胞哦
T_cellsMat=as.data.frame(GetAssayData( Tcell, slot='counts',assay='RNA'))

B_cellMat=B_cellMat[,sample(1:ncol(B_cellMat),800)]
T_cellsMat=T_cellsMat[,sample(1:ncol(T_cellsMat),800)]
dat=cbind(epiMat,B_cellMat,T_cellsMat)


groupinfo=data.frame(v1=colnames(dat),
                     v2=c(rep('epi',ncol(epiMat)),
                          rep('spike-B_cell',300),
                          rep('ref-B_cell',500),
                          rep('spike-T_cells',300),
                          rep('ref-T_cells',500)))
loc2 <- which(groupinfo$v1 %in% data3$cell.names)
groupinfo$v2[loc2] <- 'tumor'
unique(groupinfo$v2)

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
dir.create('F:/rds-liu/2024.6.20/infercnv')
setwd('F:/rds-liu/2024.6.20/infercnv')
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,
            row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,
            row.names = F)

#copykat
library(Seurat)
library(gplots)
library(ggplot2)

expFile='expFile.txt' 
groupFiles='groupFiles.txt'  
geneFile='geneFile.txt' 

library(copykat)
# library(devtools)
# install_github("navinlabcode/copykat")
exp.File <- read.delim( expFile , row.names=1)
exp.File[1:4,1:4]
groupFiles <- read.delim( groupFiles , header=FALSE)
head(groupFiles)

head( colnames(exp.File))
colnames(exp.File)=gsub('^X','', colnames(exp.File))
colnames(exp.File)=gsub('[.]','-', colnames(exp.File))
head( colnames(exp.File))

table(groupFiles$V2)
normalCells <-  colnames(exp.File)[grepl("ref",groupFiles$V2)]
head(normalCells)
Sys.time()
res <- copykat(rawmat=exp.File,  # 准备好的表达量矩阵 
               id.type="S", # 选择 symbol，因为我们的表达量矩阵 里面是它
               ngene.chr=5, 
               win.size=25, 
               KS.cut=0.1, 
               sam.name="test",  # 输出的一系列文件的前缀
               distance="euclidean", 
               norm.cell.names=normalCells, # 正常细胞的名字
               n.cores=1,
               output.seg="FLASE")
Sys.time()
save(res,file = 'ref-of-copykat.Rdata')





