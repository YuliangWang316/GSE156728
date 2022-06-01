library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/zhangLab/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/zhangLab/")
# count<-data.frame()
metadata<-data.frame()
# gene<-character()
c<-readRDS(a[1])
# count<-c@assays$data@listData$norm_exprs
# count<-as.matrix(count)
# count<-as.data.frame(count)
metadata<-as.data.frame(c@colData@listData)
d<-colnames(metadata)
# gene<-rownames(count)
# d<-str_sub(gene[1],1,4)
# if(d == "ENSG"){
#   e<-bitr(geneID = gene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
#   gene<-e$SYMBOL
# }
# remove(d,e)
# for (i in c(2:length(a))) {
#   b<-readRDS(a[i])
#   count_B<-b@assays$data@listData$norm_exprs
#   count_B<-as.matrix(count_B)
#   count_B<-as.data.frame(count_B)
#   count_B_meatadata<-as.data.frame(b@colData@listData)
#   genename<-rownames(count_B)
#   d<-str_sub(genename[1],1,4)
#   if(d == "ENSG"){
#     e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
#     genename<-e$SYMBOL
#   }
#   f<-sort(genename)
#   if(f[1] == "1"){
#     g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
#     genename<-g$SYMBOL
#   }
#   gene<-intersect(gene,genename)
#   count_B<-count_B[gene,]
#   count<-count[gene,]
#   count<-cbind(count,count_B)
#   metadata<-rbind(metadata,count_B_meatadata)
#   remove(count_B,count_B_meatadata,b,d)
# }
# 
# 
# remove(count_B,count_B_meatadata,b,d)
# remove(c,e,g,a,f,gene,genename,i)

for (i in c(2:length(a))) {
  b<-readRDS(a[i])
  count_B_meatadata<-as.data.frame(b@colData@listData)
  # count_B_meatadata<-count_B_meatadata[,d]
  # metadata<-rbind(metadata,count_B_meatadata)
  
}
write.table(metadata,"zhang.txt",sep = "\t")
