library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
count<-data.frame()
metadata<-data.frame()
gene<-character()
c<-readRDS(a[1])
count<-c@assays$data@listData$norm_exprs
count<-as.matrix(count)
count<-as.data.frame(count)
metadata<-as.data.frame(c@colData@listData)
count<-count[,metadata$cellID]
colnames(count)<-metadata$cellID.uniq
gene<-rownames(count)
for (i in c(2:23,25:length(a))) {
  b<-readRDS(a[i])
  count_B<-b@assays$data@listData$norm_exprs
  count_B<-as.matrix(count_B)
  count_B<-as.data.frame(count_B)
  count_B_metadata<-as.data.frame(b@colData@listData)
  count_B<-count_B[,count_B_metadata$cellID]
  colnames(count_B)<-count_B_metadata$cellID.uniq
  genename<-rownames(count_B)
  d<-str_sub(genename[1],1,4)
  if(d == "ENSG"){
    e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
    genename<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    genename<-g$SYMBOL
  }
  gene<-intersect(gene,genename)
  count_B<-count_B[gene,]
  count<-count[gene,]
  count<-cbind(count,count_B)
  metadata<-rbind(metadata,count_B_metadata)
  remove(count_B,count_B_metadata,b,d)
}
b<-readRDS(a[24])
count_B<-b@assays@data$norm_exprs
count_B<-as.matrix(count_B)
count_B<-as.data.frame(count_B)
count_B_metadata<-as.data.frame(b@colData@listData)
count_B<-count_B[,count_B_metadata$cellID]
colnames(count_B)<-count_B_metadata$cellID.uniq
genename<-rownames(count_B)
d<-str_sub(genename[1],1,4)
if(d == "ENSG"){
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
  genename<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  genename<-g$SYMBOL
}
gene<-intersect(gene,genename)
count_B<-count_B[gene,]
count<-count[gene,]
count<-cbind(count,count_B)
metadata<-rbind(metadata,count_B_metadata)
remove(count_B,count_B_metadata,b,d)
remove(c,e,g,a,f,gene,genename,i)
