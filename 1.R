library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
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
count_CD4<-count
metadata_CD4<-metadata
remove(count,metadata)

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
count_CD8<-count
metadata_CD8<-metadata
remove(count,metadata)

metadata<-rbind(metadata_CD4,metadata_CD8)
count<-cbind(count_CD4,count_CD8)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
f<-readRDS("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/integration/int.CD8.S35.meta.tb.rds")
g<-rbind(e,f)

remove(count_CD4,count_CD8,metadata_CD4,metadata_CD8,e,f)
h<-intersect(g$cellID.uniq,metadata$cellID.uniq)
g<-as.data.frame(g)
rownames(g)<-g$cellID.uniq
rownames(metadata)<-metadata$cellID.uniq
g_new<-g[h,]
metadata_new<-metadata[h,]
remove(metadata,g,h)
metadata<-cbind(g_new,metadata_new)
remove(g_new,metadata_new)
metadata<-metadata[colnames(count),]

Tumor_metadata<-metadata[which(metadata$loc == "T"),]
Normal_metadata<-metadata[which(metadata$loc == "N"),]
PBMC_metadata<-metadata[which(metadata$loc == "P"),]
L_metadata<-metadata[which(metadata$loc == "L"),]

Treg_Tumor_metadata<-Tumor_metadata[which(Tumor_metadata$meta.cluster == "CD4.c18.Treg.RTKN2" |Tumor_metadata$meta.cluster == "CD4.c19.Treg.S1PR1"|Tumor_metadata$meta.cluster == "CD4.c20.Treg.TNFRSF9"|Tumor_metadata$meta.cluster == "CD4.c21.Treg.OAS1"),]
Treg_Normal_metadata<-Normal_metadata[which(Normal_metadata$meta.cluster == "CD4.c18.Treg.RTKN2" |Normal_metadata$meta.cluster == "CD4.c19.Treg.S1PR1"|Normal_metadata$meta.cluster == "CD4.c20.Treg.TNFRSF9"|Normal_metadata$meta.cluster == "CD4.c21.Treg.OAS1"),]
Treg_PBMC_metadata<-PBMC_metadata[which(PBMC_metadata$meta.cluster == "CD4.c18.Treg.RTKN2" |PBMC_metadata$meta.cluster == "CD4.c19.Treg.S1PR1"|PBMC_metadata$meta.cluster == "CD4.c20.Treg.TNFRSF9"|PBMC_metadata$meta.cluster == "CD4.c21.Treg.OAS1"),]

Treg_Tumor_count<-count[,rownames(Treg_Tumor_metadata)]
Treg_Normal_count<-count[,rownames(Treg_Normal_metadata)]
Treg_PBMC_count<-count[,rownames(Treg_PBMC_metadata)]

# JMJD1C_Treg_Tumor_count<-Treg_Tumor_count[which(rownames(Treg_Tumor_count) == "JMJD1C"),]
# JMJD1C_Treg_Normal_count<-Treg_Normal_count[which(rownames(Treg_Normal_count) == "JMJD1C"),]
# JMJD1C_Treg_PBMC_count<-Treg_PBMC_count[which(rownames(Treg_PBMC_count) == "JMJD1C"),]
# 
# IFNG_Treg_Tumor_count<-Treg_Tumor_count[which(rownames(Treg_Tumor_count) == "IFNG"),]
patient<-unique(Treg_Tumor_metadata$patient.uid)
z<-data.frame(0,0)
colnames(z)<-c("JMJD1C","IFNG")#
for (i in 1:length(patient)) {
  a<-Treg_Tumor_metadata[which(Treg_Tumor_metadata$patient.uid == patient[i]),]
  if(length(rownames(a))==1){
    b<-as.data.frame(Treg_Tumor_count[,rownames(a)])
    rownames(b)<-rownames(Treg_Tumor_count)
    JMJD1C<-as.data.frame(t(b[which(rownames(b) == "JMJD1C"),])) #
    IFNG<-as.data.frame(t(b[which(rownames(b) == "IFNG"),]))    #
    colnames(JMJD1C)<-"JMJD1C"
    colnames(IFNG)<-"IFNG"
    c<-data.frame(sum(JMJD1C$JMJD1C),sum(IFNG$IFNG))             #
    colnames(c)<-c("JMJD1C","IFNG")                              #
    z<-rbind(z,c)
    remove(a,b,JMJD1C,IFNG,c)
  } else {
    b<-as.data.frame(Treg_Tumor_count[,rownames(a)])
    JMJD1C<-as.data.frame(t(b[which(rownames(b) == "JMJD1C"),])) #
    IFNG<-as.data.frame(t(b[which(rownames(b) == "IFNG"),]))    #
    colnames(JMJD1C)<-"JMJD1C"
    colnames(IFNG)<-"IFNG"
    c<-data.frame(sum(JMJD1C$JMJD1C),sum(IFNG$IFNG))             #
    colnames(c)<-c("JMJD1C","IFNG")                              #
    z<-rbind(z,c)
    remove(a,b,JMJD1C,IFNG,c)
  }
    
  }
remove(i,patient)
z<-z[-1,]
z_new<-z[which(z[,1] != "NA" & z[,2] != "NA"),]
z_new<-z_new[which(z_new$JMJD1C >0 & z_new$IFNG),]
library(ggplot2)
library(ggpubr)
ggplot(data = z_new,aes(x=JMJD1C,y=IFNG)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = z_new,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
write.table(scaldata,file = "D:/GSE156728/1C_TregIFNG.txt",sep = "\t")