a<-read.table("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4_metadata.txt",sep = "\t")
b<-read.table("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8_metadata.txt",sep = "\t")
c<-read.table("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4_count.txt",sep = "\t",header = TRUE,row.names = 1)
d<-read.table("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8_count.txt",sep = "\t",header = TRUE,row.names = 1)
metadata<-rbind(a,b)
count<-cbind(c,d)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
f<-readRDS("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/integration/int.CD8.S35.meta.tb.rds")
g<-rbind(e,f)
remove(a,b,c,d,e,f)
h<-intersect(g$cellID.uniq,metadata$cellID.uniq)
rownames(g)<-g$cellID.uniq
rownames(metadata)<-metadata$cellID.uniq
g_new<-g[h,]
metadata_new<-metadata[h,]
remove(metadata,g,h)
metadata<-cbind(g_new,metadata_new)
remove(g_new,metadata_new)
count_new<-count[,rownames(metadata)]
