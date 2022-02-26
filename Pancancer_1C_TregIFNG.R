library(dplyr)
library(Seurat)
library(patchwork)

BC_CD4.data <- read.table("D:/GSE156728/GSE156728_BC_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
BCL_CD4.data <- read.table("D:/GSE156728/GSE156728_BCL_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
ESCA_CD4.data <- read.table("D:/GSE156728/GSE156728_ESCA_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
MM_CD4.data <- read.table("D:/GSE156728/GSE156728_MM_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
PACA_CD4.data <- read.table("D:/GSE156728/GSE156728_PACA_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
RC_CD4.data <- read.table("D:/GSE156728/GSE156728_RC_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
THCA_CD4.data <- read.table("D:/GSE156728/GSE156728_THCA_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
UCEC_CD4.data <- read.table("D:/GSE156728/GSE156728_UCEC_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
OV_CD4.data <- read.table("D:/GSE156728/GSM4743199_OV_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
FTC_CD4.data <- read.table("D:/GSE156728/GSM4743231_FTC_10X.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)
CHOL_CD4.data <- read.table("D:/GSE156728/GSM4743237_CHOL_SS2.CD4.counts.txt",sep = "\t",header = TRUE,row.names = 1)

for (i in 1:length(colnames(BC_CD4.data))) {
  colnames(BC_CD4.data)[i] <- paste(colnames(BC_CD4.data)[i],"BC_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(BCL_CD4.data))) {
  colnames(BCL_CD4.data)[i] <- paste(colnames(BCL_CD4.data)[i],"BCL_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(ESCA_CD4.data))) {
  colnames(ESCA_CD4.data)[i] <- paste(colnames(ESCA_CD4.data)[i],"ESCA_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(MM_CD4.data))) {
  colnames(MM_CD4.data)[i] <- paste(colnames(MM_CD4.data)[i],"MM_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(PACA_CD4.data))) {
  colnames(PACA_CD4.data)[i] <- paste(colnames(PACA_CD4.data)[i],"PACA_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(RC_CD4.data))) {
  colnames(RC_CD4.data)[i] <- paste(colnames(RC_CD4.data)[i],"RC_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(THCA_CD4.data))) {
  colnames(THCA_CD4.data)[i] <- paste(colnames(THCA_CD4.data)[i],"THCA_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(UCEC_CD4.data))) {
  colnames(UCEC_CD4.data)[i] <- paste(colnames(UCEC_CD4.data)[i],"UCEC_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(OV_CD4.data))) {
  colnames(OV_CD4.data)[i] <- paste(colnames(OV_CD4.data)[i],"OV_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(FTC_CD4.data))) {
  colnames(FTC_CD4.data)[i] <- paste(colnames(FTC_CD4.data)[i],"FTC_CD4",i,sep = "-")  
}
for (i in 1:length(colnames(CHOL_CD4.data))) {
  colnames(CHOL_CD4.data)[i] <- paste(colnames(CHOL_CD4.data)[i],"CHOL_CD4",i,sep = "-")  
}

A<-cbind(BC_CD4.data,ESCA_CD4.data,FTC_CD4.data,OV_CD4.data,PACA_CD4.data,RC_CD4.data,THCA_CD4.data,UCEC_CD4.data)
B<-cbind(BCL_CD4.data,MM_CD4.data)
#C<-intersect(rownames(A),rownames(B))
#B_new<-B[C,]
#D<-cbind(A,B_new)
#E<-intersect(rownames(D),rownames(CHOL_CD4.data))
#D_new<-D[E,]
#CHOL_CD4.data_new<-CHOL_CD4.data[E,]
#F<-cbind(D_new,CHOL_CD4.data_new)
#remove(BC_CD4.data,BCL_CD4.data,CHOL_CD4.data,CHOL_CD4.data_new,ESCA_CD4.data,FTC_CD4.data,MM_CD4.data,OV_CD4.data,PACA_CD4.data,RC_CD4.data,THCA_CD4.data,UCEC_CD4.data,C,E,i,B_new,D_new)
metadata<-read.table("D:/GSE156728/GSE156728_metadata.txt",sep = "\t",header = TRUE,row.names = 1)
G<-as.data.frame(colnames(A))
colnames(G)<-"A"
library(tidyr)
H<-separate(G,A,into = c("A","B","C"),sep = "-")
I<-intersect(H$A,rownames(metadata))
metadata_new<-metadata[I,]
J<-metadata_new
rownames(J)<-colnames(A)

pbmc <- CreateSeuratObject(counts = A, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = J )
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:31)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- RunUMAP(pbmc, dims = 1:31)
DimPlot(pbmc, reduction = "umap",label = TRUE)
VlnPlot(pbmc,features = "FOXP3",pt.size = 0,sort = TRUE)
FeaturePlot(pbmc,features = "FOXP3",label = TRUE)
Treg<-subset(pbmc,idents = c("21","6","4","1","18","10","14","20","17"))
Idents(Treg)<-Treg@meta.data$patient
DimPlot(Treg, reduction = "umap")

#a<-read.table("D:/GSE156728/geneset (1).txt",sep = "\t",header = TRUE)
#b<-read.table("D:/GSE156728/geneset.txt",sep = "\t",header = TRUE)
#a_new<-a[2:60,]
#b_new<-b[2:92,]
#Treg<-AddModuleScore(Treg,features = a_new,name = "STAT3")
#Treg<-AddModuleScore(Treg,features = b_new,name = "PI3K")


Idents(Treg)<-Treg@meta.data$loc
Treg_T<-subset(Treg,idents = "T")
Idents(Treg_T)<-Treg_T@meta.data$patient
patient<-unique(Treg_T@meta.data$patient)
for (i in 1:length(patient)) {
  assign(paste("Treg_",i,"_TIL",sep = ""),subset(Treg_T,idents = patient[i]))
  assign("a",slot(get(paste("Treg_",i,"_TIL",sep = "")),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@data))))
  assign("c",select(b,one_of("JMJD1C")))
  assign("d",select(b,one_of("IFNG")))
  assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(mean(c$JMJD1C),mean(d$IFNG)))
  
}
scaldata<-scaldata_1_TIL
for (i in 2:length(patient)) {
  scaldata<-rbind(scaldata,get(paste("scaldata_",i,"_TIL",sep = "")))
}
colnames(scaldata)<-c("JMJD1C","IFNG")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=IFNG)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
write.table(scaldata,file = "D:/GSE156728/1C_TregIFNG.txt",sep = "\t")
