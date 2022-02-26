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

BC_CD8.data <- read.table("D:/GSE156728/GSE156728_BC_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
BCL_CD8.data <- read.table("D:/GSE156728/GSE156728_BCL_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
ESCA_CD8.data <- read.table("D:/GSE156728/GSE156728_ESCA_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
MM_CD8.data <- read.table("D:/GSE156728/GSE156728_MM_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
PACA_CD8.data <- read.table("D:/GSE156728/GSE156728_PACA_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
RC_CD8.data <- read.table("D:/GSE156728/GSE156728_RC_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
THCA_CD8.data <- read.table("D:/GSE156728/GSE156728_THCA_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
UCEC_CD8.data <- read.table("D:/GSE156728/GSE156728_UCEC_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
OV_CD8.data <- read.table("D:/GSE156728/GSM4743199_OV_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
FTC_CD8.data <- read.table("D:/GSE156728/GSM4743231_FTC_10X.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)
CHOL_CD8.data <- read.table("D:/GSE156728/GSM4743237_CHOL_SS2.CD8.counts.txt",sep = "\t",header = TRUE,row.names = 1)

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

for (i in 1:length(colnames(BC_CD8.data))) {
  colnames(BC_CD8.data)[i] <- paste(colnames(BC_CD8.data)[i],"BC_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(BCL_CD8.data))) {
  colnames(BCL_CD8.data)[i] <- paste(colnames(BCL_CD8.data)[i],"BCL_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(ESCA_CD8.data))) {
  colnames(ESCA_CD8.data)[i] <- paste(colnames(ESCA_CD8.data)[i],"ESCA_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(MM_CD8.data))) {
  colnames(MM_CD8.data)[i] <- paste(colnames(MM_CD8.data)[i],"MM_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(PACA_CD8.data))) {
  colnames(PACA_CD8.data)[i] <- paste(colnames(PACA_CD8.data)[i],"PACA_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(RC_CD8.data))) {
  colnames(RC_CD8.data)[i] <- paste(colnames(RC_CD8.data)[i],"RC_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(THCA_CD8.data))) {
  colnames(THCA_CD8.data)[i] <- paste(colnames(THCA_CD8.data)[i],"THCA_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(UCEC_CD8.data))) {
  colnames(UCEC_CD8.data)[i] <- paste(colnames(UCEC_CD8.data)[i],"UCEC_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(OV_CD8.data))) {
  colnames(OV_CD8.data)[i] <- paste(colnames(OV_CD8.data)[i],"OV_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(FTC_CD8.data))) {
  colnames(FTC_CD8.data)[i] <- paste(colnames(FTC_CD8.data)[i],"FTC_CD8",i,sep = "-")  
}
for (i in 1:length(colnames(CHOL_CD8.data))) {
  colnames(CHOL_CD8.data)[i] <- paste(colnames(CHOL_CD8.data)[i],"CHOL_CD8",i,sep = "-")  
}

A<-cbind(BC_CD4.data,ESCA_CD4.data,FTC_CD4.data,OV_CD4.data,PACA_CD4.data,RC_CD4.data,THCA_CD4.data,UCEC_CD4.data,
         BC_CD8.data,ESCA_CD8.data,FTC_CD8.data,OV_CD8.data,PACA_CD8.data,RC_CD8.data,THCA_CD8.data,UCEC_CD8.data)
B<-cbind(BCL_CD4.data,MM_CD4.data,
         BCL_CD8.data,MM_CD8.data)
remove(BC_CD4.data,ESCA_CD4.data,FTC_CD4.data,OV_CD4.data,PACA_CD4.data,RC_CD4.data,THCA_CD4.data,UCEC_CD4.data,
       BC_CD8.data,ESCA_CD8.data,FTC_CD8.data,OV_CD8.data,PACA_CD8.data,RC_CD8.data,THCA_CD8.data,UCEC_CD8.data,
       BCL_CD4.data,MM_CD4.data,
       BCL_CD8.data,MM_CD8.data)
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
K<-separate(H,B,into = c("B","D"),sep = "_")
I<-intersect(K$A,rownames(metadata))
metadata_new<-metadata[I,]
rownames(K)<-K$A
K_new<-K[rownames(metadata_new),]
J<-cbind(metadata_new,K_new)
rownames(J)<-colnames(A)
remove(metadata,metadata_new,G,H,K,K_new)
remove(I)

pbmc <- CreateSeuratObject(counts = A, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = J )
remove(A,J)
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

Idents(pbmc)<-pbmc@meta.data$D
pbmc_CD8<-subset(pbmc,idents = "CD8")
pbmc_CD4<-subset(pbmc,idents = "CD4")
pbmc_CD4 <- RunPCA(pbmc_CD4, features = VariableFeatures(object = pbmc_CD4))
ElbowPlot(pbmc_CD4)
pbmc_CD4 <- FindNeighbors(pbmc_CD4, dims = 1:31)
pbmc_CD4 <- FindClusters(pbmc_CD4, resolution = 1)
pbmc_CD4 <- RunUMAP(pbmc_CD4, dims = 1:31)
DimPlot(pbmc_CD4, reduction = "umap",label = TRUE)
VlnPlot(pbmc_CD4,features = "FOXP3",sort = TRUE,pt.size = 0)
FeaturePlot(pbmc_CD4,features = "FOXP3",label = TRUE)
Idents(pbmc_CD4)<-pbmc_CD4@meta.data$meta.cluster
DimPlot(pbmc_CD4)
VlnPlot(pbmc_CD4,features = "FOXP3",sort = TRUE,pt.size = 0)
FeaturePlot(pbmc_CD4,features = "FOXP3",label = TRUE)
Idents(pbmc)<-pbmc@meta.data$meta.cluster
VlnPlot(pbmc,features = "FOXP3",sort = TRUE,pt.size = 0)
FeaturePlot(pbmc,features = "FOXP3",label = TRUE)
remove(pbmc)
Treg<-subset(pbmc_CD4,idents = c("CD4.c21.Treg.OAS1","CD4.c20.Treg.TNFRSF9","CD4.c19.Treg.S1PR1","CD4.c18.Treg.RTKN2","CD4.c23.Mix.NME1","CD4.c24.Mix.NME2"))
Idents(pbmc_CD8)<-pbmc_CD8@meta.data$loc
pbmc_CD8_T<-subset(pbmc_CD8,idents = "T")
remove(pbmc_CD8,pbmc_CD4)
Idents(Treg)<-Treg@meta.data$loc
Treg_T<-subset(Treg,idents = "T")
remove(Treg)
Idents(pbmc_CD8_T)<-pbmc_CD8_T@meta.data$patient
Idents(Treg_T)<-Treg_T@meta.data$patient
patient<-unique(Treg_T@meta.data$patient)

for (i in 1:length(patient)) {
  assign(paste("Treg_",i,"_TIL",sep = ""),subset(Treg_T,idents = patient[i]))
  assign("a",slot(get(paste("Treg_",i,"_TIL",sep = "")),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@data))))
  assign("c",select(b,one_of("JMJD1C")))
  assign(paste("pbmc_CD8T_",i,"_TIL",sep = ""),subset(pbmc_CD8_T,idents = patient[i]))
  assign("d",length(colnames(get(paste("pbmc_CD8T_",i,"_TIL",sep = ""))))/length(colnames(get(paste("Treg_",i,"_TIL",sep = "")))))
  assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(sum(c$JMJD1C),d))
  
}
scaldata<-scaldata_1_TIL
for (i in 2:length(patient)) {
  scaldata<-rbind(scaldata,get(paste("scaldata_",i,"_TIL",sep = "")))
}
colnames(scaldata)<-c("JMJD1C","CD8_Treg")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=CD8_Treg)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
write.table(scaldata,file = "D:/GSE156728/1C_TregIFNG.txt",sep = "\t")
