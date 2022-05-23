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
remove(BC_CD4.data,BCL_CD4.data,CHOL_CD4.data,ESCA_CD4.data,FTC_CD4.data,MM_CD4.data,OV_CD4.data,PACA_CD4.data,RC_CD4.data,THCA_CD4.data,UCEC_CD4.data,i)
G<-as.data.frame(colnames(A))
colnames(G)<-"A"
library(tidyr)
H<-separate(G,A,into = c("A","B","C"),sep = "-")
I<-intersect(H$A,rownames(metadata))
metadata_new<-metadata[I,]
J<-metadata_new
rownames(J)<-colnames(A)

pbmc <- CreateSeuratObject(counts = A, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = J )
remove(B,G,H,metadata,metadata_new,I)
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
#DimPlot(pbmc, reduction = "umap",label = TRUE)

Idents(pbmc)<-pbmc@meta.data$meta.cluster
b<-read.table("D:/GSE156728/Reactome_NEGTIVE_PI3K_AKT.txt",sep = "\t",header = TRUE)
b_new<-b[2:length(rownames(b)),]
b_new<-intersect(b_new,rownames(pbmc))
pbmc<-AddModuleScore(pbmc,features = b_new,name = "Reactome_NEGTIVE_PI3K_AKT")

Treg<-subset(pbmc,idents = c("CD4.c21.Treg.OAS1","CD4.c20.Treg.TNFRSF9","CD4.c19.Treg.S1PR1","CD4.c18.Treg.RTKN2","CD4.c23.Mix.NME1","CD4.c24.Mix.NME2"))
Idents(Treg)<-Treg@meta.data$loc
Treg_T<-subset(Treg,idents = "T")



Idents(Treg_T)<-Treg_T@meta.data$patient
patient<-unique(Treg_T@meta.data$patient)
for (i in 1:length(patient)) {
  assign(paste("Treg_",i,"_TIL",sep = ""),subset(Treg_T,idents = patient[i]))
  assign("a",slot(get(paste("Treg_",i,"_TIL",sep = "")),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@data))))
  assign("c",select(b,one_of("JMJD1C")))
  assign("d",slot(get(paste("Treg_",i,"_TIL",sep = "")),"meta.data"))
  PI3K<-select(d,starts_with("Reactome_NEGTIVE_PI3K_AKT"))
  d$PI3K<-apply(PI3K,MARGIN = 1,median)
  assign("e",as.data.frame(d$PI3K))
  colnames(e)<-"PI3K"
  assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(mean(c$JMJD1C),sum(e$PI3K)))
  
}
scaldata<-scaldata_1_TIL
for (i in 2:length(patient)) {
  scaldata<-rbind(scaldata,get(paste("scaldata_",i,"_TIL",sep = "")))
}
colnames(scaldata)<-c("JMJD1C","PI3K_signaling")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=PI3K_signaling)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

write.table(scaldata,file = "1C_PI3Kneg_median_median.txt",sep = "\t")

remove(pbmc_CD8T_1_TIL,pbmc_CD8T_10_TIL,pbmc_CD8T_11_TIL,pbmc_CD8T_12_TIL,pbmc_CD8T_13_TIL,pbmc_CD8T_14_TIL,pbmc_CD8T_15_TIL,pbmc_CD8T_16_TIL,pbmc_CD8T_17_TIL,pbmc_CD8T_18_TIL,pbmc_CD8T_19_TIL,pbmc_CD8T_2_TIL,
       pbmc_CD8T_20_TIL,pbmc_CD8T_21_TIL,pbmc_CD8T_22_TIL,pbmc_CD8T_23_TIL,pbmc_CD8T_24_TIL,pbmc_CD8T_25_TIL,pbmc_CD8T_26_TIL,pbmc_CD8T_27_TIL,pbmc_CD8T_28_TIL,pbmc_CD8T_29_TIL,pbmc_CD8T_3_TIL,pbmc_CD8T_30_TIL,pbmc_CD8T_31_TIL,pbmc_CD8T_32_TIL,pbmc_CD8T_33_TIL,pbmc_CD8T_34_TIL,pbmc_CD8T_35_TIL,pbmc_CD8T_36_TIL,pbmc_CD8T_37_TIL,pbmc_CD8T_38_TIL,pbmc_CD8T_39_TIL,
       pbmc_CD8T_4_TIL,pbmc_CD8T_40_TIL,pbmc_CD8T_41_TIL,pbmc_CD8T_42_TIL,pbmc_CD8T_5_TIL,pbmc_CD8T_6_TIL,pbmc_CD8T_7_TIL,pbmc_CD8T_8_TIL,pbmc_CD8T_9_TIL)
remove(Treg_1_TIL,Treg_10_TIL,Treg_11_TIL,Treg_12_TIL,Treg_13_TIL,Treg_14_TIL,Treg_15_TIL,Treg_16_TIL,Treg_17_TIL,Treg_18_TIL,Treg_19_TIL,Treg_2_TIL,
       Treg_20_TIL,Treg_21_TIL,Treg_22_TIL,Treg_23_TIL,Treg_24_TIL,Treg_25_TIL,Treg_26_TIL,Treg_27_TIL,Treg_28_TIL,Treg_29_TIL,Treg_3_TIL,Treg_30_TIL,Treg_31_TIL,Treg_32_TIL,Treg_33_TIL,Treg_34_TIL,Treg_35_TIL,Treg_36_TIL,Treg_37_TIL,Treg_38_TIL,Treg_39_TIL,
       Treg_4_TIL,Treg_40_TIL,Treg_41_TIL,Treg_42_TIL,Treg_5_TIL,Treg_6_TIL,Treg_7_TIL,Treg_8_TIL,Treg_9_TIL,Treg_43_TIL,Treg_44_TIL,Treg_45_TIL,Treg_46_TIL,Treg_47_TIL,Treg_48_TIL)
remove(scaldata_1_TIL,scaldata_10_TIL,scaldata_11_TIL,scaldata_12_TIL,scaldata_13_TIL,scaldata_14_TIL,scaldata_15_TIL,scaldata_16_TIL,scaldata_17_TIL,scaldata_18_TIL,scaldata_19_TIL,scaldata_2_TIL,
       scaldata_20_TIL,scaldata_21_TIL,scaldata_22_TIL,scaldata_23_TIL,scaldata_24_TIL,scaldata_25_TIL,scaldata_26_TIL,scaldata_27_TIL,scaldata_28_TIL,scaldata_29_TIL,scaldata_3_TIL,scaldata_30_TIL,scaldata_31_TIL,scaldata_32_TIL,scaldata_33_TIL,scaldata_34_TIL,scaldata_35_TIL,scaldata_36_TIL,scaldata_37_TIL,scaldata_38_TIL,scaldata_39_TIL,
       scaldata_4_TIL,scaldata_40_TIL,scaldata_41_TIL,scaldata_42_TIL,scaldata_5_TIL,scaldata_6_TIL,scaldata_7_TIL,scaldata_8_TIL,scaldata_9_TIL,scaldata_43_TIL,scaldata_44_TIL,scaldata_45_TIL,scaldata_46_TIL,scaldata_47_TIL,scaldata_48_TIL)
remove(a,b,c,d,e,f,i,patient)
remove(PI3K,STAT3,scaldata,Treg,Treg_T,b_new)
