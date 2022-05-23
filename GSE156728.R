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

C<-intersect(rownames(A),rownames(B))
B_new<-B[C,]
D<-cbind(A,B_new)
#E<-intersect(rownames(D),rownames(CHOL_CD4.data))
#D_new<-D[E,]
#CHOL_CD4.data_new<-CHOL_CD4.data[E,]
#G<-cbind(D_new,CHOL_CD4.data_new)
#H<-intersect(rownames(G),rownames(CHOL_CD8.data))
#G_new<-G[H,]
#CHOL_CD8.data_new<-CHOL_CD8.data[H,]
#I<-cbind(G_new,CHOL_CD8.data_new)
remove(A,B,B_new,CHOL_CD4.data,CHOL_CD4.data_new,CHOL_CD8.data,CHOL_CD8.data_new,D_new,G,G_new,C,E,H,i)
metadata<-read.table("D:/GSE156728/GSE156728_metadata.txt",sep = "\t",header = TRUE,row.names = 1)
G<-as.data.frame(colnames(D))
colnames(G)<-"A"
library(tidyr)
H<-separate(G,A,into = c("A","B","C"),sep = "-")
K<-separate(H,B,into = c("B","D"),sep = "_")
J<-intersect(K$A,rownames(metadata))
metadata_new<-metadata[J,]
rownames(K)<-K$A
K_new<-K[rownames(metadata_new),]
L<-cbind(metadata_new,K_new)
rownames(L)<-colnames(D)
remove(metadata,metadata_new,G,H,K,K_new)
remove(J)


for (i in 1:26) {
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HNSCC_",i,"_PBMC",sep = "")))
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  g<-get(paste("HNSCC_",i,"_PBMC.data",sep = ""))
  for (j in 1:length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))) {
    colnames(g)[j]<-paste(colnames(g[j]),"PBMC",i,j,sep = "-")
  }
  assign(paste("HNSCC_",i,"_PBMC.data",sep = ""),g)
  assign(paste("HNSCC_",i,"_PBMC.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))),rep(paste("PBMC",i,sep = ""),length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))),rep("PBMC",length(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = "")))))))
  a<-get(paste("HNSCC_",i,"_PBMC.metadata",sep = ""))
  b<-as.data.frame(colnames(get(paste("HNSCC_",i,"_PBMC.data",sep = ""))))
  colnames(a)<-c("barcode","group","subgroup")
  rownames(a)<-b[,1]
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),CreateSeuratObject(counts = g, project = paste("HNSCC_",i,"_PBMC",sep = ""),meta.data = a))
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),NormalizeData(get(paste("HNSCC_",i,"_PBMC",sep = ""))))
  assign(paste("HNSCC_",i,"_PBMC",sep = ""),FindVariableFeatures(get(paste("HNSCC_",i,"_PBMC",sep = "")), selection.method = "vst", nfeatures = 2000))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),Read10X(paste("D:/Jmjd1c_Treg_Tumor/GSE139324/HNSCC_",i,"_TIL",sep = "")))
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),as.data.frame(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  h<-get(paste("HNSCC_",i,"_TIL.data",sep = ""))
  for (k in 1:length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))) {
    colnames(h)[k]<-paste(colnames(h[k]),"TIL",i,k,sep = "-")
  }
  assign(paste("HNSCC_",i,"_TIL.data",sep = ""),h)
  assign(paste("HNSCC_",i,"_TIL.metadata",sep = ""),data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))),rep(paste("TIL",i,sep = ""),length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))),rep("TIL",length(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = "")))))))
  c<-get(paste("HNSCC_",i,"_TIL.metadata",sep = ""))
  d<-as.data.frame(colnames(get(paste("HNSCC_",i,"_TIL.data",sep = ""))))
  colnames(c)<-c("barcode","group","subgroup")
  rownames(c)<-d[,1]
  assign(paste("HNSCC_",i,"_TIL",sep = ""),CreateSeuratObject(counts = h, project = paste("HNSCC_",i,"_TIL",sep = ""),meta.data = c))
  assign(paste("HNSCC_",i,"_TIL",sep = ""),NormalizeData(get(paste("HNSCC_",i,"_TIL",sep = ""))))
  assign(paste("HNSCC_",i,"_TIL",sep = ""),FindVariableFeatures(get(paste("HNSCC_",i,"_TIL",sep = "")), selection.method = "vst", nfeatures = 2000))
}

e<-paste("HNSCC_",1,"_PBMC.metadata",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_PBMC.metadata",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_PBMC.metadata<-rbind(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)
colnames(HNSCC_PBMC.metadata)<-c("barcodes","group","subgroup")
rownames(HNSCC_PBMC.metadata)<-HNSCC_PBMC.metadata[,1]

e<-paste("HNSCC_",1,"_PBMC.data",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_PBMC.data",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_PBMC.data<-cbind(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)

e<-paste("HNSCC_",1,"_TIL.metadata",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_TIL.metadata",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_TIL.metadata<-rbind(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)
colnames(HNSCC_TIL.metadata)<-c("barcodes","group","subgroup")
rownames(HNSCC_TIL.metadata)<-HNSCC_TIL.metadata[,1]

e<-paste("HNSCC_",1,"_TIL.data",sep="")

for (i in 2:26) {
  f<-paste("HNSCC_",i,"_TIL.data",sep="") 
  e<-paste(e,f,sep = ", ")
}

HNSCC_TIL.data<-cbind(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)


HNSCC.metadata<-rbind(HNSCC_PBMC.metadata,HNSCC_TIL.metadata)
HNSCC.data<-cbind(HNSCC_PBMC.data,HNSCC_TIL.data)

remove(HNSCC_1_PBMC.data, HNSCC_2_PBMC.data, HNSCC_3_PBMC.data, HNSCC_4_PBMC.data, HNSCC_5_PBMC.data, HNSCC_6_PBMC.data, HNSCC_7_PBMC.data, HNSCC_8_PBMC.data, HNSCC_9_PBMC.data, HNSCC_10_PBMC.data, HNSCC_11_PBMC.data, HNSCC_12_PBMC.data, HNSCC_13_PBMC.data, HNSCC_14_PBMC.data, HNSCC_15_PBMC.data, HNSCC_16_PBMC.data, HNSCC_17_PBMC.data, HNSCC_18_PBMC.data, HNSCC_19_PBMC.data, HNSCC_20_PBMC.data, HNSCC_21_PBMC.data, HNSCC_22_PBMC.data, HNSCC_23_PBMC.data, HNSCC_24_PBMC.data, HNSCC_25_PBMC.data, HNSCC_26_PBMC.data)
remove(HNSCC_1_PBMC.metadata, HNSCC_2_PBMC.metadata, HNSCC_3_PBMC.metadata, HNSCC_4_PBMC.metadata, HNSCC_5_PBMC.metadata, HNSCC_6_PBMC.metadata, HNSCC_7_PBMC.metadata, HNSCC_8_PBMC.metadata, HNSCC_9_PBMC.metadata, HNSCC_10_PBMC.metadata, HNSCC_11_PBMC.metadata, HNSCC_12_PBMC.metadata, HNSCC_13_PBMC.metadata, HNSCC_14_PBMC.metadata, HNSCC_15_PBMC.metadata, HNSCC_16_PBMC.metadata, HNSCC_17_PBMC.metadata, HNSCC_18_PBMC.metadata, HNSCC_19_PBMC.metadata, HNSCC_20_PBMC.metadata, HNSCC_21_PBMC.metadata, HNSCC_22_PBMC.metadata, HNSCC_23_PBMC.metadata, HNSCC_24_PBMC.metadata, HNSCC_25_PBMC.metadata, HNSCC_26_PBMC.metadata)
remove(HNSCC_1_PBMC, HNSCC_2_PBMC, HNSCC_3_PBMC, HNSCC_4_PBMC, HNSCC_5_PBMC, HNSCC_6_PBMC, HNSCC_7_PBMC, HNSCC_8_PBMC, HNSCC_9_PBMC, HNSCC_10_PBMC, HNSCC_11_PBMC, HNSCC_12_PBMC, HNSCC_13_PBMC, HNSCC_14_PBMC, HNSCC_15_PBMC, HNSCC_16_PBMC, HNSCC_17_PBMC, HNSCC_18_PBMC, HNSCC_19_PBMC, HNSCC_20_PBMC, HNSCC_21_PBMC, HNSCC_22_PBMC, HNSCC_23_PBMC, HNSCC_24_PBMC, HNSCC_25_PBMC, HNSCC_26_PBMC)
remove(HNSCC_1_TIL, HNSCC_2_TIL, HNSCC_3_TIL, HNSCC_4_TIL, HNSCC_5_TIL, HNSCC_6_TIL, HNSCC_7_TIL, HNSCC_8_TIL, HNSCC_9_TIL, HNSCC_10_TIL, HNSCC_11_TIL, HNSCC_12_TIL, HNSCC_13_TIL, HNSCC_14_TIL, HNSCC_15_TIL, HNSCC_16_TIL, HNSCC_17_TIL, HNSCC_18_TIL, HNSCC_19_TIL, HNSCC_20_TIL, HNSCC_21_TIL, HNSCC_22_TIL, HNSCC_23_TIL, HNSCC_24_TIL, HNSCC_25_TIL, HNSCC_26_TIL)
remove(HNSCC_1_TIL.data, HNSCC_2_TIL.data, HNSCC_3_TIL.data, HNSCC_4_TIL.data, HNSCC_5_TIL.data, HNSCC_6_TIL.data, HNSCC_7_TIL.data, HNSCC_8_TIL.data, HNSCC_9_TIL.data, HNSCC_10_TIL.data, HNSCC_11_TIL.data, HNSCC_12_TIL.data, HNSCC_13_TIL.data, HNSCC_14_TIL.data, HNSCC_15_TIL.data, HNSCC_16_TIL.data, HNSCC_17_TIL.data, HNSCC_18_TIL.data, HNSCC_19_TIL.data, HNSCC_20_TIL.data, HNSCC_21_TIL.data, HNSCC_22_TIL.data, HNSCC_23_TIL.data, HNSCC_24_TIL.data, HNSCC_25_TIL.data, HNSCC_26_TIL.data)
remove(HNSCC_1_TIL.metadata, HNSCC_2_TIL.metadata, HNSCC_3_TIL.metadata, HNSCC_4_TIL.metadata, HNSCC_5_TIL.metadata, HNSCC_6_TIL.metadata, HNSCC_7_TIL.metadata, HNSCC_8_TIL.metadata, HNSCC_9_TIL.metadata, HNSCC_10_TIL.metadata, HNSCC_11_TIL.metadata, HNSCC_12_TIL.metadata, HNSCC_13_TIL.metadata, HNSCC_14_TIL.metadata, HNSCC_15_TIL.metadata, HNSCC_16_TIL.metadata, HNSCC_17_TIL.metadata, HNSCC_18_TIL.metadata, HNSCC_19_TIL.metadata, HNSCC_20_TIL.metadata, HNSCC_21_TIL.metadata, HNSCC_22_TIL.metadata, HNSCC_23_TIL.metadata, HNSCC_24_TIL.metadata, HNSCC_25_TIL.metadata, HNSCC_26_TIL.metadata)
remove(a,b,c,d,e,f,g,h,i,j,k)
remove(HNSCC_PBMC.data,HNSCC_PBMC.metadata,HNSCC_TIL.data,HNSCC_TIL.metadata)

L<-L[,-10]
L<-L[,-9]
L<-L[,-8]
L<-L[,-7]
L<-L[,-6]
L<-L[,-3]
L<-L[,-1]

M<-intersect(rownames(HNSCC.data),rownames(D))
HNSCC.data_new<-HNSCC.data[M,]
D_new<-D[M,]
N<-cbind(HNSCC.data_new,D_new)
remove(HNSCC.data,HNSCC.data_new,D,M,D_new)
HNSCC.metadata<-HNSCC.metadata[,-1]
colnames(L)<-c("group","subgroup","meta.cluster")
HNSCC.metadata$meta.cluster<-""
O<-rbind(L,HNSCC.metadata)
remove(L,HNSCC.metadata)

O<-O[colnames(N),]

pbmc <- CreateSeuratObject(counts = N, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = O )
remove(N,O)
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

#Idents(pbmc)<-pbmc@meta.data$meta.cluster
VlnPlot(pbmc,features = "FOXP3",sort = TRUE,pt.size = 0)
FeaturePlot(pbmc,features = "FOXP3",label = TRUE)
Idents(pbmc)<-pbmc@meta.data$group
DimPlot(pbmc, reduction = "umap",label = TRUE)
Idents(pbmc)<-pbmc@meta.data$seurat_clusters
DotPlot(pbmc, features = "FOXP3") + RotatedAxis()
#FOXP3<-subset(pbmc,idents = c("CD4.c21.Treg.OAS1","CD4.c20.Treg.TNFRSF9","CD4.c19.Treg.S1PR1","CD4.c18.Treg.RTKN2","CD4.c23.Mix.NME1","CD4.c24.Mix.NME2"))
#Treg<-subset(pbmc,idents = c("CD4.c21.Treg.OAS1","CD4.c20.Treg.TNFRSF9","CD4.c19.Treg.S1PR1","CD4.c18.Treg.RTKN2"))
Treg<-subset(pbmc,idents = c("3","8","20"))
#Idents(FOXP3)<-FOXP3@meta.data$loc
Idents(Treg)<-Treg@meta.data$meta.cluster
DimPlot(Treg)
Treg<-subset(Treg,idents = c("","CD4.c21.Treg.OAS1","CD4.c20.Treg.TNFRSF9","CD4.c19.Treg.S1PR1","CD4.c18.Treg.RTKN2"))
DimPlot(Treg)
FeaturePlot(Treg,features = "FOXP3")
Treg<-FindNeighbors(Treg,dims = 1:31)
Treg<-FindClusters(Treg,resolution = 1)
VlnPlot(Treg,features = "FOXP3",pt.size = 0,sort = TRUE)
FeaturePlot(Treg,features = "FOXP3",label = TRUE)
DimPlot(Treg)
Treg<-subset(Treg,idents = c("0","1","2","3","5","6","7","8","9","12","13","14"))
FeaturePlot(Treg,features = "FOXP3")
Treg<-FindNeighbors(Treg,dims = 1:31)
Treg<-FindClusters(Treg,resolution = 1)
VlnPlot(Treg,features = "FOXP3",pt.size = 0,sort = TRUE)
FeaturePlot(Treg,features = "FOXP3",label = TRUE)
Treg<-subset(Treg,idents = c("0","1","2","3","5","6","8","9","12"))
FeaturePlot(Treg,features = "FOXP3")
Treg<-FindNeighbors(Treg,dims = 1:31)
Treg<-FindClusters(Treg,resolution = 1)
VlnPlot(Treg,features = "FOXP3",pt.size = 0,sort = TRUE)


Idents(Treg)<-Treg@meta.data$subgroup
VlnPlot(Treg,features = "JMJD1C",sort = TRUE)
#VlnPlot(FOXP3,features = "JMJD1C",sort = TRUE,pt.size = 0)
#VlnPlot(Treg,features = "IFNG",sort = TRUE)
#VlnPlot(FOXP3,features = "IFNG",sort = TRUE,pt.size = 0)
a<-read.table("D:/GSE156728/STAT3.txt",sep = "\t",header = TRUE)
b<-read.table("D:/GSE156728/PI3K.txt",sep = "\t",header = TRUE)
a_new<-a[2:length(rownames(a)),]
b_new<-b[2:length(rownames(b)),]
Treg<-AddModuleScore(Treg,features = a_new,name = "STAT3")
Treg<-AddModuleScore(Treg,features = b_new,name = "PI3K")

Treg_T<-subset(Treg,idents = c("T","TIL"))
Idents(Treg_T)<-Treg_T@meta.data$group
patient<-unique(Treg_T@meta.data$group)
remove(a,b,a_new,b_new)
for (i in 1:length(patient)) {
  assign(paste("Treg_",i,"_TIL",sep = ""),subset(Treg_T,idents = patient[i]))
  assign("a",slot(get(paste("Treg_",i,"_TIL",sep = "")),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@scale.data))))
  assign("c",select(b,one_of("JMJD1C")))
  assign("d",slot(get(paste("Treg_",i,"_TIL",sep = "")),"meta.data"))
  STAT3<-select(d,starts_with("STAT3"))
  d$STAT3<-apply(STAT3,MARGIN = 1,median)
  PI3K<-select(d,starts_with("PI3K"))
  d$PI3K<-apply(PI3K,MARGIN = 1,median)
  assign("e",as.data.frame(d$PI3K))
  colnames(e)<-"PI3K"
  assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(mean(c$JMJD1C),mean(e$PI3K)))
  
}
scaldata<-scaldata_1_TIL
for (i in 2:length(patient)) {
  scaldata<-rbind(scaldata,get(paste("scaldata_",i,"_TIL",sep = "")))
}
colnames(scaldata)<-c("JMJD1C","PI3K_signaling")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=PI3K_signaling)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
write.table(scaldata,"Treg_1C_Treg_PI3Ksignalingmedian_scaledata_mean.txt",sep = "\t")
remove(a,b,c,d,scaldata,scaldata_1_TIL,scaldata_10_TIL,scaldata_11_TIL,scaldata_12_TIL,scaldata_13_TIL,scaldata_14_TIL,scaldata_15_TIL,scaldata_16_TIL,scaldata_17_TIL,scaldata_18_TIL,scaldata_19_TIL)
remove(scaldata_2_TIL,scaldata_20_TIL,scaldata_21_TIL,scaldata_22_TIL,scaldata_23_TIL,scaldata_24_TIL,scaldata_25_TIL,scaldata_26_TIL,scaldata_27_TIL,scaldata_28_TIL,scaldata_29_TIL,scaldata_3_TIL,scaldata_30_TIL)
remove(scaldata_31_TIL,scaldata_32_TIL,scaldata_33_TIL,scaldata_34_TIL,scaldata_35_TIL,scaldata_36_TIL,scaldata_37_TIL,scaldata_38_TIL,scaldata_39_TIL,scaldata_4_TIL,scaldata_40_TIL,scaldata_41_TIL,scaldata_42_TIL,scaldata_43_TIL,scaldata_44_TIL,scaldata_45_TIL,scaldata_46_TIL,scaldata_47_TIL,scaldata_48_TIL,scaldata_49_TIL,scaldata_5_TIL,scaldata_6_TIL,scaldata_7_TIL,scaldata_8_TIL,scaldata_9_TIL)
remove(scaldata_50_TIL,scaldata_51_TIL,scaldata_52_TIL,scaldata_53_TIL,scaldata_54_TIL,scaldata_55_TIL,scaldata_56_TIL,scaldata_57_TIL,scaldata_58_TIL,scaldata_59_TIL,scaldata_60_TIL,scaldata_61_TIL,scaldata_62_TIL,scaldata_63_TIL,scaldata_64_TIL,scaldata_65_TIL,scaldata_66_TIL,scaldata_67_TIL,scaldata_68_TIL,scaldata_69_TIL,scaldata_70_TIL,scaldata_71_TIL,scaldata_72_TIL,scaldata_73_TIL,Treg_1_TIL,Treg_10_TIL,Treg_11_TIL,Treg_12_TIL,Treg_13_TIL,Treg_14_TIL,Treg_15_TIL,Treg_16_TIL,Treg_17_TIL,Treg_18_TIL,Treg_19_TIL)
remove(Treg_2_TIL,Treg_3_TIL,Treg_4_TIL,Treg_5_TIL,Treg_6_TIL,Treg_7_TIL,Treg_8_TIL,Treg_9_TIL,Treg_20_TIL,Treg_21_TIL,Treg_22_TIL,Treg_23_TIL,Treg_24_TIL,Treg_25_TIL,Treg_26_TIL,Treg_27_TIL,Treg_28_TIL,Treg_29_TIL,Treg_30_TIL,Treg_31_TIL,Treg_32_TIL,Treg_33_TIL,Treg_34_TIL,Treg_35_TIL,Treg_36_TIL,Treg_37_TIL,Treg_38_TIL,Treg_39_TIL,Treg_40_TIL,Treg_41_TIL,Treg_42_TIL,Treg_43_TIL,Treg_44_TIL,Treg_45_TIL,Treg_46_TIL,Treg_47_TIL,Treg_48_TIL,Treg_49_TIL,Treg_50_TIL,Treg_51_TIL,Treg_52_TIL,Treg_53_TIL,Treg_54_TIL,Treg_55_TIL,Treg_56_TIL,Treg_57_TIL,Treg_58_TIL,Treg_59_TIL)
remove(Treg_60_TIL,Treg_61_TIL,Treg_62_TIL,Treg_63_TIL,Treg_64_TIL,Treg_65_TIL,Treg_66_TIL,Treg_67_TIL,Treg_68_TIL,Treg_69_TIL,Treg_70_TIL,Treg_71_TIL,Treg_72_TIL,Treg_73_TIL,i,all.genes,patient)
remove(e,PI3K,STAT3)
