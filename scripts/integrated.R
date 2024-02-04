library(Seurat)
setwd('/share2/pub/zhangyr/zhangyr/single-RNA/R_proj/NSCLC')



## 2021_NC
p2 = Seurat::Read10X('./data/external/scRNA/2021_NC/P2_Del19/')
colnames(p2) = paste0(colnames(p2),'_p2')
p2 = CreateSeuratObject(p2)

p3 = Seurat::Read10X('./data/external/scRNA/2021_NC/P3_L858R/')
colnames(p3) = paste0(colnames(p3),'_p3')
p3 = CreateSeuratObject(p3)

p7 = Seurat::Read10X('./data/external/scRNA/2021_NC/P7_Del19/')
colnames(p7) = paste0(colnames(p7),'_p7')
p7 = CreateSeuratObject(p7)

NC = merge(x = p2,y = c(p3,p7))
print(dim(NC))
NC[["percent.mt"]] = PercentageFeatureSet(NC, pattern = "^mt-")
NC = subset(NC, subset = nFeature_RNA > 200 & percent.mt < 20)
print(dim(NC))

## 2022_JIC
JIC = readRDS('./results/Step5_validation/data/2022_JIC/2022_JIC_stageI.rds')

## our
data = SeuratDisk::LoadH5Seurat('./data/Shanghai_Lung_merged.v1.03.h5seurat')
data = data[,data$Pathology == 'LUAD']
data[['nCount_RNA']] = colSums(x = data,slot = 'counts')
data[['nFeature_RNA']] = colSums(x = GetAssayData(data,slot = 'counts')>0)
data[['percent.mt']] = PercentageFeatureSet(data,pattern = '^MT-')
data = subset(data, subset = nFeature_RNA>200 & percent.mt < 20)

library(stringr)

meta = data@meta.data
meta_temp = meta
meta_temp$Stage = str_trim(meta_temp$Stage)
meta_temp$Mutation = str_trim(meta_temp$Mutation)
length(unique(meta_temp$Sample_Id))
remove_sample = c()
#remove_sample = c(remove_sample,as.character(unique(meta_temp$Sample_Id[which(meta_temp$Stage == 'IIIA')])))
remove_sample = c(remove_sample,as.character(unique(meta_temp$Sample_Id[which(meta_temp$Smoking == 1)])))
remove_sample = c(remove_sample,as.character(unique(meta_temp$Sample_Id[which(meta_temp$Mutation == '19 790')])))
remove_sample = c(remove_sample,as.character(unique(meta_temp$Sample_Id[which(meta_temp$Mutation == 'L861Q')])))
length(setdiff(unique(meta_temp$Sample_Id),remove_sample))
meta_temp = meta_temp[which(meta$Sample_Id %in% setdiff(unique(meta_temp$Sample_Id),remove_sample)),]
meta_temp$Sample_Id = as.character(meta_temp$Sample_Id)
data = data[,rownames(meta_temp)]
data@meta.data = meta_temp
print(dim(data))


merged = merge(x=data, y = (JIC,NC)))
print(dim(merged))