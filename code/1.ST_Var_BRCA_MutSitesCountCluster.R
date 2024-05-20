#载入所需的R包；
library(Seurat)#
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r) #
###############################读取空转数据################################
data_dir<- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599"
list.files(data_dir) #查看目录中的文件；
file_name<- "filtered_feature_bc_matrix.h5"
brain<- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice = "anterior1")
brain@project.name <- "anterior1"
Idents(brain) <- "anterior1"
brain$orig.ident <- "anterior1"
for (i in colnames((brain@images$anterior1@coordinates))) {
  brain@images$anterior1@coordinates[[i]] <- as.integer(brain@images$anterior1@coordinates[[i]])
}
###############################数据的预处理################################
##数据标准化、识别高变基因,20秒
brain<- SCTransform(brain, assay = "Spatial", verbose = FALSE) #用SCTransform对数据进行标准化, 同时检测高变基因, 输出结果储存在 SCT assay中；#我们发现不同spot的mRNA分子数差异很大，特别是当组织中的细胞密度存在差异时，因此需要对数据进行标准化。由于细胞组织存在异质性，不同组织区域的细胞密度可能不同；因此如果采用常规单细胞转录组数据的LogNormalize标准化方法可能存在问题；这里Seurat作者推荐sctransform方法。
##展示高变基因
p5<- VariableFeaturePlot(brain,cols = c( "gray60", "red"))
top10<- head(VariableFeatures(brain),10) 
p6<- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0) 
#p6#解释图纵轴：残差方差（residual variance），将观测值与模型预测值之间的差异的平方和除以自由度。即基因在每个细胞中的表达与表达均值的平方和除以减一后的细胞数目
###############################降维、无监督聚类################################
brain<- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain<- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain<- FindClusters(brain, verbose = FALSE)
brain<- RunUMAP(brain, reduction = "pca", dims = 1:30)
#绘制UMAP分群图；#在切片图像上映射分群信息；
p9 <- DimPlot(brain, reduction = "umap", label = TRUE,label.size = 12)+theme(text = element_text(size=24))
p10 <- SpatialDimPlot(brain, label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 0)+ 
  theme(text = element_text(size=24))
options(repr.plot.width=16, repr.plot.height=8)#让图片展示变大
SpatialDimPlot(brain,alpha =0) + NoLegend() + p9 + p10

saveRDS(brain, "E:/数据库/STMut/数据处理/中间文件/ST_seuratObject.rds")
###########################################################################
##########################先跑出p10以前的代码，然后##########################
###########################################################################
##########################基因突变位点数突变谱聚类###############################
############################################################################
#参考 GSE167036_Cell_subsets_definedby_mutGenes.R
file_path <- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/Step4_VariantCalling"
##########################0-读取scsnv注释矩阵###############################
anno<- read.table(paste0(file_path,"/GSM6177599.variants.annovar.Step4.3.1.hg38_multianno.csv"), header=T,sep=",",as.is=T,fill=T)
anno<-anno[,c(1,2,7)]
anno<- anno[nchar(anno[,3]) > 1,]# 保留第三列字符长度大于1的行
anno<-anno[!grepl("\\|", anno[,3]), ]
library(tidyr)
anno2 <- anno %>% as_tibble() %>% 
  separate_rows(Gene.refGene, sep = ";")

file_list <- list.files(path = file_path, pattern = "genotype.tsv$", full.names = TRUE)
spot_gene_mutcount<-data.frame()
##########################1-读取scsnv突变矩阵###############################
for (file in file_list) {
  df2<- read.table(file, header=T,sep="\t",as.is=T,comment.char = "")
  df2 <- df2[which(df2$ALT_expected!="."),]
  df2 <- df2[which(df2$Base_observed!=df2$REF),]
  colnames(df2)[1]<-"Chr"
  df2<-merge(df2,anno2)
  ##spot的基因突变次数
  df3<-df2 %>%
    dplyr::select(CB,Chr,Start,Gene.refGene) %>%
    dplyr::distinct() %>%
    dplyr::count(CB,Gene.refGene)
  spot_gene_mutcount <- rbind(spot_gene_mutcount,df3)
}
library(reshape2)
##########################2-scsnv突变矩阵预处理###############################
spot_gene_mutcount_matrix <- dcast(spot_gene_mutcount, CB ~ Gene.refGene, value.var = "n")
##参考单细胞处理流程https://www.jianshu.com/p/728125be7c53
##参考单细胞处理流程代码GSE162708_infercnv.R
spot_gene_mutcount_matrix[,1] <- paste0(spot_gene_mutcount_matrix[,1],"-1")
#基于空转整合spot
spot_gene_mutcount_matrix <- merge(data.frame(CB = rownames(brain@meta.data)), spot_gene_mutcount_matrix, all.x = T)#spot从2298增加到2340
rownames(spot_gene_mutcount_matrix) <- spot_gene_mutcount_matrix[,1]
spot_gene_mutcount_matrix[is.na(spot_gene_mutcount_matrix)] <- 0

write.table(spot_gene_mutcount_matrix, "E:/数据库/STMut/数据处理/中间文件/spot_gene_mutcount_matrix.txt", sep = "\t", quote = F)

spot_gene_mutcount_matrix[, -1][spot_gene_mutcount_matrix[, -1] > 0] <- 1

write.table(spot_gene_mutcount_matrix, "E:/数据库/STMut/数据处理/中间文件/spot_gene_01_matrix.txt", sep = "\t", quote = F)
##质量控制
###不过滤变异的细胞数目非常少的基因
#gene_count <- rowSums(spot_gene_mutcount_matrix[, -1] > 0)# 计算每个基因在多少个细胞中表达
#filtered_gene_data <- spot_gene_mutcount_matrix[gene_count >= 100, ]# 筛选出在100个或更多细胞中表达的基因
seurat_object <- CreateSeuratObject(counts = t(spot_gene_mutcount_matrix[,-1]))
VlnPlot(object = seurat_object, features = c("nCount_RNA", "nFeature_RNA"))
#不打算过滤变异数目非常多或者非常少的细胞,因为表达过滤了低质量的细胞
#seurat_object <- subset(seurat_object,subset=nFeature_RNA>10 & nFeature_RNA<300 )
# 可选：进行预处理和标准化
seurat_object <- NormalizeData(seurat_object, normalization.method = 'LogNormalize', scale.factor = 100)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_object),10)#画出不带标签或带标签高变突变基因点图
plot1 <- VariableFeaturePlot(seurat_object)
top10
plot2 <- LabelPoints(plot = plot1,points = top10,repel = TRUE)
plot2
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
##########################3-scsnv突变矩阵降维聚类###############################
#6. 降维聚类
#（基于前面得到的high variable基因的scale矩阵）
seurat_object <- RunPCA(seurat_object, npcs = 50, verbose = FALSE)
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

seurat_object <- RunUMAP(seurat_object, dims = 1:30)
#seurat_object <- RunTSNE(seurat_object, dims = 1:30)
DimPlot(seurat_object , reduction = "umap", pt.size=0.5)

#输出01突变谱的seurat对象
saveRDS(seurat_object, "E:/数据库/STMut/数据处理/中间文件/Mut01_seuratObject.rds")

#输出位点数突变谱的seurat对象
saveRDS(seurat_object, "E:/数据库/STMut/数据处理/中间文件/MutSitesCount_seuratObject.rds")


##########################4-空转展示scsnv突变聚类类别###############################
#***********************************************************************************
#
seurat_object_1 <- readRDS("E:/数据库/STMut/数据处理/中间文件/Mut01_seuratObject.rds")
seurat_object_2 <- readRDS("E:/数据库/STMut/数据处理/中间文件/MutSitesCount_seuratObject.rds")
brain <- readRDS("E:/数据库/STMut/数据处理/中间文件/ST_seuratObject.rds")

head(seurat_object_1@meta.data)
head(seurat_object_2@meta.data)
aa <- seurat_object_1@meta.data
bb <- seurat_object_2@meta.data
cc <- brain@meta.data
#***********************************************************************************

##取出空转metadata的行名构成数据框
ST_rownames <- data.frame(CB = rownames(brain@meta.data))
mut_cluster <- data.frame(CB = rownames(seurat_object_1@meta.data), seurat_object_1@meta.data[,c(4,5)])
colnames(mut_cluster)[ncol(mut_cluster)] <- "mut_cluster"
mut_cluster <- merge(ST_rownames, mut_cluster, all.x=T)
mut_cluster$mut_cluster <- as.character(mut_cluster$mut_cluster)
mut_cluster[is.na(mut_cluster)]<-"WT"
rownames(mut_cluster) <- mut_cluster$CB
brain@meta.data <- cbind(brain@meta.data, mut_cluster[,-1])

nCount_nFeatureSNV <- data.frame(CB = rownames(seurat_object_2@meta.data), seurat_object_2@meta.data[,c(2,3)])
brain@meta.data <- cbind(brain@meta.data, nCount_nFeatureSNV[,-1])

saveRDS(brain, "E:/数据库/STMut/数据处理/中间文件/ST_Mut_SeuratObject.rds")

SpatialDimPlot(brain, group.by="mut_cluster",label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 0)+ 
  theme(text = element_text(size=24))


##########################5-空转展示spot突变位点数###############################
# SpatialFeaturePlot(brain, features="nCount_RNA",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 0)+ 
#   theme(text = element_text(size=24))+theme(legend.position = "right")
# ##########################6-空转展示spot突变基因数###############################
# SpatialFeaturePlot(brain, features="nFeature_RNA",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 0)+ 
#   theme(text = element_text(size=24))+theme(legend.position = "right")
# ##########################6.0-空转展示spot的mRNA分子数###############################
# SpatialFeaturePlot(brain, features="nCount_Spatial",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 0)+ 
#   theme(text = element_text(size=24))+theme(legend.position = "right")
# ##########################6.1-空转展示spot的mRNA分子数###############################
# SpatialFeaturePlot(brain, features="nFeature_Spatial",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 0)+ 
#   theme(text = element_text(size=24))+theme(legend.position = "right")
# ##########################7-空转展示scsnv突变基因###############################
# spot_gene_mutcount_matrix_1<-merge(ST_rownames,spot_gene_mutcount_matrix,all.x=T)
# spot_gene_mutcount_matrix_1[is.na(spot_gene_mutcount_matrix_1)]<-0
# brain@meta.data<-cbind(brain@meta.data,spot_gene_mutcount_matrix_1[,-1])
# SpatialFeaturePlot(brain, features = top10)#展示空转基因突变计数

#********************************************************
brain <- readRDS("E:/数据库/STMut/数据处理/中间文件/ST_Mut_SeuratObject.rds")

#111111这里导出小提琴图，放到数据库里面
nCount_Spatial <- VlnPlot(object = brain, features = c("nCount_Spatial"), group.by = "orig.ident")+theme(legend.position = "none",axis.text.x = element_blank())
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_Spatial_viloin.png", width = 337, height = 280)
# print(nCount_Spatial)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_Spatial_viloin_blow.png", width = 674, height = 674)
VlnPlot(object = brain, features = c("nCount_Spatial"), group.by = "orig.ident")+theme(legend.position = "none",axis.text.x = element_blank())
dev.off()

nFeature_Spatial <- VlnPlot(object = brain, features = c("nFeature_Spatial"), group.by = "orig.ident")+theme(legend.position = "none",axis.text.x = element_blank())
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_Spatial_viloin.png", width = 337, height = 280)
# print(nFeature_Spatial)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_Spatial_viloin_blow.png", width = 674, height = 674)
VlnPlot(object = brain, features = c("nFeature_Spatial"), group.by = "orig.ident")+theme(legend.position = "none",axis.text.x = element_blank())
dev.off()

#这里导出小提琴图，放到数据库里面
nCount_RNA <- VlnPlot(object = brain, features = c("nCount_RNA"), group.by = "orig.ident")+theme(legend.position = "none",axis.text.x = element_blank())+labs(title="nCount_SNV")
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_RNA_viloin.png", width = 337, height = 280)
# print(nCount_RNA)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_SNV_viloin_blow.png", width = 674, height = 674)
print(nCount_RNA)
# 进行绘图操作
dev.off()

nFeature_RNA <- VlnPlot(object = brain, features = c("nFeature_RNA"), group.by = "orig.ident")+theme(legend.position = "none",axis.text.x = element_blank())+labs(title="nFeature_SNV")
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_RNA_viloin.png", width = 337, height = 280)
# print(nFeature_RNA)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_SNV_viloin_blow.png", width = 674, height = 674)
print(nFeature_RNA)
# 进行绘图操作
dev.off()

#222222这里导出nCount_RNA、nFeature_RNA、nCount_Spatial、nFeature_Spatial映射的切片，放到数据库里
#展示spot突变位点数
nCount_RNA <- SpatialFeaturePlot(brain, features="nCount_RNA",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 1)+ 
  theme(text = element_text(size=8))+theme(legend.position = "right")+labs(fill="nCount_SNV")
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_RNA_slice.png", width = 337, height = 280)
# print(nCount_RNA)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_SNV_slice_blow.png", width = 674, height = 674)
print(nCount_RNA)
# 进行绘图操作
dev.off()

#展示spot突变基因数
nFeature_RNA <- SpatialFeaturePlot(brain, features="nFeature_RNA",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 1)+ 
  theme(text = element_text(size=8))+theme(legend.position = "right")+labs(fill="nFeature_SNV")
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_RNA_slice.png", width = 337, height = 280)
# print(nFeature_RNA)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_SNV_slice_blow.png", width = 674, height = 674)
print(nFeature_RNA)
# 进行绘图操作
dev.off()
#展示spot的mRNA分子数
nCount_Spatial <- SpatialFeaturePlot(brain, features="nCount_Spatial",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 1)+ 
  theme(text = element_text(size=8))+theme(legend.position = "right")
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_Spatial_slice.png", width = 337, height = 280)
# print(nCount_Spatial)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nCount_Spatial_slice_blow.png", width = 674, height = 674)
print(nCount_Spatial)
# 进行绘图操作
dev.off()
#展示spot的基因数目
nFeature_Spatial <- SpatialFeaturePlot(brain, features="nFeature_Spatial",pt.size.factor = 3,alpha = c(0.8, 1),image.alpha = 1)+ 
  theme(text = element_text(size=8))+theme(legend.position = "right")
# png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_Spatial_slice.png", width = 337, height = 280)
# print(nFeature_Spatial)
# # 进行绘图操作
# dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/nFeature_Spatial_slice_blow.png", width = 674, height = 674)
print(nFeature_Spatial)
# 进行绘图操作
dev.off()

#这里导出cluster图---------------------------------------------------------------------------------------------------------------------
seurat_object <- seurat_object_1
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/cluster_blow.png", width = 674, height = 674)
DimPlot(seurat_object, reduction = "umap", pt.size=0.5)
dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/cluster_slice_blow.png", width = 674, height = 674)
SpatialDimPlot(brain, group.by="mut_cluster",label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 1)+ 
  theme(text = element_text(size=24))
dev.off()

#********************************************************
#*获取切片图像坐标和key
key(object = brain@images$anterior1)
head(GetTissueCoordinates(brain))
head(Idents(brain))

Slice_Coor_cluster <- brain@meta.data[,c(2,3,4,5,9,10,11)]
spot_corrdinates <- GetTissueCoordinates(brain)
spot_corrdinates$CB <- rownames(spot_corrdinates)

Slice_Coor_Cluster_integrate <- cbind(Slice_Coor_cluster, spot_corrdinates)

#突变谱
library(tidyverse)
spot_gene_mutcount_matrix <- read.table("E:/数据库/STMut/数据处理/中间文件/spot_gene_mutcount_matrix.txt", sep = "\t", row.names = 1)
long_data <- spot_gene_mutcount_matrix %>% 
             gather(key = "Gene", value = "MutCount", -CB)

# 定义一个函数，将多列值按逗号分隔开
combine_values <- function(x) {
  paste(x, collapse = ",")
}
combine_values_1 <- function(x) {
  paste(x, collapse = ";")
}

#保证CB的顺序一致，所以先排序
Slice_Coor_Cluster_integrate <- Slice_Coor_Cluster_integrate[order(Slice_Coor_Cluster_integrate$mut_cluster, Slice_Coor_Cluster_integrate$CB), ]
# 按照group分组，将value1、value2、value3按逗号分隔开，不同分组用分号隔开
# 将第一列group的值也按逗号分隔开
Slice_Coor_Cluster_integrate$mut_cluster_1 <- Slice_Coor_Cluster_integrate$mut_cluster
Slice_Coor_Cluster_integrate$mut_cluster_1 <- as.character(Slice_Coor_Cluster_integrate$mut_cluster_1)
result_1 <- aggregate(. ~ mut_cluster_1, data = Slice_Coor_Cluster_integrate, FUN = combine_values)

result_2 <- data.frame(
  nCount_Spatial = paste(result_1$nCount_Spatial, collapse = ";"),
  nFeature_Spatial = paste(result_1$nFeature_Spatial, collapse = ";"),
  nCount_SCT = paste(result_1$nCount_SCT, collapse = ";"),
  nFeature_SCT = paste(result_1$nFeature_SCT, collapse = ";"),
  nCount_RNA = paste(result_1$nCount_RNA, collapse = ";"),
  nFeature_RNA = paste(result_1$nFeature_RNA, collapse = ";"),
  imagerow = paste(result_1$imagerow, collapse = ";"),
  imagecol = paste(result_1$imagecol, collapse = ";"),
  CB = paste(result_1$CB, collapse = ";"),
  mut_cluster = paste(result_1$mut_cluster, collapse = ";")
)
result_3 <- as.data.frame(t(result_2))
colnames(result_3)[1] <- "Column_value"
result_3$Column_name <- rownames(result_3)

mut_cluster <- data.frame(CB = rownames(brain@meta.data), mut_cluster = brain@meta.data$mut_cluster) 
long_data_1 <- merge(long_data, mut_cluster, all = T)
long_data_1$mut_cluster_1 <- long_data_1$mut_cluster
long_data_2 <- long_data_1[order(long_data_1$mut_cluster_1, long_data_1$CB), ]
long_data_2$mut_cluster_1 <- as.character(long_data_2$mut_cluster_1)
#按照gene和mut_cluster分组，将其他列的值用逗号分隔粘贴起来
long_data_3 <- long_data_2 %>%
  group_by(Gene, mut_cluster_1) %>%
  summarize(CB = paste(CB, collapse = ","),
            mut_cluster = paste(mut_cluster, collapse = ","),
            MutCount = paste(MutCount, collapse = ","))
long_data_3 <- long_data_3[,-2]
long_data_4 <- aggregate(. ~ Gene, data = long_data_3, FUN = combine_values_1)
long_data_5 <- long_data_4[,c(1,4)]
colnames(long_data_5) <- c("Column_name", "Column_value")

result_4 <- rbind(result_3, long_data_5)
write.table(result_4, "E:/数据库/STMut/数据处理/结果表格/Slice_Coor_Cluster_Mut.txt", row.names = F, quote = F, sep = "\t")
#********************************************************
###########################################################################
##########################先跑出p10以前的代码，然后##########################
###########################################################################
##########################基因突变位点数突变谱过滤后聚类###############################
############################################################################
#参考 GSE167036_Cell_subsets_definedby_mutGenes.R
# file_path <- "F:/Desktop/ST_data/spatial_variation/data/melanoma"
# ##########################0-读取scsnv注释矩阵###############################
# anno<- read.table(paste0(file_path,"/GSM5420750.variants.annovar.Step4.3.1.hg38_multianno.csv"), header=T,sep=",",as.is=T,fill=T)
# anno<-anno[,c(1,2,7)]
# anno<- anno[nchar(anno[,3]) > 1,]# 保留第三列字符长度大于1的行
# anno<-anno[!grepl("\\|", anno[,3]), ]
# library(tidyr)
# anno2 <-anno %>% as_tibble() %>% 
#   separate_rows(Gene.refGene, sep = ";")
# 
# file_list <- list.files(path = file_path, pattern = "genotype.tsv$", full.names = TRUE)
# spot_gene_mutcount<-data.frame()
# ##########################1-读取scsnv突变矩阵###############################
# for (file in file_list) {
#   df2<- read.table(file, header=T,sep="\t",as.is=T,comment.char = "")
#   df2 <- df2[which(df2$ALT_expected!="."),]
#   colnames(df2)[1]<-"Chr"
#   df2<-merge(df2,anno2)
#   ##spot的基因突变次数
#   df3<-df2 %>%
#     dplyr::select(CB,Chr,Start,Gene.refGene) %>%
#     dplyr::distinct() %>%
#     dplyr::count(CB,Gene.refGene)
#   spot_gene_mutcount<-rbind(spot_gene_mutcount,df3)
# }
# library(reshape2)
# ##########################2-scsnv突变矩阵预处理###############################
# spot_gene_mutcount_matrix<-dcast(spot_gene_mutcount, CB ~ Gene.refGene, value.var = "n")
# ##参考单细胞处理流程https://www.jianshu.com/p/728125be7c53
# ##参考单细胞处理流程代码GSE162708_infercnv.R
# spot_gene_mutcount_matrix[,1]<-paste0(spot_gene_mutcount_matrix[,1],"-1")
# #基于空转整合spot
# spot_gene_mutcount_matrix<-merge(data.frame(CB=rownames(brain@meta.data)),spot_gene_mutcount_matrix,all.x=T)
# rownames(spot_gene_mutcount_matrix)<-spot_gene_mutcount_matrix[,1]
# spot_gene_mutcount_matrix[is.na(spot_gene_mutcount_matrix)]<-0
# ##质量控制
# ###过滤变异的细胞数目非常少的基因
# gene_count <- rowSums(spot_gene_mutcount_matrix[, -1] > 0)# 计算每个基因在多少个细胞中表达
# filtered_gene_data <- spot_gene_mutcount_matrix[gene_count >= 100, ]# 筛选出在100个或更多细胞中表达的基因
# seurat_object <- CreateSeuratObject(counts = t(filtered_gene_data[,-1]))
# VlnPlot(object = seurat_object, features = c("nCount_RNA", "nFeature_RNA"))
# #不打算过滤变异数目非常多或者非常少的细胞,因为表达过滤了低质量的细胞
# #seurat_object <- subset(seurat_object,subset=nFeature_RNA>10 & nFeature_RNA<300 )
# # 可选：进行预处理和标准化
# seurat_object <- NormalizeData(seurat_object,normalization.method = 'LogNormalize', scale.factor = 100)
# seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# top10 <- head(VariableFeatures(seurat_object),10)#画出不带标签或带标签高变突变基因点图
# plot1 <- VariableFeaturePlot(seurat_object)
# top10
# plot2 <- LabelPoints(plot = plot1,points = top10,repel = TRUE)
# plot2
# seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
# ##########################3-scsnv突变矩阵降维聚类###############################
# #6. 降维聚类
# #（基于前面得到的high variable基因的scale矩阵）
# seurat_object <- RunPCA(seurat_object, npcs = 50, verbose = FALSE)
# seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
# seurat_object <- FindClusters(seurat_object, resolution = 0.5)
# 
# seurat_object <- RunUMAP(seurat_object, dims = 1:30)
# seurat_object <- RunTSNE(seurat_object, dims = 1:30)
# DimPlot(seurat_object , reduction = "umap", pt.size=0.5)
# ##########################4-空转展示scsnv突变聚类类别###############################
# ##取出空转metadata的行名构成数据框
# ST_rownames <- data.frame(CB=rownames(brain@meta.data))
# mut_cluster <- data.frame(CB=rownames(seurat_object@meta.data),mut_cluster=seurat_object@meta.data$seurat_clusters)
# mut_cluster <- merge(ST_rownames,mut_cluster,all.x=T)
# mut_cluster$mut_cluster<-as.character(mut_cluster$mut_cluster)
# mut_cluster[is.na(mut_cluster)] <- "WT"
# brain@meta.data <- cbind(brain@meta.data, mut_cluster[,-1])
# colnames(brain@meta.data)[ncol(brain@meta.data)] <- "mut_cluster"
# SpatialDimPlot(brain, group.by="mut_cluster",label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 0)+ 
#   theme(text = element_text(size=24))
# ##########################5-空转展示scsnv突变基因###############################
# spot_gene_mutcount_matrix_1 <- merge(ST_rownames, spot_gene_mutcount_matrix, all.x=T)
# spot_gene_mutcount_matrix_1[is.na(spot_gene_mutcount_matrix_1)] <- 0
# brain@meta.data <- cbind(brain@meta.data, spot_gene_mutcount_matrix_1[,-1])
# SpatialFeaturePlot(brain, features = top10)#展示空转基因突变计数
###########################################################################
##########################先跑出p10以前的代码，然后##########################
###########################################################################

#************************************************************************************************************************************
##利用 DoHeatmap 命令可以可视化marker基因的表达
#min.pct:基因在最少多大比例的细胞中检测到才用于后续的检测分析；logfc.threshold:高于给定的最小的变化倍数（log2）的基因用于后续统计分析。
#数值越大，计算越快，差异基因可能会少；seurat v3中计算的差异倍数是loge,也就是ln
seurat_object <- seurat_object_1

ST_rownames <- data.frame(CB = rownames(seurat_object_2@meta.data))
mut_cluster <- data.frame(CB = rownames(seurat_object_2@meta.data), seurat_object_2@meta.data[,5])
colnames(mut_cluster)[ncol(mut_cluster)] <- "mut_cluster"
mut_cluster <- merge(ST_rownames, mut_cluster, all.x=T)
mut_cluster$mut_cluster <- as.character(mut_cluster$mut_cluster)
mut_cluster[is.na(mut_cluster)]<-"WT"
rownames(mut_cluster) <- mut_cluster$CB
seurat_object@meta.data <- cbind(seurat_object@meta.data, mut_cluster[,-1])

View(seurat_object@meta.data)
View(seurat_object_2@meta.data)

saveRDS(seurat_object, "E:/数据库/STMut/数据处理/中间文件/MutSites_01_SeuratObject.rds")

seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25)
seurat_object.markers[,c(1,2,5)] <- round(seurat_object.markers[,c(1,2,5)], 4)
seurat_object.markers$disease <- "BRCA"
seurat_object.markers$dataset <- "GSE203612"
seurat_object.markers$sample <- "GSM6177599"
write.table(seurat_object.markers, "E:/数据库/STMut/数据处理/结果表格/BRCA_GSE203612_markers.txt", sep = "\t", row.names = F, quote = F)
#查看每个亚群的前两个
library(tidyverse)
seurat_object.markers %>% 
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
#每个细胞群top 5 marker gene 热图
top5 <- seurat_object.markers %>% 
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/top5marker_heatmap.png", width = 337, height = 280)
DoHeatmap(seurat_object, features = top5$gene) + NoLegend()
dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/top5marker_heatmap_blow.png", width = 674, height = 674)
DoHeatmap(seurat_object, features = top5$gene) + NoLegend()
dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/top5marker_dotplot.png", width = 337, height = 280)
DotPlot(seurat_object, features = unique(top5$gene), scale = F) + RotatedAxis() + guides(shape=guide_legend(title = NULL),
                                                                                         color=guide_legend(title = 'Average Mutation'),
                                                                                         size=guide_legend(title = 'Percent Mutated'),
)
dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/top5marker_dotplot_blow.png", width = 674, height = 674)
DotPlot(seurat_object, features = unique(top5$gene), scale = F) + RotatedAxis() + guides(shape=guide_legend(title = NULL),
                                                                                         color=guide_legend(title = 'Average Mutation'),
                                                                                         size=guide_legend(title = 'Percent Mutated'),
)
dev.off()

#功能富集
x1 <- seurat_object.markers[which(seurat_object.markers$p_val<0.05 & seurat_object.markers$avg_log2FC>0.2),7]

library(org.Hs.eg.db) #人类注释数据???
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(tidyverse)
library(GOplot)

#GO
erich.go.BP = enrichGO(gene = x1,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/GO_enrichment.png", width = 1200, height = 674)
dotplot(erich.go.BP, showCategory = 10)
dev.off()

# erich.go.MF = enrichGO(gene = x1,
#                        OrgDb = org.Hs.eg.db,
#                        keyType = "SYMBOL",
#                        ont = "MF",
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff = 0.05)
# dotplot(erich.go.MF, showCategory = 10)

write.csv(go,file="E:/数据库/STMut/数据处理/结果表格/go_0.05_ENTREZID_0.02.csv")

#KEGG
#富集没有结果
#--> No gene can be mapped....
#--> Expected input gene ID: 
#--> return NULL...
hg <- bitr(x1, fromType="SYMBOL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")#ID转换
kegg <- enrichKEGG(hg$ENTREZID, 
                   organism = 'hsa', 
                   keyType = 'kegg', 
                   pvalueCutoff = 1, 
                   pAdjustMethod = 'BH', 
                   minGSSize = 1, 
                   maxGSSize = 500, 
                   qvalueCutoff = 1, 
                   use_internal_data = FALSE)#进行KEGG富集
dotplot(kegg, showCategory = 10)

write.csv(kegg,file="E:\\lcw\\RNA_editing\\result\\kegg_P_0.5_ENTREZID_0.02.csv")
#***********************************************************************************

##########################7-空转展示scsnv突变基因###############################
# spot_gene_mutcount_matrix_1 <- merge(ST_rownames, spot_gene_mutcount_matrix, all.x=T)
# spot_gene_mutcount_matrix_1[is.na(spot_gene_mutcount_matrix_1)] <- 0
# brain@meta.data <- cbind(brain@meta.data,spot_gene_mutcount_matrix_1[,-1])
# SpatialFeaturePlot(brain, features = top10)#展示空转基因突变计数
