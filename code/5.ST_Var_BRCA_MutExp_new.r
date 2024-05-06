###########################################################################
##########################获得簇间差异突变基因##########################
###########################################################################
###########################################################################
############################################################################
#载入所需的R包；
library(Seurat)#
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r) #
###############################读取空转数据################################
# data_dir<- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599"
# list.files(data_dir) #查看目录中的文件；
# file_name<- "filtered_feature_bc_matrix.h5"
# brain<- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice = "anterior1")
# brain@project.name <- "anterior1"
# Idents(brain) <- "anterior1"
# brain$orig.ident <- "anterior1"
# for (i in colnames((brain@images$anterior1@coordinates))) {
#   brain@images$anterior1@coordinates[[i]] <- as.integer(brain@images$anterior1@coordinates[[i]])
# }
# # ###############################数据的预处理################################
# ##数据标准化、识别高变基因,20秒
# brain<- SCTransform(brain, assay = "Spatial", verbose = FALSE) #用SCTransform对数据进行标准化, 同时检测高变基因, 输出结果储存在 SCT assay中；#我们发现不同spot的mRNA分子数差异很大，特别是当组织中的细胞密度存在差异时，因此需要对数据进行标准化。由于细胞组织存在异质性，不同组织区域的细胞密度可能不同；因此如果采用常规单细胞转录组数据的LogNormalize标准化方法可能存在问题；这里Seurat作者推荐sctransform方法。
# ##展示高变基因
# p5<- VariableFeaturePlot(brain,cols = c( "gray60", "red"))
# top10<- head(VariableFeatures(brain),10)
# p6<- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0)
# p6#解释图纵轴：残差方差（residual variance），将观测值与模型预测值之间的差异的平方和除以自由度。即基因在每个细胞中的表达与表达均值的平方和除以减一后的细胞数目
# ###############################降维、无监督聚类################################
# brain<- RunPCA(brain, assay = "SCT", verbose = FALSE)
# brain<- FindNeighbors(brain, reduction = "pca", dims = 1:30)
# brain<- FindClusters(brain, verbose = FALSE)
# brain<- RunUMAP(brain, reduction = "pca", dims = 1:30)
# #绘制UMAP分群图；#在切片图像上映射分群信息；
# p9<- DimPlot(brain, reduction = "umap", label = TRUE,label.size = 12)+theme(text = element_text(size=24))
# p10<- SpatialDimPlot(brain, label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 0)+
#   theme(text = element_text(size=24))
# options(repr.plot.width=16, repr.plot.height=8)#让图片展示变大
# SpatialDimPlot(brain,alpha =0)+NoLegend()+p9+p10
# # 
# # ###########################################################################
# # ##########################先跑出p10以前的代码，然后##########################
# # ###########################################################################
# # ##########################基因突变0|1突变谱聚类###############################
# # ############################################################################
# # #参考 GSE167036_Cell_subsets_definedby_mutGenes.R
# file_path <- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/Step4_VariantCalling"
# # ##########################0-读取scsnv注释矩阵###############################
# anno<- read.table(paste0(file_path,"/GSM6177599.variants.annovar.Step4.3.1.hg38_multianno.csv"), header=T,sep=",",as.is=T,fill=T)
# anno<-anno[,c(1,2,7)]
# anno<- anno[nchar(anno[,3]) > 1,]# 保留第三列字符长度大于1的行
# anno<-anno[!grepl("\\|", anno[,3]), ]
# library(tidyr)
# anno2 <-anno %>% as_tibble() %>%
#   separate_rows(Gene.refGene, sep = ";")
# file_list <- list.files(path = file_path, pattern = "genotype.tsv$", full.names = TRUE)
# spot_gene_mutcount<-data.frame()
# # ##########################1-读取scsnv突变矩阵###############################
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
# # ##########################2-scsnv突变矩阵预处理###############################
# spot_gene_mutcount_matrix<-dcast(spot_gene_mutcount, CB ~ Gene.refGene, value.var = "n")
# ##参考单细胞处理流程https://www.jianshu.com/p/728125be7c53
# ##参考单细胞处理流程代码GSE162708_infercnv.R
# spot_gene_mutcount_matrix[,1]<-paste0(spot_gene_mutcount_matrix[,1],"-1")
# #基于空转整合spot
# spot_gene_mutcount_matrix<-merge(data.frame(CB=rownames(brain@meta.data)),spot_gene_mutcount_matrix,all.x=T)
# rownames(spot_gene_mutcount_matrix)<-spot_gene_mutcount_matrix[,1]
# spot_gene_mutcount_matrix[is.na(spot_gene_mutcount_matrix)]<-0
# spot_gene_mutcount_matrix[, -1][spot_gene_mutcount_matrix[, -1] > 0] <- 1
# ##质量控制
# ###不过滤变异的细胞数目非常少的基因
# #gene_count <- rowSums(spot_gene_mutcount_matrix[, -1] > 0)# 计算每个基因在多少个细胞中表达
# #filtered_gene_data <- spot_gene_mutcount_matrix[gene_count >= 100, ]# 筛选出在100个或更多细胞中表达的基因
# seurat_object <- CreateSeuratObject(counts = t(spot_gene_mutcount_matrix[,-1]))
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
# ST_rownames<-data.frame(CB=rownames(brain@meta.data))
# mut_cluster<-data.frame(CB=rownames(seurat_object@meta.data),seurat_object@meta.data[,-1])
# colnames(mut_cluster)[ncol(mut_cluster)]<-"mut_cluster"
# mut_cluster<-merge(ST_rownames,mut_cluster,all.x=T)
# mut_cluster$mut_cluster<-as.character(mut_cluster$mut_cluster)
# mut_cluster[is.na(mut_cluster)]<-"WT"
# brain@meta.data<-cbind(brain@meta.data,mut_cluster[,-1])
# SpatialDimPlot(brain, group.by="mut_cluster",label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 0)+
#   theme(text = element_text(size=24))
###############################5-spot类别差异变异基因###############################

brain <- readRDS("E:/数据库/STMut/数据处理/中间文件/ST_Mut_SeuratObject.rds")
spot_gene_mutcount <- read.table("E:/数据库/STMut/数据处理/中间文件/spot_gene_mutcount_matrix.txt", sep = "\t", header = T)
spot_gene_mutcount_matrix_raw <- spot_gene_mutcount

###基因突变位点数突变谱下###
library(reshape2)
# spot_gene_mutcount_matrix_raw <- dcast(spot_gene_mutcount, CB ~ Gene.refGene, value.var = "n")
# spot_gene_mutcount_matrix_raw[,1] <- paste0(spot_gene_mutcount_matrix_raw[,1],"-1")
# spot_gene_mutcount_matrix_raw <- merge(data.frame(CB=rownames(brain@meta.data)),spot_gene_mutcount_matrix_raw,all.x=T)
rownames(spot_gene_mutcount_matrix_raw) <- spot_gene_mutcount_matrix_raw[,1]
spot_gene_mutcount_matrix_raw[is.na(spot_gene_mutcount_matrix_raw)] <- 0
spot_gene_mutcount_matrix_raw <- as.data.frame(t(spot_gene_mutcount_matrix_raw[,-1]))
###基因突变位点数突变谱上###

###############################6-spot类别变异调控基因差异表达###############################
diff.results_cluster <- data.frame()
inter_G <- intersect(rownames(spot_gene_mutcount_matrix_raw),rownames(brain@assays$Spatial$counts))
#Exp_data_all <- as.data.frame(brain@assays$Spatial$counts)
Exp_data_all <- as.matrix(brain@assays$Spatial$counts)
for(cluster_jj in unique(brain@meta.data$mut_cluster)){
  #计算数据框Mut_data的每行和数据框Exp_data的对应行（相同行号）的相关性和显著性
  ##cluster_jj的spot名
  cluster_jj_spot <- rownames(brain@meta.data[brain@meta.data$mut_cluster %in% cluster_jj,])
  Mut_data <- data.frame(t(spot_gene_mutcount_matrix_raw[inter_G,cluster_jj_spot]))
  Exp_data <- data.frame(t(Exp_data_all[inter_G,cluster_jj_spot]))
  # 秩和检验
  diff_test <- mapply(function(x, y) {
    if(length(which(x>0))*length(which(x==0))>0){
      p_value <- wilcox.test(y[which(x>0)], y[which(x==0)])$p.value
      fc <- mean(y[which(x>0)])/mean(y[which(x==0)])
    }else{
      p_value <- 1
      fc <- NA
    }
    list(fc = fc, p_value = p_value)
  }, Mut_data, Exp_data)
  diff.results <- as.data.frame(t(diff_test))
  diff.results$mut_G <- rownames(diff.results)
  diff.results$mut_cluster <- cluster_jj
  diff.results_cluster <- rbind(diff.results_cluster,diff.results)
}
#fc和p_value的类型是list，转换成数据框
diff.results_cluster$fc <- as.character(diff.results_cluster$fc)
diff.results_cluster$p_value <- as.character(diff.results_cluster$p_value)

diff.results_cluster <- diff.results_cluster[grepl("^\\d",diff.results_cluster$fc),]
diff.results_cluster$fc <- as.numeric(unlist(diff.results_cluster$fc))
diff.results_cluster$p_value <- as.numeric(unlist(diff.results_cluster$p_value))
# 按照mut_cluster列分组，选出每组中abs(fc)*(-log(p_value))最大的前16行
df_selected_lab <- diff.results_cluster %>%
  group_by(mut_cluster) %>%
  arrange(desc(abs(fc)*(-log(p_value)))) %>%
  slice_head(n = 16) %>%
  ungroup()
# shuju

# 画分组的图
library(reshape2)
Exp_muttype_cluster <- data.frame()
for(cluster_jj in unique(df_selected_lab$mut_cluster)){
  #select数据框Mut_data的每行和数据框Exp_data的对应行（相同行号）的data
  ##cluster_jj的spot名
  cluster_jj_spot <- rownames(brain@meta.data[brain@meta.data$mut_cluster %in% cluster_jj,])
  ##cluster_jj的gene名
  cluster_jj_gene <- df_selected_lab[df_selected_lab$mut_cluster %in% cluster_jj,]$mut_G
  Mut_data <- data.frame(t(spot_gene_mutcount_matrix_raw[cluster_jj_gene,cluster_jj_spot]))
  Mut_data$spot <- rownames(Mut_data)
  Mut_data_melt <- melt(Mut_data, id.vars = "spot", variable.name="mut_G",value.name="mut_type")
  Mut_data_melt$mut_type <- ifelse(Mut_data_melt$mut_type>0,1,0)
  Exp_data <- data.frame(t(Exp_data_all[cluster_jj_gene,cluster_jj_spot]))
  Exp_data$spot <- rownames(Exp_data)
  Exp_data_melt <- melt(Exp_data, id.vars = "spot", variable.name="mut_G",value.name="G_Exp")
  Exp_muttype_cluster_bridge <- cbind(Exp_data_melt,mut_type=Mut_data_melt$mut_type)
  Exp_muttype_cluster_bridge$cluster <- cluster_jj
  Exp_muttype_cluster <- rbind(Exp_muttype_cluster,Exp_muttype_cluster_bridge)
}

library(ggpubr)
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/similarity_boxplot.png", width = 1225, height = 650)
ggplot(Exp_muttype_cluster,aes(mut_G,G_Exp))+
  geom_boxplot( aes(color = as.character(mut_type)))+
  scale_color_manual(values = c("#00AFBB","#E7B800"))+
  stat_compare_means(aes(group = mut_type),label = "p.signif")+
  facet_wrap(~ as.factor(cluster),ncol = 3,scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Mutant gene", y = "Gene expression", color = "Mutant type")
dev.off()

diff.results_cluster[,c(1,2)] <- round(diff.results_cluster[,c(1,2)], 4)
write.table(diff.results_cluster, "E:/数据库/STMut/数据处理/结果表格/similarity_table.txt", sep = "\t", quote = F, row.names = F)

