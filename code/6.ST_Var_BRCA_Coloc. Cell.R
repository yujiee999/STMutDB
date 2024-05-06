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
# ###############################数据的预处理################################
# ##数据标准化、识别高变基因,20秒
# brain<- SCTransform(brain, assay = "Spatial", verbose = FALSE) #用SCTransform对数据进行标准化, 同时检测高变基因, 输出结果储存在 SCT assay中；#我们发现不同spot的mRNA分子数差异很大，特别是当组织中的细胞密度存在差异时，因此需要对数据进行标准化。由于细胞组织存在异质性，不同组织区域的细胞密度可能不同；因此如果采用常规单细胞转录组数据的LogNormalize标准化方法可能存在问题；这里Seurat作者推荐sctransform方法。
# ##展示高变基因
# p5<- VariableFeaturePlot(brain,cols = c( "gray60", "red"))
# top10<- head(VariableFeatures(brain),10) 
# p6<- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0) 
# #p6#解释图纵轴：残差方差（residual variance），将观测值与模型预测值之间的差异的平方和除以自由度。即基因在每个细胞中的表达与表达均值的平方和除以减一后的细胞数目
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
# 
# ###########################################################################
# ##########################先跑出p10以前的代码，然后##########################
# ###########################################################################
# ##########################基因突变0|1突变谱聚类###############################
# ############################################################################
# #参考 GSE167036_Cell_subsets_definedby_mutGenes.R
# file_path <- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/Step4_VariantCalling"
# ##########################0-读取scsnv注释矩阵###############################
# anno<- read.table(paste0(file_path,"/GSM6177599.variants.annovar.Step4.3.1.hg38_multianno.csv"), header=T,sep=",",as.is=T,fill=T)
# anno<-anno[,c(1,2,7)]
# anno<- anno[nchar(anno[,3]) > 1,]# 保留第三列字符长度大于1的行
# anno<-anno[!grepl("\\|", anno[,3]), ]
# library(tidyr)
# anno2 <-anno %>% as_tibble() %>% 
#   separate_rows(Gene.refGene, sep = ";")
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
# spot_gene_mutcount_matrix <- dcast(spot_gene_mutcount, CB ~ Gene.refGene, value.var = "n")
# ##参考单细胞处理流程https://www.jianshu.com/p/728125be7c53
# ##参考单细胞处理流程代码GSE162708_infercnv.R
# spot_gene_mutcount_matrix[,1] <- paste0(spot_gene_mutcount_matrix[,1],"-1")
# #基于空转整合spot
# spot_gene_mutcount_matrix <- merge(data.frame(CB=rownames(brain@meta.data)),spot_gene_mutcount_matrix,all.x=T)
# rownames(spot_gene_mutcount_matrix)<-spot_gene_mutcount_matrix[,1]
# spot_gene_mutcount_matrix[is.na(spot_gene_mutcount_matrix)]<-0
# #spot_gene_mutcount_matrix[, -1][spot_gene_mutcount_matrix[, -1] > 0] <- 1
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
# #seurat_object <- RunTSNE(seurat_object, dims = 1:30)
# DimPlot(seurat_object , reduction = "umap", pt.size=0.5)
# ##########################4-空转展示scsnv突变聚类类别###############################
# ##取出空转metadata的行名构成数据框
# ST_rownames <- data.frame(CB=rownames(brain@meta.data))#2340
# mut_cluster <- data.frame(CB=rownames(seurat_object@meta.data),seurat_object@meta.data[,-1])
# colnames(mut_cluster)[ncol(mut_cluster)] <- "mut_cluster"
# mut_cluster <- merge(ST_rownames,mut_cluster,all.x=T)
# mut_cluster$mut_cluster <- as.character(mut_cluster$mut_cluster)
# mut_cluster[is.na(mut_cluster)] <- "WT"
# brain@meta.data <- cbind(brain@meta.data, mut_cluster[,-1])
# SpatialDimPlot(brain, group.by="mut_cluster",label = TRUE, pt.size.factor = 2,label.size = 12,alpha = c(0.8, 1),image.alpha = 1)+ 
#   theme(text = element_text(size=24))

brain <- readRDS("E:/数据库/STMut/数据处理/中间文件/ST_Mut_SeuratObject.rds")
spot_gene_mutcount_matrix <- read.table("E:/数据库/STMut/数据处理/中间文件/spot_gene_mutcount_matrix.txt", sep = "\t", header = T)
##########################5-展示scsnv突变聚类的每个类别的共现细胞突变###############################
library("gridExtra")

#SORC数据库的去卷积结果
#spot_cell <- read.table("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/GSM6177599_slice1_Deconvolution.txt",sep="\t",header = T)
#spot_cell_bridge <- t(spot_cell)

#CROST数据库的去卷积结果
spot_cell <- read.table("E:/数据库/STMut/数据处理/BRCA/paired_singlecell/download_page_VISDS000449_celltype_deco.csv", sep = ",", header = T, row.names = 1, check.names = F)
spot_cell_bridge <- spot_cell

rownames(spot_cell_bridge) <- gsub("\\.", "-", rownames(spot_cell_bridge))
spot_gene_mutcount_matrix_bridge <- spot_gene_mutcount_matrix
Cell_Colocalization_Mut_p_cluster <- data.frame()
for(mut_cluster_i in unique(brain@meta.data$mut_cluster)){
  aa <- rownames(brain@meta.data[brain@meta.data$mut_cluster==mut_cluster_i,])
  spot_cell <- spot_cell_bridge[aa,]
  spot_gene_mutcount_matrix <- spot_gene_mutcount_matrix_bridge[rownames(brain@meta.data[brain@meta.data$mut_cluster==mut_cluster_i,]),]
  ###############################8.1-计算细胞和变异基因共现性###############################
  ##产生所有共现的细胞和基因 组合
  all_combinations_df <- expand.grid(cell_A = colnames(spot_cell), Mut_G = colnames(spot_gene_mutcount_matrix[,-1]))
  start_time <- Sys.time()
  library(foreach)
  library(doParallel)
  num_cores <- 7  # 设置为你的计算机核心数
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  ##计算共现的细胞、基因组合的显著性
  Colocalization_Mut_p <- data.frame()
  Colocalization_Mut_p <- foreach(i = 1:nrow(all_combinations_df), .combine='rbind') %dopar% {## 使用 foreach 并行化循环
    # 计算 Colocalization_Mut、Fisher_p，填充 Colocalization_Mut_p 的现有代码
    A_spot <- rownames(spot_cell)[which(spot_cell[,all_combinations_df$cell_A[i]]*20>=1)]#出现细胞A的spot名
    G_spot <- rownames(spot_gene_mutcount_matrix)[which(spot_gene_mutcount_matrix[,all_combinations_df$Mut_G[i]]>0)]#出现基因G突变的spot名
    NoG_spot <- setdiff(rownames(spot_gene_mutcount_matrix),G_spot)#出现基因G非突变的spot名
    Colocalization_Mut <- length(intersect(A_spot,G_spot))#细胞A、基因G突变 共现时的spot个数
    Colocalization_NoMut <- length(setdiff(A_spot,G_spot))#细胞A出现,基因G非突变的spot个数
    NoColocalization_Mut <- length(setdiff(G_spot,A_spot))#细胞A非出现，基因G突变 的spot个数
    NoColocalization_NoMut <- length(setdiff(NoG_spot,A_spot))#细胞A非出现，基因G非突变 的spot个数
    Fisher_p <- 1-phyper(Colocalization_Mut-1,#两个数的交集减一
                       Colocalization_Mut+NoColocalization_Mut,#其中一个数
                       Colocalization_NoMut+NoColocalization_NoMut,#背景减其中一个数
                       Colocalization_Mut+Colocalization_NoMut)#另一个数
    Colocalization_Mut_p[i, 1] <- all_combinations_df$cell_A[i]
    Colocalization_Mut_p[i, 2] <- all_combinations_df$Mut_G[i]
    Colocalization_Mut_p[i, 3] <- Colocalization_Mut
    Colocalization_Mut_p[i, 4] <- Fisher_p
    Colocalization_Mut_p[i, ]
  }
  # 关闭并行化
  stopCluster(cl)
  end_time <- Sys.time()
  # 计算循环运行时间
  loop_time <- end_time - start_time
  print(loop_time)
  
  ###############################8.2-细胞、变异基因显著共现热图###############################
  Cell_Colocalization_Mut_p <- Colocalization_Mut_p
  Cell_Colocalization_Mut_p$Co <- ifelse(Cell_Colocalization_Mut_p[,4]<=0.05,1,0)
  colnames(Cell_Colocalization_Mut_p)[1:4] <- c("cell_A","Mut_G","Colocalization","Fisher_p")
  # 按照col分组求和并按照从大到小排序
  df_sum_col <- Cell_Colocalization_Mut_p %>%
    group_by(Mut_G) %>%
    summarise(sum_value = sum(Co)) %>%
    arrange(desc(sum_value))
  
  # 按照row分组求和并按照从大到小排序
  df_sum_row <- Cell_Colocalization_Mut_p %>%
    group_by(cell_A) %>%
    summarise(sum_value = sum(Co)) %>%
    arrange(desc(sum_value))
  
  Cell_Colocalization_Mut_p_1 <- Cell_Colocalization_Mut_p[Cell_Colocalization_Mut_p$Mut_G %in% df_sum_col$Mut_G[1:12],]
  Cell_Colocalization_Mut_p_1$cell_A <- factor(Cell_Colocalization_Mut_p_1$cell_A,levels=df_sum_row$cell_A)
  Cell_Colocalization_Mut_p_1$Mut_G <- factor(Cell_Colocalization_Mut_p_1$Mut_G,levels=rev(df_sum_col$Mut_G))
  Cell_Colocalization_Mut_p_2 <- Cell_Colocalization_Mut_p_1 %>% mutate(text = case_when(
    Fisher_p < 0.01 ~ "*"))
  p_plot <- ggplot(Cell_Colocalization_Mut_p_2, aes(x = cell_A, y = Mut_G, fill = Colocalization)) +
    geom_tile()+
    #geom_tile(aes(width = Colocalization/100, height = Colocalization/100)) +
    scale_fill_gradient(low = "white", high = "#0F4400") +
    geom_text(aes(label=text),col ="black",size = 5) +
    labs(title = paste0("mut_cluster: ",mut_cluster_i), x = "Spot.Cell", y = "Spot.Mut_Gene",fill = "No.Spot") +
    theme_set(theme_bw()) +
    theme(title=element_text(size=16),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14), 
          axis.text.y = element_text(size=12,color = 'black'),
          axis.text.x = element_text(size=12,color = 'black',angle = 90, hjust = 1))
  
  co_local <- "E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/"
  png(paste0(co_local, "co-localization",mut_cluster_i , ".png"), width = 1225, height = 650)
  print(p_plot)
  dev.off()
  
  eval(parse(text=paste("p",mut_cluster_i,"<-p_plot",sep="")))
  Cell_Colocalization_Mut_p_2$cluster <- mut_cluster_i
  Cell_Colocalization_Mut_p_cluster <- rbind(Cell_Colocalization_Mut_p_cluster,Cell_Colocalization_Mut_p_2)
}
eval(parse(text=paste("plist<-list(",paste("p",0:mut_cluster_i,sep="",collapse = ","),")",sep="")))
library("grid")
p1 <- grid.arrange(grobs = plist, ncol = 3)
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/co-localization.png", width = 1225, height = 650)
print(p1)
dev.off()

coloc <- Cell_Colocalization_Mut_p_cluster[,-c(5,6)]
coloc$Fisher_p <- round(coloc[,4], 4)
write.table(coloc, "E:/数据库/STMut/数据处理/结果表格/coloc.txt", sep = "\t", quote = F, row.names = F)


