library(Seurat)
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r) #
library(tidyverse)
################################################提取细胞表达谱################################################
# counts = Read10X("E:/数据库/STMut/数据处理/BRCA/paired_singlecell/BRCA/BRCA_GSE167036.expression/GSM5091212/")
# #dim(counts) #36601*7770
# #seurat_object <- CreateSeuratObject(counts = counts$`Gene Expression`)
# seurat_object <- CreateSeuratObject(counts = counts)
# cell_exp_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
# cell_exp_G <- seurat_object@assays[["RNA"]]@features %>% rownames()
# cell_exp_spot <- seurat_object@assays[["RNA"]]@cells %>% rownames()
# colnames(cell_exp_mat) <- cell_exp_spot
# rownames(cell_exp_mat) <- cell_exp_G
#cell_exp_mat<-cell_exp_mat[,1:500]###选择部分的细胞作为参考
#*************************************************************************************************************
#*只选取免疫细胞作为参考
library(data.table)
sc_data <- fread("E:/数据库/STMut/数据处理/BRCA/paired_singlecell/SC_data.txt")
sc_celltype <- read.table("E:/数据库/STMut/数据处理/BRCA/paired_singlecell/SC_CellType.txt", header = T, sep = "\t", row.names = 1)
sc_celltype_immune <- sc_celltype[-which(sc_celltype$cell_type_final == "Tumor"),]$cell_name
if(length(sc_celltype_immune)>=1000){
  sc_celltype_immune_1000 <- sc_celltype_immune[1:1000]
}else{
  sc_celltype_immune_1000 <- sc_celltype_immune
}

rownames(sc_data) <- sc_data$V1
seurat_object <- CreateSeuratObject(counts = sc_data)
seurat_object <- subset(seurat_object, cells = sc_celltype_immune_1000)
cell_exp_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
cell_exp_G <- seurat_object@assays[["RNA"]]@features %>% rownames()
cell_exp_spot <- seurat_object@assays[["RNA"]]@cells %>% rownames()
colnames(cell_exp_mat) <- cell_exp_spot
rownames(cell_exp_mat) <- cell_exp_G
#*************************************************************************************************************
################################################提取spot表达谱################################################
data_dir <- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/"
list.files(data_dir) #查看目录中的文件；
file_name <- "filtered_feature_bc_matrix.h5"
brain <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice = "anterior1")
brain@project.name <- "anterior1"
Idents(brain) <- "anterior1"
brain$orig.ident <- "anterior1"
for (i in colnames((brain@images$anterior1@coordinates))) {
  brain@images$anterior1@coordinates[[i]] <- as.integer(brain@images$anterior1@coordinates[[i]])
}
spot_exp_mat <- brain@assays[["Spatial"]]@layers[["counts"]]
spot_exp_G <- brain@assays[["Spatial"]]@features %>% rownames()
spot_exp_spot <- brain@assays[["Spatial"]]@cells %>% rownames()
colnames(spot_exp_mat) <- spot_exp_spot
rownames(spot_exp_mat) <- spot_exp_G
################################################整合细胞和spot的表达谱################################################
inter_G <- intersect(rownames(spot_exp_mat),rownames(cell_exp_mat))
exp_mat_1 <- cbind(spot_exp_mat[inter_G,],cell_exp_mat[inter_G,])
cell_type <- data.frame(Type=c(rep("spot",ncol(spot_exp_mat)),rep("cell",ncol(cell_exp_mat))),row.names = colnames(exp_mat_1))
################################################infercnv################################################
#
# inter_G_1 <- as.data.frame(inter_G)
# gencode_anno <- read.table("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/gencode.v43.annotation.txt", sep = "\t", header = T)
# geneOrderingFile <- gencode_anno[gencode_anno$gene_name %in% inter_G_1$inter_G, c(5,1,2,3)]
# write.table(geneOrderingFile, "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/geneOrderingFile.txt", sep = "\t", row.names = F, quote = F, col.names = F)
#
library(infercnv)
#第一步，根据上述的三个文件创建inferCNV对象
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exp_mat_1, # 可以直接提供矩阵对象
                                    annotations_file=cell_type,
                                    delim="\t",
                                    gene_order_file="E:/数据库/STMut/suminghai/slice1_gene_ordering_file.txt",
                                    ref_group_names="cell",
                                    min_max_counts_per_cell=c(0,+Inf))
#这一步的一个关键参数是ref_group_name, 用于设置参考组。假如你并不知道哪个组是正常，哪个组不正常，
#那么设置为ref_group_name=NULL, 那么inferCNV会以全局平均值作为基线，这适用于有足够细胞存在差异的情况。
#此外，这里的raw_count_matrix是排除了低质量细胞的count矩阵。

#第二步，运行标准的inferCNV流程。
setwd("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/inferCNV")
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # 输出文件夹
                             cluster_by_groups=F,   # 聚类
                             denoise=T, #去噪
                             HMM=F,# 是否基于HMM预测CNV
                             #no_plot = TRUE,
                             no_prelim_plot = FALSE,
                             num_threads = 12)

###########
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
#infercnv_obj = readRDS("C:/Users/DELL/Documents/output_dir/run.final.infercnv_obj")
infercnv_obj = readRDS("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/inferCNV/output_dir_1/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
gn <- rownames(expr)
geneFile <- read.table("E:/数据库/STMut/suminghai/slice1_gene_ordering_file.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile) = geneFile$V1
sub_geneFile <- geneFile[intersect(gn,geneFile$V1),]
expr = expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]
#聚类，7类，提取结果
set.seed(20240403)
kmeans.result <- kmeans(t(expr), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB = rownames(kmeans_df)
#kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)

#定义热图的注释，及配色
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #类别数
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

#下面是绘图
pdf("try1_1.pdf",width = 25,height = 20)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/cnv.png", width = 700, height = 400)
print(ht)
dev.off()
###############################读取空转数据################################
data_dir <- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599"
list.files(data_dir) #查看目录中的文件；
file_name <- "filtered_feature_bc_matrix.h5"
brain <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice = "anterior1")

brain@project.name <- "anterior1"
Idents(brain) <- "anterior1"
brain$orig.ident <- "anterior1"
for (i in colnames((brain@images$anterior1@coordinates))) {
  brain@images$anterior1@coordinates[[i]] <- as.integer(brain@images$anterior1@coordinates[[i]])
}
brain <- subset(brain, subset = nCount_Spatial > 0)#spot的mrna分子数
###############################数据的预处理################################
##数据标准化、识别高变基因,20秒
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) #用SCTransform对数据进行标准化, 同时检测高变基因, 输出结果储存在 SCT assay中；#我们发现不同spot的mRNA分子数差异很大，特别是当组织中的细胞密度存在差异时，因此需要对数据进行标准化。由于细胞组织存在异质性，不同组织区域的细胞密度可能不同；因此如果采用常规单细胞转录组数据的LogNormalize标准化方法可能存在问题；这里Seurat作者推荐sctransform方法。
##展示高变基因
p5 <- VariableFeaturePlot(brain,cols = c( "gray60", "red"))
top10 <- head(VariableFeatures(brain),10) 
p6 <- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0) 
p6#解释图纵轴：残差方差（residual variance），将观测值与模型预测值之间的差异的平方和除以自由度。即基因在每个细胞中的表达与表达均值的平方和除以减一后的细胞数目
###############################降维、无监督聚类################################
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
brain@meta.data$CNV_class <- as.character(kmeans_df[rownames(brain@meta.data),]$kmeans_class)
brain@meta.data$CNV_class[is.na(brain@meta.data$CNV_class)] <- "other"
#展示CNV定义的肿瘤spot
cnv_spot <- SpatialDimPlot(brain, group.by="CNV_class",label = TRUE, pt.size.factor = 3,label.size = 12,alpha = c(0.8, 1),image.alpha = 0)+ 
  theme(text = element_text(size=24))

png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/cnv_spot.png", width = 670, height = 360)
print(cnv_spot)
dev.off()
