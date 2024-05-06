setwd("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/inferCNV/output_dir")
##获得CNV_cluster的spot_cluster关系
infercnv_obj = readRDS("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/inferCNV/output_dir/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
gn <- rownames(expr)
geneFile <- read.table("E:/数据库/STMut/suminghai/slice1_gene_ordering_file.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <- geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]
#聚类，7类，提取结果
set.seed(20240313)
kmeans.result <- kmeans(t(expr), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
#kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
#rownames(kmeans_df_s) <- substr(rownames(kmeans_df_s),1,16)
rownames(kmeans_df_s) <- gsub("-1", "", rownames(kmeans_df_s))
head(kmeans_df_s)
library(dplyr)
##参考阿伟数据库
##对于空转CNV谱，对每行的CNV进行Z-score标准化并平方，计算每个spot的CNV_score
cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    #resacle(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

data <- read.table("infercnv.observations.txt", header = T)
cnv_score <- cnvScore(data)
rownames(cnv_score)<-substr(rownames(cnv_score),1,16)
##对于参考单细胞CNV谱，对每行的CNV进行Z-score标准化并平方，计算每个细胞的CNV_score
ref_data <- read.table("infercnv.references.txt", header=T)
ref_cnv_score <- cnvScore(ref_data)
##计算每个CNV_cluster中spot拷贝数和参考细胞拷贝数差异程度FC*(-log(p))
cluster_CNVsim <- data.frame(cluster=character(),sim=numeric())
for(i in unique(kmeans_df_s$kmeans_class)){
  cluster_cnv <- cnv_score[intersect(rownames(cnv_score),rownames(kmeans_df_s)[which(kmeans_df_s$kmeans_class ==i)]),1]
  if(length(cluster_cnv)>0){
    result <- wilcox.test(cluster_cnv, ref_cnv_score$`colSums(data * data)`, alternative = "two.sided")
    FC <- mean(cluster_cnv)/mean(ref_cnv_score$`colSums(data * data)`)
    rowindex <- nrow(cluster_CNVsim)+1
    cluster_CNVsim[i,1] <- i
    cluster_CNVsim[i,2] <- -log(result$p.value)*FC
  }
}
library(Seurat)
###############################读取空转数据################################
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
brain <- subset(brain, subset = nCount_Spatial > 0)
###############################数据的预处理################################
##数据标准化、识别高变基因,20秒
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) #用SCTransform对数据进行标准化, 同时检测高变基因, 输出结果储存在 SCT assay中；#我们发现不同spot的mRNA分子数差异很大，特别是当组织中的细胞密度存在差异时，因此需要对数据进行标准化。由于细胞组织存在异质性，不同组织区域的细胞密度可能不同；因此如果采用常规单细胞转录组数据的LogNormalize标准化方法可能存在问题；这里Seurat作者推荐sctransform方法。
##展示高变基因
p5 <- VariableFeaturePlot(brain,cols = c( "gray60", "red"))
top10 <- head(VariableFeatures(brain),10) 
p6 <- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0) 
#p6#解释图纵轴：残差方差（residual variance），将观测值与模型预测值之间的差异的平方和除以自由度。即基因在每个细胞中的表达与表达均值的平方和除以减一后的细胞数目
###############################降维、无监督聚类################################
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
brain@meta.data$CNV_class <- as.character(kmeans_df[rownames(brain@meta.data),]$kmeans_class)
brain@meta.data$CNV_class[is.na(brain@meta.data$CNV_class)] <- "other"
#**************************************************************************
# other <- brain@meta.data$CNV_class == "other"
# tmp_other <- brain@meta.data
#brain@meta.data <- brain@meta.data[-which(brain@meta.data$CNV_class == "other"),]
#**************************************************************************
# 创建一个颜色渐变函数
my_ramp <- colorRamp(c("grey", "red"))
# 定义一个向量
my_values <- cluster_CNVsim$sim
# 根据向量中的值生成对应的颜色
get_color <- function(value) {
  color_index <- (value - min(my_values)) / (max(my_values) - min(my_values))
  return(rgb(my_ramp(color_index), maxColorValue = 255))
}
# 获取向量中每个值对应的颜色
my_colors <- sapply(my_values, get_color)
##并画出CNV决定spot簇的核心区域可能性图
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/cnv_slice.png", width = 670, height = 360)
SpatialDimPlot(brain,group.by="CNV_class",label = TRUE, pt.size.factor = 3,label.size = 6,alpha = c(0.8, 1),image.alpha = 1)+scale_fill_manual(values=my_colors)+
  theme(text = element_text(size=24))
dev.off()


