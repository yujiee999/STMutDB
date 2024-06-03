library(Seurat)
# 中间文件
input_inter_path <- "E:/数据库/STMut/数据处理/结果/inter"

# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_inter_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    tmp_CNVpath <- paste0(input_inter_path, "/", dataset_name, "/", sample_name, "/output_dir")
    setwd(tmp_CNVpath)
    ##获得CNV_cluster的spot_cluster关系
    infercnv_obj = readRDS(paste0(tmp_CNVpath, "/run.final.infercnv_obj"))
    expr <- infercnv_obj@expr.data
    gn <- rownames(expr)
    geneFile <- read.table("E:/数据库/STMut/suminghai/slice1_gene_ordering_file.txt",header = F,sep = "\t",stringsAsFactors = F)
    rownames(geneFile)=geneFile$V1
    sub_geneFile <- geneFile[intersect(gn,geneFile$V1),]
    expr=expr[intersect(gn,geneFile$V1),]
    head(sub_geneFile,4)
    expr[1:4,1:4]
    #聚类，7类，提取结果
    set.seed(20240603)
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
    
    brain <- readRDS(paste0(input_inter_path, "/", dataset_name, "/", sample_name, "/ST_seuratObject.rds"))
    brain@meta.data$CNV_class <- as.character(kmeans_df[rownames(brain@meta.data),]$kmeans_class)
    brain@meta.data$CNV_class[is.na(brain@meta.data$CNV_class)] <- "other"
    
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
    png(paste0("E:/数据库/STMut/数据处理/结果/result/", dataset_name, "/", sample_name, "/cnv_slice.png"), width = 674, height = 674)
    SpatialDimPlot(brain,group.by="CNV_class",label = TRUE, pt.size.factor = 3,label.size = 6,alpha = c(0.8, 1),image.alpha = 1)+scale_fill_manual(values=my_colors)+
      theme(text = element_text(size=24))
    dev.off()
    
  }
}
