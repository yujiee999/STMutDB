#载入所需的R包；
library(Seurat)#
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r) #

input_path <- "E:/数据库/STMut/数据处理/结果"
output_path <- "E:/数据库/STMut/数据处理/结果/result"

#在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  for(gsm_dir in gsm_dirs) {
    # 在每个 GSM 文件夹下找到 Step4_VariantCalling 文件夹
    step4_dir <- list.dirs(path = gsm_dir, full.names = TRUE)
    step4_dir <- step4_dir[grepl("Step4_VariantCalling", basename(step4_dir))]
    # 1.得到每个spot的突变负荷
    file_list <- list.files(path = step4_dir, pattern = "genotype.tsv$", full.names = TRUE)
    spot_mutALL <- data.frame()
    for (file in file_list) {
      df2 <- read.table(file, header=T,sep="\t",as.is=T,comment.char = "")
      df2 <- df2[which(df2$ALT_expected!="."),]
      df2 <- df2[which(df2$Base_observed!=df2$REF),]
      df3 <- data.frame(table(unique(df2[,c(1,2,8)])$CB))
      colnames(df3) <- c("CB","MutPerSpot")
      spot_mutALL <- rbind(spot_mutALL,df3)
    }
    head(spot_mutALL)#每个spot的突变位点数
    
    # 2.得到每个spot簇注释，并画出突变负荷图
    # 这里直接读取中间文件，得到空转数据+突变类别的seurat对象
    inter_path <- gsub(input_path, "E:/数据库/STMut/数据处理/结果/inter", gsm_dir)
    brain <- readRDS(paste0(inter_path,"/" ,"ST_Mut_SeuratObject.rds"))
    #画图
    spot_cluster <- data.frame(CB=sub("-1$", "",rownames(brain@meta.data)),cluster=brain@meta.data$mut_cluster)
    spot_cluster_mut <- merge(spot_mutALL,spot_cluster)
    p3 <- ggplot(spot_cluster_mut, aes(x = cluster, y = MutPerSpot))+
      geom_jitter(width = 0.3,color = "black", height = 0, size = 3,shape = 21,fill = "gray",stroke = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", width = 0.7, color = "red") +
      labs(title = "Plot of Mut by Cluster special", x = "Cluster", y = "Muts. perSpot")+
      theme_classic()
    #输出结果的路径
    result_path <- gsub(input_path, output_path, gsm_dir)
    
    png(paste0(result_path, "/mutation_load.png"), width = 674, height = 674)
    print(p3)
    dev.off()
    
    # 3.画出突变负荷决定spot簇的核心区域可能性图
    cluster_medianMut <- spot_cluster_mut %>% 
      group_by(cluster) %>%
      summarise(medianMut=median(MutPerSpot)) %>%
      as.data.frame()
    
    if(length(unique(cluster_medianMut$medianMut)) == 1){
      cluster_medianMut <- spot_cluster_mut %>% 
        group_by(cluster) %>%
        summarise(medianMut=mean(MutPerSpot)) %>%
        as.data.frame()
    }
    # 创建一个颜色渐变函数
    my_ramp <- colorRamp(c("grey", "red"))
    # 定义一个向量
    my_values <- cluster_medianMut$medianMut
    # 根据向量中的值生成对应的颜色
    get_color <- function(value) {
      color_index <- (value - min(my_values)) / (max(my_values) - min(my_values))
      return(rgb(my_ramp(color_index), maxColorValue = 255))
    }
    # 获取向量中每个值对应的颜色
    my_colors <- sapply(my_values, get_color)
    p4 <- SpatialDimPlot(brain,group.by="mut_cluster",label = TRUE, pt.size.factor = 3,label.size = 6,alpha = c(0.8, 1),image.alpha = 1)+scale_fill_manual(values=my_colors)+
      theme(text = element_text(size=24))
    
    png(paste0(result_path, "/mutation_load_slice.png"), width = 674, height = 674)
    print(p4)
    dev.off()
    
  }
}



