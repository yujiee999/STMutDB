#********************************************************
#*获取切片图像坐标和key
key(object = brain@images$anterior1)
head(GetTissueCoordinates(brain))
head(Idents(brain))

library(Seurat)#
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r)
library(reshape2)

input_path <- "E:/数据库/STMut/数据处理/结果/inter"
cnv_path <- "G:/database/STMut/data_process/inter/"
# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

# 定义一个函数，将多列值按逗号分隔开
combine_values <- function(x) {
  paste(x, collapse = ",")
}
combine_values_1 <- function(x) {
  paste(x, collapse = ";")
}

for (i in 1:40) {
  gse_dir <- gse_dirs[i]
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(kk in 1:length(gsm_dirs)) {
    
    gsm_dir <- gsm_dirs[kk]
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    brain <- readRDS(paste0(gsm_dir, "/ST_Mut_SeuratObject.rds"))
    #1.提取每个spot的空间坐标+突变的cluster
    Slice_Coor_cluster <- brain@meta.data[,c(2,3,4,5,9,10,11)]
    spot_corrdinates <- GetTissueCoordinates(brain)
    spot_corrdinates$CB <- rownames(spot_corrdinates)
    
    Slice_Coor_Cluster_integrate <- cbind(Slice_Coor_cluster, spot_corrdinates)
    
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
      # mrna_cluster = paste(result_1$mRNA_cluster, collapse = ";"),
      # cnv_cluster = paste(result_1$CNV_cluster, collapse = ";")
    )
    result_3 <- as.data.frame(t(result_2))
    colnames(result_3)[1] <- "Column_value"
    result_3$Column_name <- rownames(result_3)
    
    #2.突变基因的突变位点数
    #突变谱
    #将宽数据转换成长数据
    library(tidyverse)
    spot_gene_mutcount_matrix <- read.table(paste0(gsm_dir, "/spot_gene_mutcount_matrix.txt"), sep = "\t", row.names = 1, check.names = FALSE)
    long_data <- spot_gene_mutcount_matrix %>% 
      gather(key = "Gene", value = "MutCount", -CB)
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
    long_data_4 <- aggregate(. ~ Gene, data = long_data_3, FUN = combine_values)
    long_data_5 <- long_data_4[,c(1,4)]
    colnames(long_data_5) <- c("Column_name", "Column_value")
    
    result_4 <- rbind(result_3, long_data_5)
    
    #3.加上位点的信息，如果这个位点在spot里出现突变，值为1，否则就为0
    curr_path <- gsub("/inter", "", gsm_dir)
    step4_dir <- list.dirs(path = curr_path, full.names = TRUE)
    step4_dir <- step4_dir[grepl("Step4_VariantCalling", basename(step4_dir))]
    
    file_list <- list.files(path = step4_dir, pattern = "genotype.tsv$", full.names = TRUE)
    spot_mutALL <- data.frame()
    for (file in file_list) {
      df2 <- read.table(file, header=T, sep="\t", as.is=T, comment.char = "")
      df2 <- df2[which(df2$ALT_expected!="."),]
      df2 <- df2[which(df2$Base_observed!=df2$REF),]
      #拼接位点
      # 对数据框的每一行进行操作，拼接字符串，并保存到新列 result
      df2$variants <- apply(df2, 1, function(row) {
        aa <- paste(row[1], "_", row[2], sep = "")
        bb <- paste(row[4], "_", row[5], sep = "")
        paste(aa, bb, sep = "_")
      })
      df2$variants <- gsub(" ", "", df2$variants)
      df3 <- df2[,c("CB", "variants")]
      spot_mutALL <- rbind(spot_mutALL,df3)
    }
    spot_mutALL$CB <- paste0(spot_mutALL$CB, "-1")
    spot_mutALL$mutOrNot <- "1"
    
    combinations <- expand.grid(rownames(brain_ST@meta.data), unique(spot_mutALL$variants))#750*3620
    names(combinations) <- c("CB", "variants")
    spot_mutALL_1 <- unique(merge(combinations, spot_mutALL, by = c("CB", "variants"), all = T))
    # 将第三列中的NA值替换为0
    spot_mutALL_1$mutOrNot[is.na(spot_mutALL_1$mutOrNot)] <- 0
    #按照位点去分组，组内内部将细胞排序,顺序要求和Slice_Coor_Cluster_integrate$CB的顺序一致
    spot_mutALL_1_sorted <- spot_mutALL_1 %>% 
                        group_by(variants) %>%
                        arrange(match(CB, Slice_Coor_Cluster_integrate_2$CB))
    spot_mutALL_2_sorted <- aggregate(. ~ variants, data = spot_mutALL_1_sorted, FUN = combine_values)
    spot_mutALL_3_sorted <- spot_mutALL_2_sorted[,-2]
    names(spot_mutALL_3_sorted) <- c("Column_name", "Column_value")
    
    result_5 <- rbind(result_4, spot_mutALL_3_sorted)
    
    result_5$sample <- sample_name
    
    write.table(result_5, paste0(gsub("E:/数据库/STMut/数据处理/结果/inter", "E:/数据库/STMut/数据处理/结果/result_table", gsm_dir), "/Slice_Coor_MutCluster.txt"), row.names = F, quote = F, sep = "\t")
    
    #4.mRNA的cluster---------------------------------------------------------------------------------------------------------
    brain_ST <- readRDS(paste0(gsm_dir, "/ST_seuratObject.rds"))
    brain_ST_cluster <- data.frame(CB = rownames(brain_ST@meta.data),mRNA_cluster = brain_ST@meta.data$seurat_clusters)
    #合并mrna_cluster
    Slice_Coor_Cluster_integrate_1 <- merge(Slice_Coor_Cluster_integrate, brain_ST_cluster, by = "CB", all = T)
    Slice_Coor_Cluster_integrate_1 <- Slice_Coor_Cluster_integrate_1[,c("CB", "imagerow", "imagecol", "mRNA_cluster")]
    #保证CB的顺序一致，所以先排序
    Slice_Coor_Cluster_integrate_1 <- Slice_Coor_Cluster_integrate_1[order(Slice_Coor_Cluster_integrate_1$mRNA_cluster, Slice_Coor_Cluster_integrate_1$CB), ]
    # 按照group分组，将value1、value2、value3按逗号分隔开，不同分组用分号隔开
    # 将第一列group的值也按逗号分隔开
    Slice_Coor_Cluster_integrate_1$mRNA_cluster_1 <- Slice_Coor_Cluster_integrate_1$mRNA_cluster
    Slice_Coor_Cluster_integrate_1$mRNA_cluster_1 <- as.character(Slice_Coor_Cluster_integrate_1$mRNA_cluster_1)
    
    result_1 <- aggregate(. ~ mRNA_cluster_1, data = Slice_Coor_Cluster_integrate_1, FUN = combine_values)
    result_2 <- data.frame(
      imagerow = paste(result_1$imagerow, collapse = ";"),
      imagecol = paste(result_1$imagecol, collapse = ";"),
      CB = paste(result_1$CB, collapse = ";"),
      mrna_cluster = paste(result_1$mRNA_cluster, collapse = ";")
    )
    result_3 <- as.data.frame(t(result_2))
    colnames(result_3)[1] <- "Column_value"
    result_3$Column_name <- rownames(result_3)
    
    result_3$sample <- sample_name
    write.table(result_3, paste0(gsub("E:/数据库/STMut/数据处理/结果/inter", "E:/数据库/STMut/数据处理/结果/result_table", gsm_dir), "/Slice_Coor_mRNACluster.txt"), row.names = F, quote = F, sep = "\t")
    
    #5.CNV的cluster---------------------------------------------------------------------------------------------------------
    tmp_CNVpath <- paste0(cnv_path, dataset_name, "/", sample_name, "/output_dir")
    infercnv_obj = readRDS(paste0(tmp_CNVpath, "/run.final.infercnv_obj"))
    expr <- infercnv_obj@expr.data
    kmeans.result <- kmeans(t(expr), 7)
    kmeans_df <- data.frame(CNV_cluster=kmeans.result$cluster)
    kmeans_df$CB=rownames(kmeans_df)
    
    Slice_Coor_Cluster_integrate_2 <- Slice_Coor_Cluster_integrate_1[,1:3]
    #合并CNV_cluster
    Slice_Coor_Cluster_integrate_3 <- merge(Slice_Coor_Cluster_integrate_2, kmeans_df, by= "CB", all.x = T)

    #排序
    Slice_Coor_Cluster_integrate_3 <- Slice_Coor_Cluster_integrate_3[order(Slice_Coor_Cluster_integrate_3$CNV_cluster, Slice_Coor_Cluster_integrate_3$CB), ]
    
    Slice_Coor_Cluster_integrate_3$CNV_cluster_1 <- Slice_Coor_Cluster_integrate_3$CNV_cluster
    Slice_Coor_Cluster_integrate_3$CNV_cluster_1 <- as.character(Slice_Coor_Cluster_integrate_3$CNV_cluster_1)
    
    result_1 <- aggregate(. ~ CNV_cluster_1, data = Slice_Coor_Cluster_integrate_3, FUN = combine_values)
    result_2 <- data.frame(
      imagerow = paste(result_1$imagerow, collapse = ";"),
      imagecol = paste(result_1$imagecol, collapse = ";"),
      CB = paste(result_1$CB, collapse = ";"),
      cnv_cluster = paste(result_1$CNV_cluster, collapse = ";")
    )
    result_3 <- as.data.frame(t(result_2))
    colnames(result_3)[1] <- "Column_value"
    result_3$Column_name <- rownames(result_3)
    
    result_3$sample <- sample_name
    write.table(result_3, paste0(gsub("E:/数据库/STMut/数据处理/结果/inter", "E:/数据库/STMut/数据处理/结果/result_table", gsm_dir), "/Slice_Coor_CNVCluster.txt"), row.names = F, quote = F, sep = "\t")
    
    #6.gene expression---------------------------------------------------------------------------------------------------------
    brain_ST@assays$SCT@counts[1:3,1:3]
    exp_matrix <- t(as.data.frame(brain_ST@assays$SCT@counts))
    melted_data <- melt(exp_matrix, varnames = c("CB", "ExpGene"), value.name = "Expression")
    #取前小于等于3000个高变基因和突变基因的并集
    expression_gene <- union(brain_ST@assays[["SCT"]]@var.features, colnames(spot_gene_mutcount_matrix)[-1])
    melted_data_filter <- melted_data[which(melted_data$ExpGene %in% expression_gene), ]
    
    Slice_Coor_Cluster_integrate_1 <- Slice_Coor_Cluster_integrate[,c("CB", "imagerow", "imagecol", "mut_cluster")]
    Slice_Coor_Cluster_integrate_2 <- merge(Slice_Coor_Cluster_integrate_1, melted_data_filter, by = "CB", all.y = T)
    
    Slice_Coor_Cluster_integrate_2$mut_cluster_1 <- Slice_Coor_Cluster_integrate_2$mut_cluster
    Slice_Coor_Cluster_integrate_2 <- Slice_Coor_Cluster_integrate_2[order(Slice_Coor_Cluster_integrate_2$mut_cluster_1, Slice_Coor_Cluster_integrate_2$CB), ]
    #Slice_Coor_Cluster_integrate_2$mut_cluster_1 <- as.character(Slice_Coor_Cluster_integrate_2$mut_cluster_1)
    Slice_Coor_Cluster_integrate_3 <- Slice_Coor_Cluster_integrate_2 %>%
                                      group_by(ExpGene) %>%
                                      summarize(CB = paste(CB, collapse = ","),
                                      imagerow = paste(imagerow, collapse = ","),
                                      imagecol = paste(imagecol, collapse = ","),
                                      Expression = paste(Expression, collapse = ","))
    
    Slice_Coor_Cluster_integrate_3$sample <- sample_name
    Slice_Coor_Cluster_integrate_3$ExpGeneFix <- gsub("\\.", "-", Slice_Coor_Cluster_integrate_3$ExpGene)
    write.table(Slice_Coor_Cluster_integrate_3, paste0(gsub("E:/数据库/STMut/数据处理/结果/inter", "E:/数据库/STMut/数据处理/结果/result_table", gsm_dir), "/Slice_Coor_GeneExp.txt"), row.names = F, quote = F, sep = "\t")
    
  }
}

#合并所有样本的数据
Slice_Coor_Cluster_Mut <- NULL
Slice_Coor_mRNACluster <- NULL
Slice_Coor_CNVCluster <- NULL
Slice_Coor_GeneExp <- NULL

result_path <- "E:/数据库/STMut/数据处理/结果/result_table"
# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = result_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (i in 1:40) {
  gse_dir <- gse_dirs[i]
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    slice_coor <- read.table(paste0(gsm_dir, "/Slice_Coor_MutCluster.txt"), sep = "\t", header = T, check.names = F)
    Slice_Coor_Cluster_Mut <- rbind(Slice_Coor_Cluster_Mut, slice_coor)
    
    Slice_Coor_mRNA <- read.table(paste0(gsm_dir, "/Slice_Coor_MutCluster.txt"), sep = "\t", header = T, check.names = F)
    Slice_Coor_mRNACluster <- rbind(Slice_Coor_mRNACluster, Slice_Coor_mRNA)
    
    Slice_Coor_CNV <- read.table(paste0(gsm_dir, "/Slice_Coor_CNVCluster.txt"), sep = "\t", header = T, check.names = F)
    Slice_Coor_CNVCluster <- rbind(Slice_Coor_CNVCluster, Slice_Coor_CNV)
    
    Slice_Coor_exp <- read.table(paste0(gsm_dir, "/Slice_Coor_GeneExp.txt"), sep = "\t", header = T, check.names = F)
    Slice_Coor_GeneExp <- rbind(Slice_Coor_GeneExp, Slice_Coor_exp)
    
  }
  
}
write.table(Slice_Coor_Cluster_Mut, "E:/数据库/STMut/数据处理/结果/result_table/hebing/Slice_Coor_MutCluster.txt", row.names = F, quote = F, sep = "\t")
write.table(Slice_Coor_mRNACluster, "E:/数据库/STMut/数据处理/结果/result_table/hebing/Slice_Coor_mRNACluster.txt", row.names = F, quote = F, sep = "\t")
write.table(Slice_Coor_CNVCluster, "E:/数据库/STMut/数据处理/结果/result_table/hebing/Slice_Coor_CNVCluster.txt", row.names = F, quote = F, sep = "\t")
write.table(Slice_Coor_GeneExp, "E:/数据库/STMut/数据处理/结果/result_table/hebing/Slice_Coor_GeneExp.txt", row.names = F, quote = F, sep = "\t")

# Slice_Coor_Cluster_Mut <- read.table("E:/数据库/STMut/数据处理/结果/result_table/hebing/Slice_Coor_Cluster_Mut.txt", sep = "\t", header = T)
# Slice_Coor_Cluster_Mut$Column_name <- gsub("\\.", "-", Slice_Coor_Cluster_Mut$Column_name)
# write.table(Slice_Coor_Cluster_Mut, "E:/数据库/STMut/数据处理/结果/result_table/hebing/Slice_Coor_Cluster_Mut.txt", row.names = F, sep = "\t", quote = F)









