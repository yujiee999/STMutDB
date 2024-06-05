######################################构建spot突变簇类别somatic突变
##参考https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247494325&idx=1&sn=94891e349ceb9b459e91e86db144b4ea&chksm=9b4baa0eac3c2318bfb4cfa66875078db89baf0eaf4543ed785ed4892532cbdd0145fb2ba1fd#rd
rm(list = ls())
library(Seurat)
input_path <- "E:/数据库/STMut/数据处理/结果"

# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    
      input_path <- "E:/数据库/STMut/数据处理/结果"
      sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    
      file_path <- paste0(input_path, "/", dataset_name, "/", sample_name, "/Step4_VariantCalling")
      csv_files <- list.files(path = file_path, pattern = "\\.csv$")
      anno <- read.table(paste0(file_path,"/", csv_files[1]), header=T,sep=",",as.is=T,fill=T)
      anno<-anno[,c(1,2,7)]
      anno<- anno[nchar(anno[,3]) > 1,]# 保留第三列字符长度大于1的行
      anno<-anno[!grepl("\\|", anno[,3]), ]
      library(tidyr)
      anno2 <-anno %>% as_tibble() %>% 
        separate_rows(Gene.refGene, sep = ";")
      
      file_list <- list.files(path = file_path, pattern = "genotype.tsv$", full.names = TRUE)
      gene_mut_spot<-data.frame()
      ##1-读取scsnv突变矩阵
      for (file in file_list) {
        df2<- read.table(file, header=T,sep="\t",as.is=T,comment.char = "")
        df2 <- df2[which(df2$ALT_expected!="."),]
        df2 <- df2[which(df2$Base_observed!=df2$REF),]
        colnames(df2)[1]<-"Chr"
        df2<-merge(df2,anno2)
        gene_mut_spot<-rbind(gene_mut_spot,df2)
      }
      gene_mut_spot$CB <- paste0(gene_mut_spot$CB,"-1")
      
      ##2-识别的spot突变簇类别
      data_dir <- readRDS(paste0("E:/数据库/STMut/数据处理/结果/inter/", dataset_name, "/", sample_name, "/ST_Mut_SeuratObject.rds"))
      spot_cluster <- data_dir@meta.data
      spot_cluster$mut_cluster <- paste0("C_",spot_cluster$mut_cluster)
      spot_cluster$CB <- rownames(spot_cluster)
      Maf_data <- merge(gene_mut_spot,spot_cluster)
      ######################################提取突变上下文已经计算96突变形式的比例
      library(SomaticSignatures)
      library(BSgenome.Hsapiens.UCSC.hg38)
      Maf_data <- Maf_data[which(Maf_data$REF !=Maf_data$ALT_expected),]
      sca_vr = VRanges(
        seqnames =  Maf_data$Chr ,
        ranges = IRanges(start = Maf_data$Start,end = Maf_data$Start+1),
        ref = Maf_data$REF,
        alt = Maf_data$ALT_expected,
        sampleNames = as.character(Maf_data$mut_cluster),
        study = as.character(Maf_data$mut_cluster))
      #从基因组根据坐标获取碱基上下文
      sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg38)
      # 对每个cluster，计算 96 突变可能性的 比例分布情况
      escc_sca_mm = motifMatrix(sca_motifs, group = "study", normalize = F)
      
      ######################################空间分布簇96个三核苷酸突变图
      #  install the latest version from GitHub
      # if(!require(devtools)) install.packages("devtools")
      # devtools::install_github("Niinleslie/MesKit")
      library("MesKit")
      library(ComplexHeatmap)
      library(viridis)
      trimatrix_sample<-list()
      trimatrix_sample[["P"]][["tri_matrix"]]<-data.frame(t(escc_sca_mm))
      colnames(trimatrix_sample[["P"]][["tri_matrix"]])<-c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
      trimatrix_sample[["P"]][["tsb.label"]] <- data.frame(Tumor_Sample_Barcode=unique(Maf_data$mut_cluster),Tumor_Sample_Label=unique(Maf_data$mut_cluster))
      
      len <- length(plotMutSigProfile(trimatrix_sample))
      for (kk in 1:len) {
        fil_name <- paste0("E:/数据库/STMut/数据处理/结果/inter/",dataset_name,"/",sample_name,"/","mut_signature_", kk, ".png")
        p1 <- plotMutSigProfile(trimatrix_sample)[[kk]]
        png(fil_name, width = 1225, height = 650)
        print(p1)
        dev.off()
      }
      # plotMutSigProfile(trimatrix_sample)[[1]]
      # plotMutSigProfile(trimatrix_sample)[[2]]
      # plotMutSigProfile(trimatrix_sample)[[3]]
      # plotMutSigProfile(trimatrix_sample)[[4]]
      # plotMutSigProfile(trimatrix_sample)[[5]]
      # plotMutSigProfile(trimatrix_sample)[[6]]
      
      ######################################空间分布簇突变集与COSMIC特征相关性######################################
      ##每簇突变集与COSMIC特征相关性
      fit_sample <- fitSignatures(trimatrix_sample,signaturesRef="exome_cosmic_v3")
      png(paste0("E:/数据库/STMut/数据处理/结果/inter/",dataset_name,"/",sample_name,"/","Cosine_similarity.png"), width = 1225, height = 650)
      p2 <- ComplexHeatmap::Heatmap(fit_sample$P$cosine.similarity, name = "Cosine similarity",col = viridis(100))
      print(p2)
      dev.off()
      
      fit_sample_cos_similarity <- as.data.frame(fit_sample$P$cosine.similarity)
      fit_sample_cos_similarity$cluster <- rownames(fit_sample_cos_similarity)
      library(tidyverse)
      # 使用gather函数将宽数据转换为长数据
      long_data <- fit_sample_cos_similarity %>% 
        gather(key = "COSMIC Signatures", value = "Cosine Similarity", -cluster)
      SBS_anno <- read.table("E:/数据库/STMut/数据处理/结果表格/SBS_anno.txt", header = T, sep = "\t", check.names = F)
      Cosmic_signature_anno <- merge(long_data, SBS_anno)
      Cosmic_signature_anno$`Cosine Similarity` <- round(Cosmic_signature_anno$`Cosine Similarity`, 4)
      
      write.table(Cosmic_signature_anno, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/Cosmic_signature_anno.txt"), sep = "\t", quote = F, row.names = F)
      ######################################每簇突变集的COSMIC特征贡献【每簇贡献大于0.08的特征######################################
      #trimatrix_sample <- triMatrix(tree.NJ)
      fit_sample <- fitSignatures(trimatrix_sample,signaturesRef="exome_cosmic_v3")
      data_df<-as.data.frame(fit_sample[["P"]][["contribution.relative"]])
      data_df$Category <- rownames(data_df)
      data_long <- tidyr::gather(data_df, key = "Signature", value = "Value", -Category)
      #删除小于0.08的特征
      data_long<- data_long[which(data_long$Value>=0.08),]
      data_long$Signature<-factor(data_long$Signature,levels = sort(unique(data_long$Signature),decreasing = T))
      my_colors <- c("grey","#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7", "#009B77","#FFB6C1", "#7FFFD4", "#FFD700", "#40E0D0")
      # 创建分组百分比堆叠图
      library(ggplot2)
      png(paste0("E:/数据库/STMut/数据处理/结果/inter/",dataset_name,"/",sample_name,"/mut_sigature_contribution.png"), width = 1225, height = 650)
      p3 <- ggplot(data_long, aes(x = Category, y = Value, fill = Signature)) +
        geom_bar(stat = "identity",color = "black") +
        scale_fill_manual(values = my_colors) +
        labs(title = "Mutational signature contribution",  y = "Mutational signature
       contribution") +
        theme_classic()+
        coord_flip()  # 翻转坐标轴
      print(p3)
      dev.off()
      
      #*************************************数据库需要的数据格式（最后把所有样本的合到一起）*****************************************************
      ###堆叠柱状图-----------------------------------------------------------------------------
      #将data_long的前两列两两组合，组成新的数据框
      length(unique(data_long$Signature))
      length(unique(unique(data_long$Category)))
      # 使用expand.grid生成所有可能的组合
      combinations <- expand.grid(unique(data_long$Category), unique(data_long$Signature))
      # 给数据框的列命名
      names(combinations) <- c("Category", "Signature")
      
      data_long_1 <- merge(combinations, data_long, all.x = T)
      data_long_1[is.na(data_long_1)] <- 0
      data_long_1$Value <- round(data_long_1$Value, 4)
      
      # 使用gsub()函数提取C_后面的数字，并将结果转换为数字类型
      data_long_1$numeric_order <- as.numeric(gsub("C_", "", data_long_1$Category))
      # 根据提取的数字对数据框进行排序
      data_long_2 <- data_long_1[order(data_long_1$numeric_order), ]
      # 如果不需要保留辅助排序列，可以将其删除
      data_long_2$numeric_order <- NULL
      
      # 使用paste()函数和sapply()函数将第一列的值转换为字符向量
      formatted_col1 <- sapply(unique(data_long_2$Category), function(x) paste0("'", x, "'"))
      # 将字符向量转换为单个字符串，格式为['C_0','C_1','C_2','C_3','C_4','C_5']
      formatted_string <- paste(formatted_col1, collapse = ", ")
      # 添加方括号[]来创建最终格式
      formatted_string <- paste0("[", formatted_string, "]")
      
      df_yaxix <- data.frame(Signature = 'yAxis', concatenated_values = formatted_string)
      
      # 按照第二列分组，并将第三列的值拼接成指定格式
      df_grouped <- data_long_2 %>%
        group_by(Signature) %>%
        summarise(concatenated_values = paste0("[", paste(Value, collapse = ","), "]"))
      
      df <- rbind(df_yaxix, df_grouped)
      df$sample <- sample_name
      df$dataset <- dataset_name
      
      write.table(df, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/MutSigContribution.txt"), sep = "\t", quote = F, row.names = F)
      
      ###热图------------------------------------------------------------------------------
      #去掉最后一列
      fit_sample_cos_similarity <- fit_sample_cos_similarity[, -ncol(fit_sample_cos_similarity)]
      
      df_1 <- sapply(rownames(fit_sample_cos_similarity), function(x) paste0("'", x, "'"))
      df_2 <- paste(df_1, collapse = ", ")
      #df_3 <- paste0("[", df_2, "]")
      mut_cluster <- df_2
      
      df_4 <- sapply(colnames(fit_sample_cos_similarity), function(x) paste0("'", x, "'"))
      df_5 <- paste(df_4, collapse = ", ")
      #df_6 <- paste0("[", df_5, "]")
      signature <- df_5
      
      # 将数据框中的所有数值列四舍五入到四位小数
      df_rounded <- fit_sample_cos_similarity
      numeric_columns <- sapply(df_rounded, is.numeric)  # 找出是数值的列
      df_rounded[numeric_columns] <- lapply(df_rounded[numeric_columns], round, 4)
      
      df_7 <- NULL
      for (i in 1:nrow(df_rounded))  {
        for (j in 1:ncol(df_rounded)) {
          tmp <- paste0("[", i-1, ",", j-1, ",", df_rounded[i,j], "],")
          df_7 <- paste0(df_7, tmp)
        }
      }
      
      SignatureSimilarities <- data.frame(
        MutCluster <- mut_cluster,
        Signature <- signature,
        Value <- df_7
      )
      SignatureSimilarities$sample <- sample_name
      SignatureSimilarities$dataset <- dataset_name
      names(SignatureSimilarities) <- c('MutCluster', 'Signature', 'Value', 'sample', 'dataset')
      
      write.table(SignatureSimilarities, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/SignatureSimilarities.txt"), sep = "\t", quote = F, row.names = F)
      
      ##柱状图-----------------------------------------------------------------------------------------
      MutProb <- trimatrix_sample$P$tri_matrix
      
      tmp_1 <- sapply(colnames(MutProb), function(x) paste0("'", x, "'"))
      tmp_2 <- paste(tmp_1, collapse = ", ")
      
      MutProb_1 <- data.frame()
      for (k in 1:nrow(MutProb)) {
        tmp_3 <- paste(MutProb[k,], collapse = ",")
        MutProb_1[k,1] <- tmp_2
        MutProb_1[k,2] <- tmp_3
        MutProb_1[k,3] <- rownames(MutProb)[k]
      }
      MutProb_1$sample <- sample_name
      MutProb_1$dataset <- dataset_name
      names(MutProb_1) <- c('Mutational_type', 'Mutation_prob', 'Mut_cluster', 'sample', 'dataset')
      
      write.table(MutProb_1, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/MutationSignatures.txt"), sep = "\t", quote = F, row.names = F)
      
      #清除变量
      # 列出所有的对象
      objs <- ls()
      
      # 保留aa变量，移除其他所有变量
      objs <- objs[which(objs != c("gse_dirs", "gse_dir", "gsm_dirs", "gsm_dir", "dataset_name", "input_path"))]
      
      # 使用rm函数移除这些变量以释放内存
      rm(list=objs)
      gc()
    }
}
