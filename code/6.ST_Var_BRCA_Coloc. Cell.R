#载入所需的R包；
library(Seurat)#
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r) #

input_path <- "E:/数据库/STMut/数据处理/结果/inter"

# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    
    brain <- readRDS(paste0(gsm_dir, "/ST_Mut_SeuratObject.rds"))
    spot_gene_mutcount_matrix <- read.table(paste0(gsm_dir, "/spot_gene_mutcount_matrix.txt"), sep = "\t", header = T)
    
    ###########################5-展示scsnv突变聚类的每个类别的共现细胞突变###############################
    library("gridExtra")
    #SORC数据库的去卷积结果
    spot_cell <- read.table("E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/GSM6177599_slice1_Deconvolution.txt",sep="\t",header = T, check.names = F)
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
      
      # co_local <- "E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/"
      # png(paste0(co_local, "co-localization",mut_cluster_i , ".png"), width = 1225, height = 650)
      # print(p_plot)
      # dev.off()
      
      eval(parse(text=paste("p",mut_cluster_i,"<-p_plot",sep="")))
      Cell_Colocalization_Mut_p_2$cluster <- mut_cluster_i
      Cell_Colocalization_Mut_p_cluster <- rbind(Cell_Colocalization_Mut_p_cluster,Cell_Colocalization_Mut_p_2)
    }
    eval(parse(text=paste("plist<-list(",paste("p",0:mut_cluster_i,sep="",collapse = ","),")",sep="")))
    library("grid")
    p1 <- grid.arrange(grobs = plist, ncol = 3)
    print(p1)
    # png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/co-localization.png", width = 1225, height = 650)
    # print(p1)
    # dev.off()
    
    coloc <- Cell_Colocalization_Mut_p_cluster[,-c(5,6)]
    coloc$Fisher_p <- round(coloc[,4], 4)
    coloc$sample <- sample_name
    coloc$dataset <- dataset_name
    write.table(coloc, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/coloc.txt"), sep = "\t", quote = F, row.names = F)
    
    #co_local <- "E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/"
    for(mut_cluster_i in unique(brain@meta.data$mut_cluster)){
      #Mutation-related co-localization 修改后的
      colocal <- ggplot(Cell_Colocalization_Mut_p_cluster[Cell_Colocalization_Mut_p_cluster$cluster%in%c(mut_cluster_i),], aes(x = cell_A, y = Mut_G, fill = Co)) +
        geom_tile()+
        #geom_tile(aes(width = Colocalization/100, height = Colocalization/100)) +
        scale_fill_gradient(low = "white", high = "#FF5500") +
        geom_text(aes(label=Colocalization),col ="black",size = 5) +
        labs(title = paste0("mut_cluster: ",mut_cluster_i), x = "Spot.Cell", y = "Spot.Mut_Gene",fill = "-logp") +
        theme_set(theme_bw()) +
        theme(title=element_text(size=16),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14), 
              axis.text.y = element_text(size=12,color = 'black'),
              axis.text.x = element_text(size=12,color = 'black',angle = 90, hjust = 1))
      # png(paste0(co_local, "co-localization",mut_cluster_i , ".png"), width = 1225, height = 650)
      # print(colocal)
      # dev.off()
    }
    
    #热图
    coloc_sorted <- coloc[order(coloc[,5]), ]
    coloc_sorted$Fisher_p <- as.numeric(coloc_sorted$Fisher_p)
    # 使用ifelse()函数进行条件赋值
    coloc_sorted$text <- ifelse(coloc_sorted$Fisher_p <= 0.05, "*", "")
    
    coloc_result <- data.frame()
    for (i in unique(coloc_sorted$cluster)) {
      tmp_coloc <- coloc_sorted[which(coloc_sorted$cluster == i),]
      tmp_coloc$log_coloc <- log(tmp_coloc$Colocalization+1, base = 10)
      tmp_coloc$log_coloc <- round(tmp_coloc$log_coloc, 4)
      
      tmp_coloc_sorted <- tmp_coloc[order(tmp_coloc[,2]), ]
      
      df_1 <- sapply(unique(tmp_coloc_sorted$Mut_G), function(x) paste0("'", x, "'"))
      df_2 <- paste(df_1, collapse = ",")
      Mut_gene <- df_2
      
      df_3 <- sapply(unique(tmp_coloc_sorted$cell_A), function(x) paste0("'", x, "'"))
      df_4 <- paste(df_3, collapse = ",")
      Cell_type <- df_4
      
      value <- ""
      text <- ""
      for (m in unique(tmp_coloc_sorted$Mut_G)) {
        tmp_data <- tmp_coloc_sorted[which(tmp_coloc_sorted$Mut_G == m),]
        
        #tmp_data_1 <- sapply(tmp_data$Colocalization, function(x) paste0("'", x, "'"))
        tmp_data_2 <- paste(tmp_data$Colocalization, collapse = ",")
        tmp_data_3 <- paste0("[", tmp_data_2, "],")
        
        value <- paste0(value, tmp_data_3)
        
        text_1 <- sapply(tmp_data$text, function(x) paste0("'", x, "'"))
        text_2 <- paste(text_1, collapse = ",")
        text_3 <- paste0("[", text_2, "],")
        
        text <- paste0(text, text_3)
      }
      
      index <- as.numeric(i)+1
      
      coloc_result[index,1] <- Mut_gene
      coloc_result[index,2] <- Cell_type
      coloc_result[index,3] <- value
      coloc_result[index,4] <- text
      coloc_result[index,5] <- i
      
    }
    
    names(coloc_result) <- c('Mut_gene', 'Cell_type', 'Value', 'remark_text', 'Mut_cluster')
    coloc_result$sample <- 'GSM6177599'
    coloc_result$dataset <- 'GSE203612'
    
    write.table(coloc_result, "E:/数据库/STMut/数据处理/结果表格/Coloc_heatmap.txt", sep = "\t", quote = F, row.names = F)
    
  }
}
