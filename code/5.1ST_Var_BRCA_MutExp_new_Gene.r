
#**************************************************************
#2024/7/2 11:14开始跑 18:30结束
#**************************************************************

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

for (yj in 1:40) {
  if(yj == 4){
    #GSE158704删除
  }else{
    
    gse_dir <- gse_dirs[yj]
    # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
    gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
    gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
    
    dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
    
    for(kkk in 1:length(gsm_dirs)) {
      
      gsm_dir <- gsm_dirs[kkk]
      sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
      
      ###############################5-spot类别差异变异基因###############################
      brain <- readRDS(paste0(gsm_dir, "/ST_Mut_SeuratObject.rds"))
      spot_gene_mutcount <- read.table(paste0(gsm_dir, "/spot_gene_mutcount_matrix.txt"), sep = "\t", header = T, check.names = F)
      spot_gene_mutcount_matrix_raw <- spot_gene_mutcount
      ###基因突变位点数突变谱下###
      library(reshape2)
      rownames(spot_gene_mutcount_matrix_raw) <- spot_gene_mutcount_matrix_raw[,1]
      spot_gene_mutcount_matrix_raw[is.na(spot_gene_mutcount_matrix_raw)] <- 0
      spot_gene_mutcount_matrix_raw <- as.data.frame(t(spot_gene_mutcount_matrix_raw[,-1]))
      ###############################6-spot类别变异调控基因差异表达###############################
      diff.results_cluster <- data.frame()
      
      #查看突变谱和表达谱的大小---------------------
      dim(spot_gene_mutcount_matrix_raw)#2888*750
      dim(brain@assays$Spatial$counts)#36601*750
      #查看基因名字，是否有-变成了.
      intersect(rownames(spot_gene_mutcount_matrix_raw), "C8orf37-AS1")#有交集，基因名没有变化
      intersect(rownames(brain@assays$Spatial$counts), "C8orf37-AS1")#有交集，基因名没有变化
      
      #取交集基因
      inter_G <- intersect(rownames(spot_gene_mutcount_matrix_raw),rownames(brain@assays$Spatial$counts))
      #Exp_data_all <- as.data.frame(brain@assays$Spatial$counts)
      Exp_data_all <- as.matrix(brain@assays$Spatial$counts)
      for(cluster_jj in unique(brain@meta.data$mut_cluster)){
        #计算数据框Mut_data的每行和数据框Exp_data的对应行（相同行号）的相关性和显著性
        ##cluster_jj的spot名
        cluster_jj_spot <- rownames(brain@meta.data[brain@meta.data$mut_cluster %in% cluster_jj,])
        #取inter_G+cluster_jj_spot
        Mut_data <- data.frame(t(spot_gene_mutcount_matrix_raw[inter_G,cluster_jj_spot]))
        Exp_data <- data.frame(t(Exp_data_all[inter_G,cluster_jj_spot]))
        
        names(Mut_data) <- inter_G
        names(Exp_data) <- inter_G
        
        #转置的时候，基因名改变了
        # intersect(colnames(Mut_data), "C8orf37-AS1")
        # intersect(colnames(Mut_data), "C8orf37-AS1")
        # intersect(colnames(Mut_data), "C8orf37.AS1")
        # intersect(colnames(Mut_data), "C8orf37.AS1")
        
        # 秩和检验
        diff_test <- mapply(function(x, y) {
          if(length(which(x>0))*length(which(x==0))>0){#如果这个基因在所有的spot中都突变或者都非突变就不做秩和检验
            p_value <- wilcox.test(y[which(x>0)], y[which(x==0)])$p.value#突变的spot的表达值和非突变的spot的表达值做了秩和检验
            fc <- mean(y[which(x>0)])/mean(y[which(x==0)])#
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
      
      df_selected_lab <- diff.results_cluster
      
      length(intersect(df_selected_lab$mut_G, rownames(spot_gene_mutcount_matrix_raw)))
      length(intersect(df_selected_lab$mut_G, "C8orf37-AS1")) 
      
      # 画分组的图
      #library(reshape2)
      Exp_muttype_cluster <- data.frame()
      for(cluster_jj in unique(df_selected_lab$mut_cluster)){
        #select数据框Mut_data的每行和数据框Exp_data的对应行（相同行号）的data
        ##cluster_jj的spot名
        cluster_jj_spot <- rownames(brain@meta.data[brain@meta.data$mut_cluster %in% cluster_jj,])
        ##cluster_jj的gene名
        cluster_jj_gene <- df_selected_lab[df_selected_lab$mut_cluster %in% cluster_jj,]$mut_G
        
        #length(intersect(cluster_jj_gene, rownames(spot_gene_mutcount_matrix_raw)))
        
        Mut_data <- data.frame(t(spot_gene_mutcount_matrix_raw[cluster_jj_gene,cluster_jj_spot]))
        names(Mut_data) <- cluster_jj_gene
        Mut_data$spot <- rownames(Mut_data)
        Mut_data_melt <- melt(Mut_data, id.vars = "spot", variable.name="mut_G",value.name="mut_type")
        Mut_data_melt$mut_type <- ifelse(Mut_data_melt$mut_type>0,1,0)
        
        Exp_data_all_bridge <- as.data.frame(Exp_data_all)##class(Exp_data_all)是"matrix" "array",需要转为data.frame然后再转置
        Exp_data <- data.frame(t(Exp_data_all_bridge[cluster_jj_gene,cluster_jj_spot]))
        names(Exp_data) <- cluster_jj_gene
        Exp_data$spot <- rownames(Exp_data)
        Exp_data_melt <- melt(Exp_data, id.vars = "spot", variable.name="mut_G",value.name="G_Exp")
        
        Exp_muttype_cluster_bridge <- cbind(Exp_data_melt,mut_type=Mut_data_melt$mut_type)
        Exp_muttype_cluster_bridge$cluster <- cluster_jj
        Exp_muttype_cluster <- rbind(Exp_muttype_cluster,Exp_muttype_cluster_bridge)
      }
      intersect(Exp_muttype_cluster$mut_G, "C8orf37.AS1")
      intersect(Exp_muttype_cluster$mut_G, "C8orf37-AS1")
      
      #**************************************数据库需要的数据（最后把所有样本的合到一起）*********************************************
      similarity_boxplot <- data.frame()
      Exp_muttype_cluster_sorted <- Exp_muttype_cluster[order(Exp_muttype_cluster[,5]), ]
      Exp_muttype_cluster_sorted$cluster_new <- paste("C_", Exp_muttype_cluster_sorted$cluster, sep = "")
      
      for (curr_gene in unique(Exp_muttype_cluster_sorted$mut_G)) {
        tmp_1 <- Exp_muttype_cluster_sorted[which(Exp_muttype_cluster_sorted$mut_G == curr_gene),]
        mut <- tmp_1[which(tmp_1$mut_type == "1"),]
        non_mut <- tmp_1[which(tmp_1$mut_type == "0"),]
        
        mut_ordered <- mut[order(mut[,5]),]
        non_mut_ordered <- non_mut[order(non_mut[,5]),]
        
        mut_cluster <- paste(paste0("'",mut_ordered$cluster_new,"'"), collapse = ",")
        mut_exp <- paste(mut_ordered$G_Exp, collapse = ",")
        nonMut_cluster <- paste(paste0("'",non_mut_ordered$cluster_new,"'"), collapse = ",")
        nonMut_exp <- paste(non_mut_ordered$G_Exp, collapse = ",")
        
        tmp_table <- data.frame(
          Mut_cluster = mut_cluster,
          Mut_exp = mut_exp,
          NonMut_cluster = nonMut_cluster,
          NonMut_exp = nonMut_exp,
          Mut_gene = curr_gene
        )
        
        similarity_boxplot <- rbind(similarity_boxplot, tmp_table)
      }
      
      similarity_boxplot$sample <- sample_name
      similarity_boxplot$dataset <- dataset_name
      
      write.table(similarity_boxplot, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/similarity_boxplot_gene.txt"), sep = "\t", quote = F, row.names = F)
      
    }
  }
  
}

