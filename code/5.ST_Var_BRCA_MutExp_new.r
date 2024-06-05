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
    
    #gsm_dir <- gsm_dirs[i]
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    
    ###############################5-spot类别差异变异基因###############################
    brain <- readRDS(paste0(gsm_dir, "/ST_Mut_SeuratObject.rds"))
    spot_gene_mutcount <- read.table(paste0(gsm_dir, "/spot_gene_mutcount_matrix.txt"), sep = "\t", header = T)
    spot_gene_mutcount_matrix_raw <- spot_gene_mutcount
    ###基因突变位点数突变谱下###
    library(reshape2)
    rownames(spot_gene_mutcount_matrix_raw) <- spot_gene_mutcount_matrix_raw[,1]
    spot_gene_mutcount_matrix_raw[is.na(spot_gene_mutcount_matrix_raw)] <- 0
    spot_gene_mutcount_matrix_raw <- as.data.frame(t(spot_gene_mutcount_matrix_raw[,-1]))
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
    png(paste0(gsm_dir, "/similarity_boxplot.png"), width = 1225, height = 650)
    p1 <- ggplot(Exp_muttype_cluster,aes(mut_G,G_Exp))+
      geom_boxplot( aes(color = as.character(mut_type)))+
      scale_color_manual(values = c("#00AFBB","#E7B800"))+
      stat_compare_means(aes(group = mut_type),label = "p.signif")+
      facet_wrap(~ as.factor(cluster),ncol = 3,scales = "free")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      labs(x = "Mutant gene", y = "Gene expression", color = "Mutant type")
    print(p1)
    dev.off()
    
    diff.results_cluster[,c(1,2)] <- round(diff.results_cluster[,c(1,2)], 4)
    write.table(diff.results_cluster, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/similarity_table.txt"), sep = "\t", quote = F, row.names = F)
    
    
    #**************************************数据库需要的数据（最后把所有样本的合到一起）*********************************************
    similarity_boxplot <- data.frame()
    Exp_muttype_cluster_sorted <- Exp_muttype_cluster[order(Exp_muttype_cluster[,5]), ]
    for (i in unique(Exp_muttype_cluster_sorted$cluster)) {
      tmp_1 <- Exp_muttype_cluster_sorted[which(Exp_muttype_cluster_sorted$cluster == i),]
      mut <- tmp_1[which(tmp_1$mut_type == "1"),]
      non_mut <- tmp_1[which(tmp_1$mut_type == "0"),]
      
      mut_ordered <- mut[order(mut[,2]),]
      non_mut_ordered <- non_mut[order(non_mut[,2]),]
      
      mut_gene <- paste(paste0("'",mut_ordered$mut_G,"'"), collapse = ",")
      mut_exp <- paste(mut_ordered$G_Exp, collapse = ",")
      nonMut_gene <- paste(paste0("'",non_mut_ordered$mut_G,"'"), collapse = ",")
      nonMut_exp <- paste(non_mut_ordered$G_Exp, collapse = ",")
      
      similarity_boxplot[i,1] <- mut_gene
      similarity_boxplot[i,2] <- mut_exp
      similarity_boxplot[i,3] <- nonMut_gene
      similarity_boxplot[i,4] <- nonMut_exp
      similarity_boxplot[i,5] <- i
    }
    names(similarity_boxplot) <- c('Mut_gene', 'Mut_exp', 'NonMut_gene', 'NonMut_exp', 'Mut_cluster')
    
    similarity_boxplot$sample <- sample_name
    similarity_boxplot$dataset <- dataset_name
    
    write.table(similarity_boxplot, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/similarity_boxplot.txt"), sep = "\t", quote = F, row.names = F)
    
  }
}
