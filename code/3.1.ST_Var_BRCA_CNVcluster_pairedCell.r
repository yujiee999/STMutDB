library(Seurat)
library(ggplot2)
library(patchwork)#
library(dplyr)
library(hdf5r) #
library(tidyverse)
library(data.table)
library(readxl)

# 中间文件
input_inter_path <- "E:/数据库/STMut/数据处理/结果/inter"

# 配对的单细胞
input_singlecell_path <- "E:/数据库/STMut/数据处理/结果/paired_single_cell"
  
# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_inter_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    #************************提取spot表达谱***************************
    # 读取空转的seurat对象
    brain <- readRDS(paste0(gsm_dir, "/ST_seuratObject.rds"))
    spot_exp_mat <- brain@assays[["Spatial"]]@layers[["counts"]]
    spot_exp_G <- brain@assays[["Spatial"]]@features %>% rownames()
    spot_exp_spot <- brain@assays[["Spatial"]]@cells %>% rownames()
    colnames(spot_exp_mat) <- spot_exp_spot
    rownames(spot_exp_mat) <- spot_exp_G
    
    #************************提取细胞表达谱***************************
    #只选取前1000个免疫细胞作为参考
    
    # 使用strsplit函数分割字符串
    split_string <- strsplit(gsm_dir, "/")
    # 提取最后一个元素
    sample_name <- tail(unlist(split_string), n=1)
    # 读取每个样本对应的单细胞数据文件
    paired_singelcell <- read_excel("E:/数据库/STMut/数据处理/结果/paired_single_cell/paired_singlecell_data.xlsx")
    # 对应的文件夹的名字
    singlecell_datst <- paired_singelcell[which(paired_singelcell$Sample == sample_name),2]
    # 对应的路径
    tmp_path <- paste0(input_singlecell_path, "/", singlecell_datst, "/")
    
    if(singlecell_datst %in% c("BRCA", "Childhood_brain_tumor", "CRC", "CRC_liver_metastases", "CSCC", "GBM", "Gsatrointestinal_stromal_tumor", "LIHC", "MIBC", "OV", "OVCA", "PDAC")){
      
      sc_data <- read.table(paste0(tmp_path, "SC_data.txt"), header = T, sep = "\t")
      sc_celltype <- read.table(paste0(tmp_path, "SC_CellType.txt"), header = T, sep = "\t", row.names = 1)
      sc_celltype_immune <- sc_celltype[-which(sc_celltype$cell_type_final == "Tumor"),]$cell_name
      if(length(sc_celltype_immune)>=1000){
        sc_celltype_immune_1000 <- sc_celltype_immune[1:1000]
      }else{
        sc_celltype_immune_1000 <- sc_celltype_immune
      }
      
      seurat_object <- CreateSeuratObject(counts = sc_data)
      seurat_object <- subset(seurat_object, cells = sc_celltype_immune_1000)
      cell_exp_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
      cell_exp_G <- seurat_object@assays[["RNA"]]@features %>% rownames()
      cell_exp_spot <- seurat_object@assays[["RNA"]]@cells %>% rownames()
      colnames(cell_exp_mat) <- cell_exp_spot
      rownames(cell_exp_mat) <- cell_exp_G
      
    }else if(singlecell_datst %in% c("HNSCC_GSE103322")){
      
      sc_data <- readRDS(paste0(tmp_path, "/", singlecell_datst, ".rds"))
      sc_celltype_immune <- rownames(sc_data@meta.data[-which(sc_data$celltype %in% c("Tumor cell", "Fibroblast", "Endothelial cell")),])
      if(length(sc_celltype_immune)>=1000){
        sc_celltype_immune_1000 <- sc_celltype_immune[1:1000]
      }else{
        sc_celltype_immune_1000 <- sc_celltype_immune
      }
      
      seurat_object <- subset(sc_data, cells = sc_celltype_immune_1000)
      cell_exp_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
      cell_exp_G <- seurat_object@assays[["RNA"]]@features %>% rownames()
      cell_exp_spot <- seurat_object@assays[["RNA"]]@cells %>% rownames()
      colnames(cell_exp_mat) <- cell_exp_spot
      rownames(cell_exp_mat) <- cell_exp_G
      
    }else if(singlecell_datst %in% c("SCC_GSE123813")){
      
      sc_data <- readRDS(paste0(tmp_path, "/", singlecell_datst, ".rds"))
      sc_celltype_immune_1000 <- rownames(sc_data)[1:1000]
      
      seurat_object <- subset(sc_data, cells = sc_celltype_immune_1000)
      cell_exp_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
      cell_exp_G <- seurat_object@assays[["RNA"]]@features %>% rownames()
      cell_exp_spot <- seurat_object@assays[["RNA"]]@cells %>% rownames()
      colnames(cell_exp_mat) <- cell_exp_spot
      rownames(cell_exp_mat) <- cell_exp_G
      
    }else{
      
      # 列出所有以.txt .csv结尾的文件
      txt_files <- list.files(path = tmp_path, pattern = "\\.txt$")
      csv_files <- list.files(path = tmp_path, pattern = "\\.csv$")
      
      sc_data <- read.table(paste0(tmp_path, "/", txt_files[1]), sep = ",",header = T, row.names = 1)
      sc_celltype <- read.csv(paste0(tmp_path, "/", csv_files[1]), row.names = 1)
      
      sc_celltype_immune <- rownames(sc_celltype[which(sc_celltype$CT != "Tumor cell"),])
      if(length(sc_celltype_immune)>=1000){
        sc_celltype_immune_1000 <- sc_celltype_immune[1:1000]
      }else{
        sc_celltype_immune_1000 <- sc_celltype_immune
      }
      
      seurat_object <- CreateSeuratObject(counts = sc_data)
      seurat_object <- subset(seurat_object, cells = sc_celltype_immune_1000)
      cell_exp_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
      cell_exp_G <- seurat_object@assays[["RNA"]]@features %>% rownames()
      cell_exp_spot <- seurat_object@assays[["RNA"]]@cells %>% rownames()
      colnames(cell_exp_mat) <- cell_exp_spot
      rownames(cell_exp_mat) <- cell_exp_G
      
    }
    
    #************************整合细胞和spot的表达谱***************************
    inter_G <- intersect(rownames(spot_exp_mat),rownames(cell_exp_mat))
    exp_mat_1 <- cbind(spot_exp_mat[inter_G,],cell_exp_mat[inter_G,])
    cell_type <- data.frame(Type=c(rep("spot",ncol(spot_exp_mat)),rep("cell",ncol(cell_exp_mat))),row.names = colnames(exp_mat_1))
    
    library(infercnv)
    #第一步，根据上述的三个文件创建inferCNV对象
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exp_mat_1, # 可以直接提供矩阵对象
                                        annotations_file=cell_type,
                                        delim="\t",
                                        gene_order_file="E:/数据库/STMut/suminghai/slice1_gene_ordering_file.txt",
                                        ref_group_names="cell",
                                        min_max_counts_per_cell=c(0,+Inf))
    
    setwd(paste0(input_inter_path, "/", dataset_name, "/", sample_name))
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
    
    
  }
  
  getwd()
  
}


library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    infercnv_obj = readRDS(paste0(input_inter_path, "/", dataset_name, "/", sample_name, "/output_dir/", "run.final.infercnv_obj"))
    expr <- infercnv_obj@expr.data
    gn <- rownames(expr)
    geneFile <- read.table("E:/数据库/STMut/suminghai/slice1_gene_ordering_file.txt",header = F,sep = "\t",stringsAsFactors = F)
    rownames(geneFile) = geneFile$V1
    sub_geneFile <- geneFile[intersect(gn,geneFile$V1),]
    expr = expr[intersect(gn,geneFile$V1),]
    head(sub_geneFile,4)
    expr[1:4,1:4]
    #聚类，7类，提取结果
    set.seed(20240603)
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
    
    output_file <- paste0("E:/数据库/STMut/数据处理/结果/result/", dataset_name, "/", sample_name, "/cnv.png")
    png(output_file, width = 674, height = 674)
    print(ht)
    dev.off()
    
  }
  
}
