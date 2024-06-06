cellcell_LR<-function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
                       pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
                       sort.by.source.priority = TRUE, color.heatmap = c("Spectral", 
                                                                         "viridis"), n.colors = 10, direction = -1, thresh = 0.05, 
                       comparison = NULL, group = NULL, remove.isolate = FALSE, 
                       max.dataset = NULL, min.dataset = NULL, min.quantile = 0, 
                       max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, 
                       color.text = NULL, title.name = NULL, font.size = 10, font.size.title = 10, 
                       show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
                       angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    }
    else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net", 
                                  sources.use = sources.use, targets.use = targets.use, 
                                  signaling = signaling, pairLR.use = pairLR.use, 
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, 
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), 
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, 
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    # df.net$pval[df.net$pval > 0.05] = 1
    # df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    # df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) * 
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2), 
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, 
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      signaling = signaling, pairLR.use = pairLR.use, 
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, 
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), 
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, 
                                       unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, 
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), 
                               each = length(levels(df.net$target))), levels(df.net$target), 
                           sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2", 
                              "source.target", "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, 
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) * 
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], 
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, 
                                dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)], 
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          }
          else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, 
    ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = interaction_name_2.order)
  }
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
                                                                    unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, 
                                                        df$target))
      df <- with(df, df[order(target, source), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      df <- with(df, df[order(source, target), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, 
                                                          df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target), ])
      }
      else {
        df <- with(df, df[order(target, source), ])
      }
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  return(df)
}
options(stringsAsFactors = FALSE)
library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)

input_path <- "E:/数据库/STMut/数据处理/结果/inter"
# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  #gse_dir <- gse_dirs[i]
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  dataset_name <- tail(unlist(strsplit(gse_dir, "/")), n=1)
  
  for(gsm_dir in gsm_dirs) {
    
    #gsm_dir <- gsm_dirs[i]
    sample_name <- tail(unlist(strsplit(gsm_dir, "/")), n=1)
    
    brain <- readRDS(paste0(gsm_dir, "/ST_Mut_SeuratObject.rds"))
    seurat_object <- readRDS(paste0(gsm_dir, "/MutSites_01_SeuratObject.rds"))
    spot_gene_mutcount_matrix <- read.table(paste0(gsm_dir, "/spot_gene_mutcount_matrix.txt"), sep = "\t", header = T, check.names = F)
    
    Brain_ST <- brain
    
    
    ##计算每簇spot突变最多的前10个基因。按照cluster分组，然后按照每个基因的表达水平从高到低排序，最后选择每个细胞群中的前10个基因作为代表性基因
    seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25)
    gene_top10_clutster <- seurat_object.markers %>% 
      group_by(cluster) %>%
      arrange(desc(pct.1)) %>%
      slice_head(n = 10)
    cluster_group_prob <- data.frame()
    ###对于spot簇C_0的突变基因gene1来说,C_0可以分为两簇
    for(cluster_jj in unique(brain@meta.data$mut_cluster)){
      #计算数据框Mut_data的每行和数据框Exp_data的对应行（相同行号）的相关性和显著性
      ##cluster_jj的spot名
      cluster_jj_spot <- rownames(brain@meta.data[brain@meta.data$mut_cluster %in% cluster_jj,])
      gene_bridge <- gene_top10_clutster[gene_top10_clutster$cluster %in% cluster_jj,]$gene
      for(gene_jj in  gene_bridge){
        cluster_jj_mut_spot <- spot_gene_mutcount_matrix[spot_gene_mutcount_matrix[,gene_jj]>0,]$CB
        cluster_jj_spot_mut <- intersect(cluster_jj_mut_spot,cluster_jj_spot)
        cluster_jj_spot_nonmut <- setdiff(cluster_jj_spot,cluster_jj_mut_spot)
        ############################### 1.数据输入，处理################################
        Brain_ST_cellchat <- brain
        #查看数据情况
        Brain_ST_cellchat@meta.data$celltype <- paste("C",Brain_ST_cellchat$mut_cluster,sep = "_")
        #基于基因突变重命名mut_cluster类别
        Brain_ST_cellchat@meta.data[cluster_jj_spot_mut,]<-"Mut"
        Brain_ST_cellchat@meta.data[cluster_jj_spot_nonmut,]<-"NonMut"
        Idents(Brain_ST_cellchat) <- "celltype"
        #head(Brain_ST_cellchat)
        #可定义颜色
        color.use <- scPalette(nlevels(Brain_ST_cellchat))
        names(color.use) <- levels(Brain_ST_cellchat)
        SpatialDimPlot(Brain_ST_cellchat, label = TRUE, label.size = 3, cols = color.use)
        ############################### 2.准备输入文件################################
        #矩阵信息
        data.input = Seurat::GetAssayData(Brain_ST_cellchat, slot = "data", assay = "SCT") 
        #meta信息
        meta = data.frame(labels = Idents(Brain_ST_cellchat), #名字自定义
                          row.names = names(Idents(Brain_ST_cellchat))) # manually create a dataframe consisting of the cell labels
        unique(meta$labels)
        # 空间图像信息
        spatial.locs = Seurat::GetTissueCoordinates(Brain_ST_cellchat, scale = NULL, 
                                                    cols = c("imagerow", "imagecol")) 
        # Scale factors and spot diameters 信息 
        scale.factors = jsonlite::fromJSON(txt = 
                                             file.path(paste0("E:/数据库/STMut/数据处理/结果/", dataset_name, "/", sample_name, "/spatial"), 'scalefactors_json.json'))
        scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                             fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
        )
        ############################### 3.创建CellChat对象################################
        cellchat <- createCellChat(object = data.input, 
                                   meta = meta, 
                                   group.by = "labels", #前面的meta ，定义的名字是labels
                                   datatype = "spatial", ###
                                   coordinates = spatial.locs, 
                                   scale.factors = scale.factors)
        ############################### 4.设置参考数据库################################
        CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
        
        # use a subset of CellChatDB for cell-cell communication analysis
        #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
        # use all CellChatDB for cell-cell communication analysis
        CellChatDB.use <- CellChatDB # simply use the default CellChatDB
        
        # set the used database in the object
        cellchat@DB <- CellChatDB.use
        ############################### 5.CellChat预处理################################
        # subset the expression data of signaling genes for saving computation cost
        cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
        future::plan("multisession", workers = 1) #笔记本可以选1
        ##识别过表达基因
        cellchat <- identifyOverExpressedGenes(cellchat)
        #识别过表达配体受体对
        cellchat <- identifyOverExpressedInteractions(cellchat)
        ############################### 6.推断细胞通讯网络################################
        cellchat <- computeCommunProb(cellchat, 
                                      type = "truncatedMean", trim = 0.1, 
                                      distance.use = TRUE, 
                                      scale.distance = 0.01, nboot=10)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        ############################### 7.计算cell-cell communication################################
        #计算每个信号通路相关的所有配体-受体相互作用的通信结果
        cellchat <- computeCommunProbPathway(cellchat)
        #计算整合的细胞类型之间通信结果
        cellchat <- aggregateNet(cellchat)
        # ############################### 8.根据使用netVisual_circle显示任意两个celltype之间的通讯次数################################
        # netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
        #                  weight.scale = T, label.edge= F, title.name = "Number of interactions")
        # ############################### 9.单个信号通路可视化########################################################################
        # pathways.show <- c("SPP1")
        # netVisual_aggregate(cellchat, signaling = pathways.show, 
        #                     edge.width.max = 18, vertex.size.max = 10, 
        #                     alpha.image = 0.2, vertex.label.cex = 3.5)#二维图
        # netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", 
        #                     edge.width.max = 8, vertex.size.max = 7, 
        #                     alpha.image = 0.4, vertex.label.cex = 7,point.size = 4)#空间的图
        # ############################### 10.计算中心性分数########################################################################
        # # Compute the network centrality scores
        # cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
        # # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
        # netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
        #                                   width = 8, height = 2.5, font.size = 10)
        ############################### 11.计算受体配体互作可能性prob,pvalue########################################################################
        ccc<-cellcell_LR(cellchat)
        ccc_1<-ccc[ccc$source%in%c("Mut","NonMut") | ccc$target%in%c("Mut","NonMut"),]
        ############################### 12.计算受体配体互作可能性差异prob_diff########################################################################
        #将字符串向量unique(ccc_1$source.target)元素分为两组，第一组元素以“Mut”开头或以“ Mut”结尾，剩下的为一组
        strings<-unique(as.character(ccc_1$source.target))
        group1 <- strings[grep("^Mut| Mut$", strings)]
        group2 <- strings[!(strings %in% group1)]
        #然后构建第一个数据框，第一列是group1的元素，第二列也是group1的元素，同时将第二列元素中以“Mut”开头或以“ Mut”结尾的字符串替换为"NonMut"。
        #构建第二个数据框，第一列是group2的元素，第二列也是group2的元素，同时将第一列元素中以“NonMut”开头或以“ NonMut”结尾的字符串替换为"Mut"。
        #最后合并两个数据框获得统一接口，包括group1.source.target和group2.source.target两列
        df1 <- data.frame(group1 = group1, group2 = group1)
        df1$group2 <- gsub("^Mut", "NonMut", df1$group2)
        df1$group2 <- gsub(" Mut$", " NonMut", df1$group2)
        df2 <- data.frame(group1 = group2, group2 = group2)
        df2$group1 <- gsub("^NonMut", "Mut", df2$group1)
        df2$group1 <- gsub(" NonMut$", " Mut", df2$group1)
        result <- rbind(df1, df2)
        result<-unique(result)
        #给突变前后两组互作数据添加统一接口列，设置给自的prob列
        colnames(result)<-c("group1.source.target","group2.source.target")
        group1_prob<-ccc_1[,c("source.target","ligand","receptor","interaction_name","interaction_name_2","pathway_name","annotation","evidence","prob")]
        colnames(group1_prob)[which(colnames(group1_prob)=="source.target")]<-"group1.source.target"
        colnames(group1_prob)[which(colnames(group1_prob)=="prob")]<-"group1.prob"
        group1_prob<-merge(result,group1_prob)
        group2_prob<-ccc_1[,c("source.target","ligand","receptor","interaction_name","interaction_name_2","pathway_name","annotation","evidence","prob")]
        colnames(group2_prob)[which(colnames(group2_prob)=="source.target")]<-"group2.source.target"
        colnames(group2_prob)[which(colnames(group2_prob)=="prob")]<-"group2.prob"
        group2_prob<-merge(result,group2_prob)
        #合并突变前后两组互作数据，并获得受体配体互作可能性差异prob_diff
        group_prob<-merge(group1_prob,group2_prob,all=T)
        group_prob[is.na(group_prob)]<-0#没有检测到的互作可能性设置为0
        group_prob$prob_diff<-group_prob$group1.prob-group_prob$group2.prob
        group_prob$mut_Gene<-gene_jj
        group_prob$mut_cluster<-paste0("C_",cluster_jj)
        cluster_group_prob <- rbind(group_prob,cluster_group_prob)
      }
    }
    cluster_group_prob[,c(10,11,12)] <- round(cluster_group_prob[,c(10,11,12)] ,4)
    write.table(cluster_group_prob, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/cluster_group_prob.txt"), sep = "\t", quote = F, row.names = F)
    
    #cluster_group_prob的interaction_type列基于group1.prob和group2.prob设置类别，group1.prob等于0并且group2.prob大于0设置为GAIN，group1.prob大于0并且group2.prob等于0设置为LOSS，group1.prob大于0并且group2.prob大于0并且group1.prob大于group2.prob设置为DECRESS，group1.prob大于0并且group2.prob大于0group1.prob小于group2.prob设置为INCRESS
    cluster_group_prob<-cluster_group_prob %>%
      mutate(
        interaction_type = case_when(
          group1.prob == 0 & group2.prob > 0 ~ "Gain",
          group1.prob > 0 & group2.prob == 0 ~ "Loss",
          group1.prob > 0 & group2.prob > 0 & group1.prob > group2.prob ~ "Decress",
          group1.prob > 0 & group2.prob > 0 & group1.prob < group2.prob ~ "Incress",
          TRUE ~ "Other"  # 其他情况设置为"OTHER"或根据需要进行调整
        )
      )
    ##分为两组展示，突变在发送信号的一组，突变在接受信号的二组
    #提取group1.source.target列元素以" Mut"结尾的行
    send_with_mut <- cluster_group_prob[grepl("^Mut", cluster_group_prob$group1.source.target), ]
    reci_with_mut <- cluster_group_prob[grepl(" Mut$", cluster_group_prob$group1.source.target), ]
    #互作的功能注释
    LR_class <- read.table("E:/数据库/STMut/数据处理/union_LR_database.txt", header=T,sep="\t",as.is=T,comment.char = "")
    
    ##发送信号的突变簇互作情况
    #每簇中突变基因影响的互作对的功能性获得丢失个数统计
    send_with_mut_static <- send_with_mut %>%
      distinct(mut_cluster,interaction_name_2,interaction_type) %>%
      dplyr::count(mut_cluster,interaction_type)
    
    ggplot()+
      geom_bar(data=send_with_mut_static,
               aes(x=interaction_type,y=n),
               stat = "identity",position = "stack",
               width=0.6,fill = "gray", color = "gray")+
      facet_grid(mut_cluster ~ ., scales = "free", space = "free")+
      theme_bw()+
      coord_flip() -> send_mut_statistic 
    
    png(paste0(gsm_dir, "/send_mut_statistic.png"), width = 1225, height = 650)
    print(send_mut_statistic)
    dev.off()
    
    ##统计每cluster的获得丢失的LR集合
    #1.对数据框send_with_mut进行操作，通过distinct()函数对mut_cluster、interaction_name和interaction_type这三列进行去重操作。
    #2.将interaction_name列按照"_"进行拆分，然后使用unnest()函数将拆分出的元素和其他列的元素构成新的行
    send_with_mut_GO_input <- send_with_mut %>%
      distinct(mut_cluster,interaction_name,interaction_type) %>%
      mutate(interaction_name = str_split(interaction_name, "_")) %>%
      unnest(cols = interaction_name)%>%
      distinct()
    
    
    library(org.Hs.eg.db) #人类注释数据???
    library(clusterProfiler)#进行GO富集和KEGG富集
    library(dplyr) #进行数据转换
    library(ggplot2)#绘图
    library(tidyverse)
    library(GOplot)
    enrich_result_all<-data.frame()
    for(mut_cluster_j in unique(send_with_mut_GO_input$mut_cluster)){
      send_with_mut_GO_input_1<-send_with_mut_GO_input[send_with_mut_GO_input$mut_cluster %in% mut_cluster_j,]
      for(type_i in unique(send_with_mut_GO_input_1$interaction_type)){
        # 创建一个geneList，包含interaction_name列的元素
        geneList <- send_with_mut_GO_input_1[send_with_mut_GO_input_1$interaction_type %in% type_i,]$interaction_name
        # 进行GO功能富集分析
        enrich_result <- enrichGO(gene = geneList, 
                                  keyType = "SYMBOL", 
                                  OrgDb = org.Hs.eg.db, 
                                  ont = "BP", 
                                  pvalueCutoff = 0.05, 
                                  pAdjustMethod = "BH", 
                                  qvalueCutoff = 0.2)
        enrich_result@result$interaction_type<-type_i
        enrich_result@result$mut_cluster<-mut_cluster_j
        enrich_result_all<-rbind(enrich_result_all,enrich_result@result)
      }
    }
    ##挑选
    enrich_result_all_top3<-enrich_result_all %>%
      group_by(mut_cluster,interaction_type) %>%
      arrange(desc(-log(qvalue))) %>%
      slice_head(n = 3)
    ##功能富集
    ggplot(enrich_result_all_top3, aes(interaction_type, Description)) +
      geom_point(aes(color=-log(qvalue), size=Count))+theme_bw()+
      theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size = 12))+
      scale_color_gradient(low='#6699CC',high='#CC3333')+
      labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
      scale_y_discrete(position = "left")+
      facet_grid(mut_cluster ~ ., scales = "free", space = "free") -> send_GO_enrich
    
    png(paste0(gsm_dir, "/send_GO_enrich.png"), width = 1225, height = 650)
    print(send_GO_enrich)
    dev.off()
    
    # library(aplot)
    # p1%>%
    #   insert_left(mut_cluster, width = .05)
    #每组展示10个突变基因扰动强度最大的3对互作的河流图。按照mut_Gene和mut_cluster分组，然后按照互作强度改变绝对值从高到低排序，最后选择前3个互作
    send_with_mut_top10 <- send_with_mut %>% 
      group_by(mut_Gene,mut_cluster) %>%
      arrange(desc(abs(prob_diff))) %>%
      slice_head(n = 3)
    library(tidyr)
    # 将group1.source.target列按照" -> "拆分为source和target两列
    send_with_mut_top10 <- separate(send_with_mut_top10, group1.source.target, into = c("source", "target"), sep = " -> ",remove = FALSE)
    # 将source列中包含"Mut"的元素替换为mut_cluster的元素
    send_with_mut_top10 <- send_with_mut_top10 %>%
      mutate(source = ifelse(grepl("Mut", source), mut_cluster, source))
    # 将target列中包含"Mut"的元素替换为mut_cluster的元素
    send_with_mut_top10 <- send_with_mut_top10 %>%
      mutate(target = ifelse(grepl("Mut", target), mut_cluster, target))
    send_with_mut_top10<-as.data.frame(send_with_mut_top10)
    # 加载ggalluvial包
    library(ggalluvial)
    #用数据框send_with_mut_top10画桑基图，第一列是mut_Gene，第二列是source，第三列是ligand，第五列是interaction_type，第六列是receptor，第七列是target
    # 绘制第一种桑基图
    ggplot(send_with_mut_top10, aes(axis1 = mut_Gene, axis2 = source, axis3 = ligand, axis4 = interaction_type, axis5 = receptor, axis6 = target)) +
      geom_alluvium(aes(fill = source)) + # 线的颜色用source填充
      geom_stratum(aes(fill = source)) + # 第二列块的填充色用source填充
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_fill_brewer(palette = "Set3") +  # 使用调色板
      theme_void()+
      scale_x_continuous(breaks = 1:6, labels = c("mut_Gene","source","ligand","interaction_type","receptor","target"))+
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(size = 12, face = "bold",vjust = 5))
    
    #用数据框send_with_mut_top10画桑基图，第一列是mut_Gene，第二列是source，第三列是ligand，第五列是interaction_type，第六列是receptor，第七列是target
    # 绘制另外一种桑基图
    ggplot(send_with_mut_top10, aes(axis1 = paste0(source," ",mut_Gene), axis3 = ligand, axis4 = interaction_type, axis5 = receptor, axis6 = target)) +
      geom_alluvium(aes(fill = paste0(source," ",mut_Gene))) + # 线的颜色用source," ",mut_Gene填充
      geom_stratum() + # 无填充色
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      #scale_fill_brewer(palette = "Set3") +  # 使用调色板
      theme_void()+
      scale_x_continuous(breaks = 1:6, labels = c("mut_Gene","source","ligand","interaction_type","receptor","target"))+
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(size = 12, face = "bold",vjust = 5))+
      labs(color = "Mutant gene") -> send_sankeyPlot
    
    png(paste0(gsm_dir, "/send_sankeyPlot.png"), width = 1225, height = 650)
    print(send_sankeyPlot)
    dev.off()
    
    send_with_mut_top10_static<-send_with_mut_top10 %>%
      distinct(mut_cluster,mut_Gene,interaction_name_2,pathway_name,interaction_type) %>%
      dplyr::count(mut_cluster,pathway_name,interaction_type)
    ggplot()+
      geom_bar(data=send_with_mut_top10_static,
               aes(x=interaction_type,y=n,fill=pathway_name),
               stat = "identity",position = "stack",width=0.6)+
      facet_grid(mut_cluster ~ ., scales = "free", space = "free")+
      theme_bw()+
      coord_flip()
    
    ##接收信号的突变簇互作情况
    #每簇中突变基因影响的互作对的功能性获得丢失个数统计
    reci_with_mut_static<-reci_with_mut %>%
      distinct(mut_cluster,as.character(interaction_name_2),interaction_type) %>%
      dplyr::count(mut_cluster,interaction_type)
    
    ggplot()+
      geom_bar(data=reci_with_mut_static,
               aes(x=interaction_type,y=n),
               stat = "identity",position = "stack",
               width=0.6,fill = "gray", color = "gray")+
      facet_grid(mut_cluster ~ ., scales = "free", space = "free")+
      theme_bw()+
      coord_flip() -> recieved_mut_statistic
    
    png(paste0(gsm_dir, "/recieved_mut_statistic.png"), width = 1225, height = 650)
    print(recieved_mut_statistic)
    dev.off()
    
    ##统计每cluster的获得丢失的LR集合
    #1.对数据框reci_with_mut进行操作，通过distinct()函数对mut_cluster、interaction_name和interaction_type这三列进行去重操作。
    #2.将interaction_name列按照"_"进行拆分，然后使用unnest()函数将拆分出的元素和其他列的元素构成新的行
    reci_with_mut_GO_input <- reci_with_mut %>%
      distinct(mut_cluster,interaction_name,interaction_type) %>%
      mutate(interaction_name = str_split(interaction_name, "_")) %>%
      unnest(cols = interaction_name)%>%
      distinct()
    
    
    enrich_result_all<-data.frame()
    for(mut_cluster_j in unique(reci_with_mut_GO_input$mut_cluster)){
      reci_with_mut_GO_input_1<-reci_with_mut_GO_input[reci_with_mut_GO_input$mut_cluster %in% mut_cluster_j,]
      for(type_i in unique(reci_with_mut_GO_input_1$interaction_type)){
        # 创建一个geneList，包含interaction_name列的元素
        geneList <- reci_with_mut_GO_input_1[reci_with_mut_GO_input_1$interaction_type %in% type_i,]$interaction_name
        # 进行GO功能富集分析
        enrich_result <- enrichGO(gene = geneList, 
                                  keyType = "SYMBOL", 
                                  OrgDb = org.Hs.eg.db, 
                                  ont = "BP", 
                                  pvalueCutoff = 0.05, 
                                  pAdjustMethod = "BH", 
                                  qvalueCutoff = 0.2)
        enrich_result@result$interaction_type<-type_i
        enrich_result@result$mut_cluster<-mut_cluster_j
        enrich_result_all<-rbind(enrich_result_all,enrich_result@result)
      }
    }
    ##挑选
    enrich_result_all_top3 <- enrich_result_all %>%
      group_by(mut_cluster,interaction_type) %>%
      arrange(desc(-log(qvalue))) %>%
      slice_head(n = 3)
    
    received_enrich_result_all_top3 <- enrich_result_all_top3
    
    ##功能富集
    ggplot(enrich_result_all_top3, aes(interaction_type, Description)) +
      geom_point(aes(color=-log(qvalue), size=Count))+theme_bw()+
      theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size = 12))+
      scale_color_gradient(low='#6699CC',high='#CC3333')+
      labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
      scale_y_discrete(position = "left")+
      facet_grid(mut_cluster ~ ., scales = "free", space = "free") -> recieved_GO_enrich
    
    png(paste0(gsm_dir, "/recieved_GO_enrich.png"), width = 1225, height = 650)
    print(recieved_GO_enrich)
    dev.off()
    
    # library(aplot)
    # p1%>%
    #   insert_left(mut_cluster, width = .05)
    #每组展示10个突变基因扰动强度最大的3对互作的河流图。按照mut_Gene和mut_cluster分组，然后按照互作强度改变绝对值从高到低排序，最后选择前3个互作
    reci_with_mut_top10 <- reci_with_mut %>% 
      group_by(mut_Gene,mut_cluster) %>%
      arrange(desc(abs(prob_diff))) %>%
      slice_head(n = 3)
    library(tidyr)
    # 将group1.source.target列按照" -> "拆分为source和target两列
    reci_with_mut_top10 <- separate(reci_with_mut_top10, group1.source.target, into = c("source", "target"), sep = " -> ",remove = FALSE)
    # 将source列中包含"Mut"的元素替换为mut_cluster的元素
    reci_with_mut_top10 <- reci_with_mut_top10 %>%
      mutate(source = ifelse(grepl("Mut", source), mut_cluster, source))
    # 将target列中包含"Mut"的元素替换为mut_cluster的元素
    reci_with_mut_top10 <- reci_with_mut_top10 %>%
      mutate(target = ifelse(grepl("Mut", target), mut_cluster, target))
    reci_with_mut_top10<-as.data.frame(reci_with_mut_top10)
    # 加载ggalluvial包
    library(ggalluvial)
    #用数据框reci_with_mut_top10画桑基图，第一列是mut_Gene，第二列是target，第三列是receptor，第五列是interaction_type，第六列是ligand，第七列是source
    # 绘制第一种桑基图
    ggplot(reci_with_mut_top10, aes(axis1 = mut_Gene, axis2 = target, axis3 = receptor, axis4 = interaction_type, axis5 = ligand, axis6 = source)) +
      geom_alluvium(aes(fill = target)) + # 线的颜色用source填充
      geom_stratum(aes(fill = target)) + # 第二列块的填充色用source填充
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_fill_brewer(palette = "Set3") +  # 使用调色板
      theme_void()+
      scale_x_continuous(breaks = 1:6, labels = c("mut_Gene","target","receptor","interaction_type","ligand","source"))+
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(size = 12, face = "bold",vjust = 5)) -> recieved_sankeyPlot
    
    png(paste0(gsm_dir, "/recieved_sankeyPlot.png"), width = 1225, height = 650)
    print(recieved_sankeyPlot)
    dev.off()
    
    
    #***********************************************数据库需要的数据（最后把所有样本的合到一起）******************************************************************
    #1.发送信号的-------------------------------------------------------------------------------------------------------
    #1.1 柱状图（每四个互作类型作为一个簇的group）
    #去掉other的，每个簇有四种互作类型
    #send_with_mut_static_noOther <- send_with_mut_static[-which(send_with_mut_static$interaction_type == "Other"),]
    #combinations <- expand.grid(unique(send_with_mut_static_noOther$mut_cluster), unique(send_with_mut_static_noOther$interaction_type))
    # send_with_mut_static_noOther_ordered <- send_with_mut_static_noOther[order(send_with_mut_static_noOther[,1]), ]
    # tmp_1 <- sapply(send_with_mut_static_noOther_ordered$interaction_type, function(x) paste0("'", x, "'"))
    # tmp_2 <- paste(tmp_1, collapse = ",")
    # tmp_3 <- paste0("[", tmp_2, "]")
    # 
    # tmp_4 <- paste(send_with_mut_static_noOther_ordered$n, collapse = ",")
    # tmp_5 <- paste0("[", tmp_4, "]")
    # 
    # tmp_6 <- sapply(send_with_mut_static_noOther_ordered$mut_cluster, function(x) paste0("'", x, "'"))
    # tmp_7 <- paste(tmp_6, collapse = ",")
    # tmp_8 <- paste0("[", tmp_7, "]")
    # 
    # send_static <- data.frame(
    #   interaction_type = tmp_3,
    #   interaction_number = tmp_5,
    #   mut_cluster = tmp_8,
    #   sample = sample_name,
    #   dataset = dataset_name,
    # )
    # 
    # write.table(send_static, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/send_statistic.txt"), sep = "\t", row.names = F, quote = F)
    # 
    #1.1 柱状图（x轴是簇，按照四种互作类型分组，图例也是四种互作类型）
    #去掉other的，每个簇有四种互作类型
    send_with_mut_static_noOther <- send_with_mut_static[-which(send_with_mut_static$interaction_type == "Other"),]
    send_with_mut_static_noOther_ordered <- send_with_mut_static_noOther[order(send_with_mut_static_noOther[,1]), ]
    tmp_1 <- sapply(unique(send_with_mut_static_noOther_ordered$mut_cluster), function(x) paste0("'", x, "'"))
    tmp_2 <- paste(tmp_1, collapse = ",")
    tmp_3 <- paste0("[", tmp_2, "]")
    
    tmp_6 <- sapply(unique(send_with_mut_static_noOther_ordered$interaction_type), function(x) paste0("'", x, "'"))
    tmp_7 <- paste(tmp_6, collapse = ",")
    
    intertype_num <- data.frame()
    #send_with_mut_static_noOther_ordered_1 <- send_with_mut_static_noOther_ordered[order(send_with_mut_static_noOther_ordered[,2]),]
    
    cluster_len <- length(unique(send_with_mut_static_noOther_ordered$mut_cluster))
    if(cluster_len*4 != nrow(send_with_mut_static_noOther_ordered)){
      combinations <- expand.grid(unique(send_with_mut_static_noOther_ordered$mut_cluster), unique(send_with_mut_static_noOther_ordered$interaction_type))
      names(combinations) <- c("mut_cluster", "interaction_type")
      send_with_mut_static_noOther_ordered_1 <- merge(combinations, send_with_mut_static_noOther_ordered, all.x = T)
      send_with_mut_static_noOther_ordered_1[is.na(send_with_mut_static_noOther_ordered_1)] <- 0
      send_with_mut_static_noOther_ordered_1 <- send_with_mut_static_noOther_ordered_1[order(send_with_mut_static_noOther_ordered_1[,2]),]
    }else{
      send_with_mut_static_noOther_ordered_1 <- send_with_mut_static_noOther_ordered[order(send_with_mut_static_noOther_ordered[,2]),]
    }
    
    for (kk in unique(send_with_mut_static_noOther_ordered_1$interaction_type)) {
      tmp_data <- send_with_mut_static_noOther_ordered_1[which(send_with_mut_static_noOther_ordered_1$interaction_type == kk),]
      #tmp_data_1 <- sapply(tmp_data$n, function(x) paste0("'", x, "'"))
      tmp_data_2 <- paste(tmp_data$n, collapse = ",")
      tmp_data_3 <- paste0("[", tmp_data_2, "]")
      
      intertype_num[kk,1] <- kk
      intertype_num[kk,2] <- tmp_data_3
    }
    
    #拼接成大的字符串
    interac_num <- paste(intertype_num$V2, collapse = ";")
    
    send_static_1 <- data.frame(
      Mut_cluster <- tmp_3,
      interact_type <- tmp_7,
      interact_num <- interac_num,
      sample = sample_name,
      dataset = dataset_name
    )
    names(send_static_1) <- c("Mut_cluster", "interaction_type", "interaction_number", "sample", "dataset")
    write.table(send_static_1, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/send_statistic.txt"), sep = "\t", row.names = F, quote = F)
    
    #1.2 气泡图
    enrich_result_all_top3_noOther <- enrich_result_all_top3[-which(enrich_result_all_top3$interaction_type == "Other"), c("Description", "qvalue", "Count", "interaction_type", "mut_cluster")]
    enrich_result_all_top3_noOther$log_q <- -log(enrich_result_all_top3_noOther$qvalue)
    
    send_dotplot <- data.frame()
    index <- 1
    for (kk in unique(enrich_result_all_top3_noOther$mut_cluster)) {
      tmp_1 <- enrich_result_all_top3_noOther[which(enrich_result_all_top3_noOther$mut_cluster == kk),]
      
      tmp_2 <- sapply(unique(enrich_result_all_top3_noOther$interaction_type), function(x) paste0("'", x, "'"))
      tmp_3 <- paste(tmp_2, collapse = ",")
      tmp_interaction_type <- paste0("[", tmp_3, "]")
      
      tmp_4 <- sapply(unique(tmp_1$Description), function(x) paste0("'", x, "'"))
      tmp_5 <- paste(tmp_4, collapse = ",")
      tmp_Description <- paste0("[", tmp_5, "]")
      
      combinations <- expand.grid(unique(tmp_1$Description), unique(enrich_result_all_top3_noOther$interaction_type))
      names(combinations) <- c("Description", "interaction_type")
      combinations$mut_cluster <- kk
      tmp_6 <- merge(combinations, tmp_1, all.x = T)
      tmp_6[is.na(tmp_6)] <- 0
      
      k <- 1
      value <- ""
      for (i in 1:length(unique(tmp_1$Description))) {
        for (j in 1:length(unique(tmp_1$interaction_type))) {
          
          tmp_7 <- paste0("{value:[", i-1, ",", j-1, ",", tmp_6[k,5], ",", tmp_6[k,6], "],symbolSize:", tmp_6[k,5], "},")
          k <- k+1
          value <- paste0(value, tmp_7)
          
        }
      }
      
      send_dotplot[index,1] <- kk
      send_dotplot[index,2] <- tmp_interaction_type
      send_dotplot[index,3] <- tmp_Description
      send_dotplot[index,4] <- value
      
      index <- index+1
    }  
    
    names(send_dotplot) <- c("Mut_cluster", "interaction_type", "Description", "enrich_value")
    send_dotplot$sample <- sample_name
    send_dotplot$dataset <- dataset_name
    
    write.table(send_dotplot, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/send_dotplot.txt"), quote = F, row.names = F, sep = "\t")
    
    #1.3 桑基图
    #导出json的数据格式
    tmp_1 <- sapply(unique(send_with_mut_top10$mut_Gene), function(x) paste0("{'name':'", x, "'}"))
    tmp_2 <- paste(tmp_1, collapse = ",")
    
    send_with_mut_top10$source_1 <- paste('source_', send_with_mut_top10$source, sep = '')
    tmp_11 <- sapply(unique(send_with_mut_top10$source_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_22 <- paste(tmp_11, collapse = ",")
    
    send_with_mut_top10$ligand_1 <- paste('ligand_', send_with_mut_top10$ligand, sep = '')
    tmp_111 <- sapply(unique(send_with_mut_top10$ligand_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_222 <- paste(tmp_111, collapse = ",")
    
    tmp_1111 <- sapply(unique(send_with_mut_top10$interaction_type), function(x) paste0("{'name':'", x, "'}"))
    tmp_2222 <- paste(tmp_1111, collapse = ",")
    
    send_with_mut_top10$receptor_1 <- paste('receptor_', send_with_mut_top10$receptor, sep = '')
    tmp_11111 <- sapply(unique(send_with_mut_top10$receptor_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_22222 <- paste(tmp_11111, collapse = ",")
    
    send_with_mut_top10$target_1 <- paste('target_', send_with_mut_top10$target, sep = '')
    tmp_111111 <- sapply(unique(send_with_mut_top10$target_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_222222 <- paste(tmp_111111, collapse = ",")
    
    nodes <- paste0(tmp_2, ",", tmp_22, ",", tmp_222, ",", tmp_2222, ",", tmp_22222, ",", tmp_222222)
    
    links <- ""
    #①mut_Gene->source的值
    mutG_source <- send_with_mut_top10[,c("mut_Gene", "source")]
    mutG_source$source_1 <- paste('source_', mutG_source$source, sep = '')
    mutG_source_1 <- mutG_source %>%
      group_by(mut_Gene, source_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(mutG_source_1)) {
      tmp_3 <- paste0("{'source': '", mutG_source_1[i,1], "',", "'target':'",mutG_source_1[i,2], "',", "'value':", mutG_source_1[i,3], "},")
      links <- paste0(links, tmp_3)
    }
    #②source->ligand的值
    source_ligand <- send_with_mut_top10[,c("source", "ligand")]
    source_ligand$source_1 <- paste('source_', source_ligand$source, sep = '')
    source_ligand$ligand_1 <- paste('ligand_', source_ligand$ligand, sep = '')
    source_ligand_1 <- source_ligand %>%
      group_by(source_1, ligand_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(source_ligand_1)) {
      tmp_4 <- paste0("{'source': '", source_ligand_1[i,1], "',", "'target':'",source_ligand_1[i,2], "',", "'value':", source_ligand_1[i,3], "},")
      links <- paste0(links, tmp_4)
    }
    #③ligand->interaction_type的值
    ligand_intertype <- send_with_mut_top10[,c("ligand", "interaction_type")]
    ligand_intertype$ligand_1 <- paste('ligand_', ligand_intertype$ligand, sep = '')
    ligand_intertype_1 <- ligand_intertype %>%
      group_by(ligand_1, interaction_type) %>%
      summarise(Count = n())
    for (i in 1:nrow(ligand_intertype_1)) {
      tmp_5 <- paste0("{'source': '", ligand_intertype_1[i,1], "',", "'target':'",ligand_intertype_1[i,2], "',", "'value':", ligand_intertype_1[i,3], "},")
      links <- paste0(links, tmp_5)
    }
    #④interaction_type->receptor的值
    intertype_receptor <- send_with_mut_top10[,c("interaction_type", "receptor")]
    intertype_receptor$receptor_1 <- paste('receptor_', intertype_receptor$receptor, sep = '')
    intertype_receptor_1 <- intertype_receptor %>%
      group_by(interaction_type, receptor_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(intertype_receptor_1)) {
      tmp_6 <- paste0("{'source': '", intertype_receptor_1[i,1], "',", "'target':'",intertype_receptor_1[i,2], "',", "'value':", intertype_receptor_1[i,3], "},")
      links <- paste0(links, tmp_6)
    }
    #⑤receptor->target的值
    receptor_target <- send_with_mut_top10[,c("receptor", "target")]
    receptor_target$receptor_1 <- paste('receptor_', receptor_target$receptor, sep = '')
    receptor_target$target_1 <- paste('target_', receptor_target$target, sep = '')
    receptor_target_1 <- receptor_target %>%
      group_by(receptor_1, target_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(receptor_target_1)) {
      tmp_7 <- paste0("{'source': '", receptor_target_1[i,1], "',", "'target':'",receptor_target_1[i,2], "',", "'value':", receptor_target_1[i,3], "},")
      links <- paste0(links, tmp_7)
    }
    
    sankey_plot <- data.frame(nodes = nodes,links = links)
    sankey_plot$sample <- sample_name
    sankey_plot$dataset <- dataset_name
    
    sankey_plot$nodes <- gsub("'", '"', sankey_plot$nodes)
    sankey_plot$links <- gsub("'", '"', sankey_plot$links)
    
    write.table(sankey_plot, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/send_sankey.txt"), sep = "\t", quote = F, row.names = F)
    
    
    #2.接收信号的--------------------------------------------------------------------------------------------------
    #2.1 柱状图（x轴是簇，按照四种互作类型分组，图例也是四种互作类型）
    #去掉other的，每个簇有四种互作类型
    reci_with_mut_static_noOther <- reci_with_mut_static[-which(reci_with_mut_static$interaction_type == "Other"),]
    reci_with_mut_static_noOther_ordered <- reci_with_mut_static_noOther[order(reci_with_mut_static_noOther[,1]), ]
    tmp_1 <- sapply(unique(reci_with_mut_static_noOther_ordered$mut_cluster), function(x) paste0("'", x, "'"))
    tmp_2 <- paste(tmp_1, collapse = ",")
    tmp_3 <- paste0("[", tmp_2, "]")
    
    tmp_6 <- sapply(unique(reci_with_mut_static_noOther_ordered$interaction_type), function(x) paste0("'", x, "'"))
    tmp_7 <- paste(tmp_6, collapse = ",")
    
    intertype_num <- data.frame()
    
    cluster_len <- length(unique(reci_with_mut_static_noOther_ordered$mut_cluster))
    if(cluster_len*4 != nrow(reci_with_mut_static_noOther_ordered)){
      combinations <- expand.grid(unique(reci_with_mut_static_noOther_ordered$mut_cluster), unique(reci_with_mut_static_noOther_ordered$interaction_type))
      names(combinations) <- c("mut_cluster", "interaction_type")
      reci_with_mut_static_noOther_ordered_1 <- merge(combinations, reci_with_mut_static_noOther_ordered, all.x = T)
      reci_with_mut_static_noOther_ordered_1[is.na(reci_with_mut_static_noOther_ordered_1)] <- 0
      reci_with_mut_static_noOther_ordered_1 <- reci_with_mut_static_noOther_ordered_1[order(reci_with_mut_static_noOther_ordered_1[,2]),]
    }else{
      reci_with_mut_static_noOther_ordered_1 <- reci_with_mut_static_noOther_ordered[order(reci_with_mut_static_noOther_ordered[,2]),]
    }
    
    
    for (kk in unique(reci_with_mut_static_noOther_ordered_1$interaction_type)) {
      tmp_data <- reci_with_mut_static_noOther_ordered_1[which(reci_with_mut_static_noOther_ordered_1$interaction_type == kk),]
      tmp_data_2 <- paste(tmp_data$n, collapse = ",")
      tmp_data_3 <- paste0("[", tmp_data_2, "]")
      
      intertype_num[kk,1] <- kk
      intertype_num[kk,2] <- tmp_data_3
    }
    
    #拼接成大的字符串
    interac_num <- paste(intertype_num$V2, collapse = ";")
    
    reci_static_1 <- data.frame(
      Mut_cluster <- tmp_3,
      interact_type <- tmp_7,
      interact_num <- interac_num,
      sample = sample_name,
      dataset = dataset_name
    )
    names(reci_static_1) <- c("Mut_cluster", "interaction_type", "interaction_number", "sample", "dataset")
    write.table(reci_static_1, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/recived_statistic.txt"), sep = "\t", row.names = F, quote = F)
    
    #2.2 气泡图
    enrich_result_all_top3_noOther <- received_enrich_result_all_top3[-which(received_enrich_result_all_top3$interaction_type == "Other"), c("Description", "qvalue", "Count", "interaction_type", "mut_cluster")]
    enrich_result_all_top3_noOther$log_q <- -log(enrich_result_all_top3_noOther$qvalue)
    
    recived_dotplot <- data.frame()
    index <- 1
    for (kk in unique(enrich_result_all_top3_noOther$mut_cluster)) {
      tmp_1 <- enrich_result_all_top3_noOther[which(enrich_result_all_top3_noOther$mut_cluster == kk),]
      
      tmp_2 <- sapply(unique(enrich_result_all_top3_noOther$interaction_type), function(x) paste0("'", x, "'"))
      tmp_3 <- paste(tmp_2, collapse = ",")
      tmp_interaction_type <- paste0("[", tmp_3, "]")
      
      tmp_4 <- sapply(unique(tmp_1$Description), function(x) paste0("'", x, "'"))
      tmp_5 <- paste(tmp_4, collapse = ",")
      tmp_Description <- paste0("[", tmp_5, "]")
      
      combinations <- expand.grid(unique(tmp_1$Description), unique(enrich_result_all_top3_noOther$interaction_type))
      names(combinations) <- c("Description", "interaction_type")
      combinations$mut_cluster <- kk
      tmp_6 <- merge(combinations, tmp_1, all.x = T)
      tmp_6[is.na(tmp_6)] <- 0
      
      k <- 1
      value <- ""
      for (i in 1:length(unique(tmp_1$Description))) {
        for (j in 1:length(unique(tmp_1$interaction_type))) {
          
          tmp_7 <- paste0("{value:[", i-1, ",", j-1, ",", tmp_6[k,5], ",", tmp_6[k,6], "],symbolSize:", tmp_6[k,5], "},")
          k <- k+1
          value <- paste0(value, tmp_7)
          
        }
      }
      
      recived_dotplot[index,1] <- kk
      recived_dotplot[index,2] <- tmp_interaction_type
      recived_dotplot[index,3] <- tmp_Description
      recived_dotplot[index,4] <- value
      
      index <- index+1
    }  
    
    names(recived_dotplot) <- c("Mut_cluster", "interaction_type", "Description", "enrich_value")
    recived_dotplot$sample <- sample_name
    recived_dotplot$dataset <- dataset_name
    
    write.table(recived_dotplot, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/recived_dotplot.txt"), quote = F, row.names = F, sep = "\t")
    
    #2.3 桑基图
    #导出json的数据格式
    tmp_1 <- sapply(unique(reci_with_mut_top10$mut_Gene), function(x) paste0("{'name':'", x, "'}"))
    tmp_2 <- paste(tmp_1, collapse = ",")
    
    reci_with_mut_top10$source_1 <- paste('source_', reci_with_mut_top10$source, sep = '')
    tmp_11 <- sapply(unique(reci_with_mut_top10$source_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_22 <- paste(tmp_11, collapse = ",")
    
    reci_with_mut_top10$ligand_1 <- paste('ligand_', reci_with_mut_top10$ligand, sep = '')
    tmp_111 <- sapply(unique(reci_with_mut_top10$ligand_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_222 <- paste(tmp_111, collapse = ",")
    
    tmp_1111 <- sapply(unique(reci_with_mut_top10$interaction_type), function(x) paste0("{'name':'", x, "'}"))
    tmp_2222 <- paste(tmp_1111, collapse = ",")
    
    reci_with_mut_top10$receptor_1 <- paste('receptor_', reci_with_mut_top10$receptor, sep = '')
    tmp_11111 <- sapply(unique(reci_with_mut_top10$receptor_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_22222 <- paste(tmp_11111, collapse = ",")
    
    reci_with_mut_top10$target_1 <- paste('target_', reci_with_mut_top10$target, sep = '')
    tmp_111111 <- sapply(unique(reci_with_mut_top10$target_1), function(x) paste0("{'name':'", x, "'}"))
    tmp_222222 <- paste(tmp_111111, collapse = ",")
    
    nodes <- paste0(tmp_2, ",", tmp_22, ",", tmp_222, ",", tmp_2222, ",", tmp_22222, ",", tmp_222222)
    
    links <- ""
    #①mut_Gene->source的值
    mutG_source <- reci_with_mut_top10[,c("mut_Gene", "source")]
    mutG_source$source_1 <- paste('source_', mutG_source$source, sep = '')
    mutG_source_1 <- mutG_source %>%
      group_by(mut_Gene, source_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(mutG_source_1)) {
      tmp_3 <- paste0("{'source': '", mutG_source_1[i,1], "',", "'target':'",mutG_source_1[i,2], "',", "'value':", mutG_source_1[i,3], "},")
      links <- paste0(links, tmp_3)
    }
    #②source->ligand的值
    source_ligand <- reci_with_mut_top10[,c("source", "ligand")]
    source_ligand$source_1 <- paste('source_', source_ligand$source, sep = '')
    source_ligand$ligand_1 <- paste('ligand_', source_ligand$ligand, sep = '')
    source_ligand_1 <- source_ligand %>%
      group_by(source_1, ligand_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(source_ligand_1)) {
      tmp_4 <- paste0("{'source': '", source_ligand_1[i,1], "',", "'target':'",source_ligand_1[i,2], "',", "'value':", source_ligand_1[i,3], "},")
      links <- paste0(links, tmp_4)
    }
    #③ligand->interaction_type的值
    ligand_intertype <- reci_with_mut_top10[,c("ligand", "interaction_type")]
    ligand_intertype$ligand_1 <- paste('ligand_', ligand_intertype$ligand, sep = '')
    ligand_intertype_1 <- ligand_intertype %>%
      group_by(ligand_1, interaction_type) %>%
      summarise(Count = n())
    for (i in 1:nrow(ligand_intertype_1)) {
      tmp_5 <- paste0("{'source': '", ligand_intertype_1[i,1], "',", "'target':'",ligand_intertype_1[i,2], "',", "'value':", ligand_intertype_1[i,3], "},")
      links <- paste0(links, tmp_5)
    }
    #④interaction_type->receptor的值
    intertype_receptor <- reci_with_mut_top10[,c("interaction_type", "receptor")]
    intertype_receptor$receptor_1 <- paste('receptor_', intertype_receptor$receptor, sep = '')
    intertype_receptor_1 <- intertype_receptor %>%
      group_by(interaction_type, receptor_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(intertype_receptor_1)) {
      tmp_6 <- paste0("{'source': '", intertype_receptor_1[i,1], "',", "'target':'",intertype_receptor_1[i,2], "',", "'value':", intertype_receptor_1[i,3], "},")
      links <- paste0(links, tmp_6)
    }
    #⑤receptor->target的值
    receptor_target <- reci_with_mut_top10[,c("receptor", "target")]
    receptor_target$receptor_1 <- paste('receptor_', receptor_target$receptor, sep = '')
    receptor_target$target_1 <- paste('target_', receptor_target$target, sep = '')
    receptor_target_1 <- receptor_target %>%
      group_by(receptor_1, target_1) %>%
      summarise(Count = n())
    for (i in 1:nrow(receptor_target_1)) {
      tmp_7 <- paste0("{'source': '", receptor_target_1[i,1], "',", "'target':'",receptor_target_1[i,2], "',", "'value':", receptor_target_1[i,3], "},")
      links <- paste0(links, tmp_7)
    }
    
    sankey_plot <- data.frame(nodes = nodes,links = links)
    sankey_plot$sample <- sample_name
    sankey_plot$dataset <- dataset_name
    
    sankey_plot$nodes <- gsub("'", '"', sankey_plot$nodes)
    sankey_plot$links <- gsub("'", '"', sankey_plot$links)
    
    write.table(sankey_plot, paste0("E:/数据库/STMut/数据处理/结果/result_table/", dataset_name, "/", sample_name, "/recived_sankey.txt"), sep = "\t", quote = F, row.names = F)
    
    
  }
}

