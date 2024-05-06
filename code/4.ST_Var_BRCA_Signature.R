######################################构建spot突变簇类别somatic突变
##参考https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247494325&idx=1&sn=94891e349ceb9b459e91e86db144b4ea&chksm=9b4baa0eac3c2318bfb4cfa66875078db89baf0eaf4543ed785ed4892532cbdd0145fb2ba1fd#rd
file_path <- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/Step4_VariantCalling"
##0-读取scsnv注释矩阵
anno<- read.table(paste0(file_path,"/GSM6177599.variants.annovar.Step4.3.1.hg38_multianno.csv"), header=T,sep=",",as.is=T,fill=T)
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
  colnames(df2)[1]<-"Chr"
  df2<-merge(df2,anno2)
  gene_mut_spot<-rbind(gene_mut_spot,df2)
}
gene_mut_spot$CB<-paste0(gene_mut_spot$CB,"-1")
##2-识别的spot突变簇类别
data_dir<- "E:/数据库/STMut/数据处理/BRCA/GSE203612/GSM6177599/"
spot_cluster<- read.table(paste(data_dir,tail(unlist(strsplit(data_dir,"/")),1),"_spot.mut_cluster.txt",sep=""), header=T,sep="\t",as.is=T,fill=T,row.names = 1)
spot_cluster$mut_cluster<-paste0("C_",spot_cluster$mut_cluster)
Maf_data<-merge(gene_mut_spot,spot_cluster)
######################################提取突变上下文已经计算96突变形式的比例
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg38)
Maf_data<-Maf_data[which(Maf_data$REF !=Maf_data$ALT_expected),]
sca_vr = VRanges(
  seqnames =  Maf_data$Chr ,
  ranges = IRanges(start = Maf_data$Start,end = Maf_data$Start+1),
  ref = Maf_data$REF,
  alt = Maf_data$ALT_expected,
  sampleNames = as.character(Maf_data$mut_cluster),
  study=as.character(Maf_data$mut_cluster))
#从基因组根据坐标获取碱基上下文
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg38)
# 对每个cluster，计算 96 突变可能性的 比例分布情况
escc_sca_mm = motifMatrix(sca_motifs, group = "study", normalize = F)

######################################空间分布簇96个三核苷酸突变图
#  install the latest version from GitHub
if(!require(devtools)) install.packages("devtools")
devtools::install_github("Niinleslie/MesKit")
library("MesKit")
library(ComplexHeatmap)
library(viridis)
trimatrix_sample<-list()
trimatrix_sample[["P"]][["tri_matrix"]]<-data.frame(t(escc_sca_mm))
colnames(trimatrix_sample[["P"]][["tri_matrix"]])<-c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
trimatrix_sample[["P"]][["tsb.label"]]<-data.frame(Tumor_Sample_Barcode=unique(Maf_data$mut_cluster),Tumor_Sample_Label=unique(Maf_data$mut_cluster))

len <- length(plotMutSigProfile(trimatrix_sample))
for (kk in 1:len) {
  fil_name <- paste0("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/mut_signature_", kk, ".png")
  p1 <- plotMutSigProfile(trimatrix_sample)[[kk]]
  png(fil_name, width = 1225, height = 650)
  print(p1)
  dev.off()
}
plotMutSigProfile(trimatrix_sample)[[1]]
plotMutSigProfile(trimatrix_sample)[[2]]
plotMutSigProfile(trimatrix_sample)[[3]]
plotMutSigProfile(trimatrix_sample)[[4]]
plotMutSigProfile(trimatrix_sample)[[5]]
plotMutSigProfile(trimatrix_sample)[[6]]
######################################空间分布簇突变集与COSMIC特征相关性######################################
##每簇突变集与COSMIC特征相关性
fit_sample <- fitSignatures(trimatrix_sample,signaturesRef="exome_cosmic_v3")
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/Cosine_similarity.png", width = 1225, height = 522)
ComplexHeatmap::Heatmap(fit_sample$P$cosine.similarity, name = "Cosine similarity",col = viridis(100))
dev.off()

fit_sample_cos_similarity <- as.data.frame(fit_sample$P$cosine.similarity)
fit_sample_cos_similarity$cluster <- rownames(fit_sample_cos_similarity)
library(tidyverse)
# 使用gather函数将宽数据转换为长数据
long_data <- fit_sample_cos_similarity %>% 
  gather(key = "COSMIC Signatures", value = "Cosine Similarity", -cluster)
SBS_anno <- read.table("E:/数据库/STMut/数据处理/结果表格/SBS_anno.txt", header = T, sep = "\t", check.names = F)
Cosmic_signature_anno <- merge(long_data, SBS_anno)

write.table(Cosmic_signature_anno, "E:/数据库/STMut/数据处理/结果表格/Cosmic_signature_anno.txt", sep = "\t", quote = F, row.names = F)
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
png("E:/数据库/STMut/数据处理/结果/BRCA/GSE203612/GSM6177599/mut_sigature_contribution.png", width = 1225, height = 650)
ggplot(data_long, aes(x = Category, y = Value, fill = Signature)) +
  geom_bar(stat = "identity",color = "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "Mutational signature contribution",  y = "Mutational signature
       contribution") +
  theme_classic()+
  coord_flip()  # 翻转坐标轴
dev.off()
