# 设置源文件夹和目标文件夹的路径
source_dir <- "E:/数据库/STMut/数据处理/datasets_series_matrix"
destination_dir <- "E:/数据库/STMut/数据处理/datasets_series_matrix2"

# 定义函数来解压指定路径下的所有.txt.gz文件到另一个指定路径
unzip_all_txt_gz_files <- function(source_dir, destination_dir) {
  # 列出source_dir下所有以.txt.gz结尾的文件
  files <- list.files(source_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  
  # 遍历每个文件，使用gunzip函数解压文件，并将解压后的文件名去掉.gz后缀
  for (file in files) {
    # 读取.gz压缩文件
    con <- gzfile(file)
    text <- readLines(con)
    close(con)
    
    # 确定目的地文件路径，即去掉原压缩文件的路径和.gz后缀的文件名
    destination_file <- gsub("\\.gz$", "", basename(file))
    destination_path <- file.path(destination_dir, destination_file)
    
    # 将解压后的内容写入目的地文件
    writeLines(text, destination_path)
  }
}

# 调用函数，解压指定路径下的所有.txt.gz文件到目标路径
unzip_all_txt_gz_files(source_dir, destination_dir)


#------------------------------------------
rm(list = ls())
input_path <- "E:/数据库/STMut/数据处理/datasets_series_matrix2"
datasets_samples_clinical_info <- data.frame()
all_files <- list.files(destination_dir)
for (single_file in all_files) {
  #获取数据集的名字
  dataset_name <- unlist(str_split(single_file, "_"))[1]
  # 逐行读取文本文件
  lines <- readLines(paste0(input_path, "/", single_file))
  # 提取以!Sample_geo_accession开头的行
  Sample_geo_accession <- lines[grep("^!Sample_geo_accession", lines)]
  # 提取以!Sample_source_name_ch1开头的行
  Sample_source_name_ch1 <- lines[grep("^!Sample_source_name_ch1", lines)]
  # 提取以!Sample_organism_ch1开头的行
  Sample_organism_ch1 <- lines[grep("^!Sample_organism_ch1", lines)]
  # 提取以!Sample_characteristics_ch1开头的行
  Sample_characteristics_ch1 <- lines[grep("^!Sample_characteristics_ch1", lines)]
  
  tmp <- as.data.frame(Sample_characteristics_ch1[1])
  tmp_Sample_characteristics <- tmp
  len <- length(Sample_characteristics_ch1)
  if(len == 1){
    
    tmp_Sample_characteristics$Sample_characteristics_ch1.2 <- NA
    tmp_Sample_characteristics$Sample_characteristics_ch1.3 <- NA
    tmp_Sample_characteristics$Sample_characteristics_ch1.4 <- NA
    tmp_Sample_characteristics$Sample_characteristics_ch1.5 <- NA
    tmp_Sample_characteristics$Sample_characteristics_ch1.6 <- NA
    tmp_Sample_characteristics$Sample_characteristics_ch1.7 <- NA
    
  }else{
    
    for (ii in 2:len) {
      tmp <- as.data.frame(Sample_characteristics_ch1[ii])
      tmp_Sample_characteristics <- cbind(tmp_Sample_characteristics, tmp)
    }
    
    if(len == 2){
      tmp_Sample_characteristics$Sample_characteristics_ch1.3 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.4 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.5 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.6 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.7 <- NA
    }
    if(len == 3){
      tmp_Sample_characteristics$Sample_characteristics_ch1.4 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.5 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.6 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.7 <- NA
    }
    if(len == 4){
      tmp_Sample_characteristics$Sample_characteristics_ch1.5 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.6 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.7 <- NA
    }
    if(len == 5){
      tmp_Sample_characteristics$Sample_characteristics_ch1.6 <- NA
      tmp_Sample_characteristics$Sample_characteristics_ch1.7 <- NA
    }
    if(len == 6){
      tmp_Sample_characteristics$Sample_characteristics_ch1.7 <- NA
    }
  }
  
  colnames(tmp_Sample_characteristics) <- c("Sample_characteristics_ch1.1", "Sample_characteristics_ch1.2", "Sample_characteristics_ch1.3", "Sample_characteristics_ch1.4", "Sample_characteristics_ch1.5", "Sample_characteristics_ch1.6", "Sample_characteristics_ch1.7")
  
  tmp_data <- data.frame(dataset = dataset_name, 
                         Sample_geo_accession = Sample_geo_accession,
                         Sample_source_name_ch1 = Sample_source_name_ch1,
                         Sample_organism_ch1 = Sample_organism_ch1
                         #Sample_characteristics_ch1 = Sample_characteristics_ch1,
                        )
  
  tmp_data_1 <- cbind(tmp_data, tmp_Sample_characteristics)
  datasets_samples_clinical_info <- rbind(datasets_samples_clinical_info, tmp_data_1)
}

#--------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(reshape2)

datasets_samples_clinical_info_1 <- datasets_samples_clinical_info %>%
                                    gather(key, value, -dataset) %>%
                                    separate_rows(value, sep = "\t") %>%
                                    as.data.frame()
                                    #spread(key, value) %>%
                                    #select(dataset, everything())

for (i in unique(datasets_samples_clinical_info_1$dataset)) {
    tmp_dataset <- datasets_samples_clinical_info_1[which(datasets_samples_clinical_info_1$dataset == i),]
    # 将数据转换为长格式，group列转为行标签
    #tmp_dataset_1 <- dcast(tmp_dataset, dataset ~ key, value.var = "value")
    tmp_dataset_1 <- tmp_dataset %>%
      pivot_wider(names_from = key, values_from = value)
    
}






library( "GEOquery" )
## 取表达矩阵和样本信息表
result <- data.frame()
myfile <- list.files("E:/数据库/STMut/数据处理/datasets_series_matrix")
for (i in myfile) {
  
  GSE_data <- getGEO(filename = paste0("E:/数据库/STMut/数据处理/datasets_series_matrix/", i), getGPL = F)
  pdata = pData( GSE_data )
  tmp_dataset <- unlist(str_split(i, "_"))[1]
  
  write.table(pdata, paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", tmp_dataset, ".txt"), sep="\t", quote = F, row.names = F)
  
}

myfile_1 <- list.files("E:/数据库/STMut/数据处理/datasets_series_matrix3")
for (j in myfile_1) {
  
  tmp_dataset <- unlist(str_split(j, ".txt"))[1]
  df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
  df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1","characteristics_ch1", "characteristics_ch1.1")]
  df_1$dataset <- tmp_dataset
  
  result <- rbind(result, df_1)
}

#GSE171351
j <- myfile_1[10]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1")]
df_1$characteristics_ch1.1 <- NA
df_1$dataset <- tmp_dataset
result <- rbind(result, df_1)

j <- myfile_1[11]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
result <- rbind(result, df_1)

j <- myfile_1[12]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
result <- rbind(result, df_1)

j <- myfile_1[13]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
result <- rbind(result, df_1)

result$characteristics_ch1.2 <- NA
result$characteristics_ch1.3 <- NA

j <- myfile_1[14]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
colnames(df_1)[6] <- "characteristics_ch1.2"
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[15]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
result$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[16]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[17]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[18]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[19]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[20]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[21]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[22]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[23]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[24]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
colnames(df_1)[6] <- "characteristics_ch1.2"
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[25]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
colnames(df_1)[6] <- "characteristics_ch1.2"
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[26]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[27]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[28]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[29]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.3", "characteristics_ch1.4")]
df_1$dataset <- tmp_dataset
colnames(df_1)[6] <- "characteristics_ch1.2"
colnames(df_1)[7] <- "characteristics_ch1.3"
result <- rbind(result, df_1)

j <- myfile_1[30]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[31]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[32]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[33]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[34]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
colnames(df_1)[6] <- "characteristics_ch1.2"
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[35]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[36]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.1 <- NA
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[37]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[38]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[39]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.3", "characteristics_ch1.4")]
df_1$dataset <- tmp_dataset
colnames(df_1)[5] <- "characteristics_ch1.1"
colnames(df_1)[6] <- "characteristics_ch1.2"
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[40]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.1 <- NA
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[41]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.1 <- NA
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[42]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[43]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
result <- rbind(result, df_1)

j <- myfile_1[44]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3")]
df_1$dataset <- tmp_dataset
result <- rbind(result, df_1)

j <- myfile_1[45]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[46]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.1 <- NA
df_1$characteristics_ch1.2 <- NA
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[47]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

j <- myfile_1[48]
tmp_dataset <- unlist(str_split(j, ".txt"))[1]
df <- read.table(paste0("E:/数据库/STMut/数据处理/datasets_series_matrix3/", j), sep = "\t", header = T, fill = T)
df_1 <- df[, c("geo_accession", "source_name_ch1", "organism_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
df_1$dataset <- tmp_dataset
df_1$characteristics_ch1.3 <- NA
result <- rbind(result, df_1)

library(readxl)
dataset_sample_tissue <- read_excel("E:/数据库/STMut/data/dataset_sample_tissue.xlsx")
sample <- unique(dataset_sample_tissue$GSM)

result_1 <- result %>% as_tibble() %>% 
  separate_rows(dataset, sep = "-")
#
result_2 <- unique(result_1[result_1$dataset %like% "^GSE", ])

result_3 <- result[which(result_2$geo_accession %in% sample),]

write.table(result_3, "E:/数据库/STMut/数据处理/结果表格/sample_basic.txt", sep = "\t", quote = F, row.names = F)


aa <- unique(result$geo_accession)
length(unique(result_2$dataset))#40
length(unique(dataset_sample_tissue$GSE))#40

length(unique(dataset_sample_tissue$GSM))#363
length(unique(result_2$geo_accession))#460
length(intersect(unique(dataset_sample_tissue$GSM), unique(result_2$geo_accession)))#346











                                    
