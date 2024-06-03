# 加载readxl包
library(readxl)
library(dplyr)
# 读取.xlsx文件
data <- read_excel("E:/数据库/STMut/data/dataset_sample_tissue.xlsx")
data[is.na(data)] <- "-"
data$Disease <- gsub("psoriasis", "Psoriasis", data$Disease)
write.table(data, "E:/数据库/STMut/数据处理/结果表格/sample_basic.txt", sep = "\t", quote = F, row.names = F)

data1 <- data %>%
		group_by(Species, PMID, Disease, Dataset, Technology, Tissue, Abbreviation) %>%
		summarize(Sample_count = n_distinct(Sample))
data1[is.na(data1)] <- "-"
write.table(data1, "E:/数据库/STMut/数据处理/结果表格/dataset_basic.txt", sep = "\t", quote = F, row.names = F)


data_disease <- data %>%
		group_by(Disease) %>%
		summarize(Sample_count = n_distinct(Sample))
data_disease <- data_disease[-which(data_disease$Disease=="Normal"),]
colnames(data_disease)[1] <- "Tissue"
data_disease <- as.data.frame(data_disease)
data_disease$type <- "Disease"

data_normal <- data[which(data$Disease=="Normal"),]
data_normal <- data_normal %>%
			group_by(Tissue) %>%
		summarize(Sample_count = n_distinct(Sample))
data_normal <- as.data.frame(data_normal)
data_normal$type <- "Normal"

data4 <- NULL
data4 <- rbind(data_disease, data_normal)
write.table(data4, "E:/数据库/STMut/数据处理/结果表格/tissue_count.txt", sep = "\t", quote = F, row.names = F)
