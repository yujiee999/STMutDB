filefolder_num <- 0
file_num <- 0
input_path <- "E:/数据库/STMut/数据处理/结果/result"
# <- "H:/database/STMut/data_process/result"
# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = input_path, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (i in 1:40) {
  gse_dir <- gse_dirs[i]
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  filefolder_num = filefolder_num + length(gsm_dirs)
  
  for(gsm_dir in gsm_dirs) {
    # 列出当前路径下的所有文件和文件夹
    files_and_folders <- list.files(gsm_dir)
    # 过滤出文件
    files <- files_and_folders[!file.info(files_and_folders)$isdir]
    # 计算文件的数量
    num_files <- length(files)
    
    #print(num_files)#每个样本应该有15个文件
    if(num_files != 15){
      print(gsm_dir)
    }
  }
  
}
print(filefolder_num)#一共338个样本





