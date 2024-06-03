# 指定待搜索的路径
path_to_search <- "E:/数据库/STMut/数据处理/结果"

# 指定目标路径
destination_path <- "E:/数据库/STMut/数据处理/结果/result"

# 在给定目录中找到所有GSE开头的文件夹
gse_dirs <- list.dirs(path = path_to_search, full.names = TRUE, recursive = FALSE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

for (gse_dir in gse_dirs) {
  # 在 GSE 文件夹下找到所有 GSM 开头的文件夹
  gsm_dirs <- list.dirs(path = gse_dir, full.names = TRUE)
  gsm_dirs <- gsm_dirs[grepl("GSM", basename(gsm_dirs))]
  
  if(gse_dir == "E:/数据库/STMut/数据处理/结果/GSE144239_ST"){
    
    for(gsm_dir in gsm_dirs) {
      # 在每个 GSM 文件夹下找到以.jpg 结尾的文件并复制到目标路径
      jpg_files <- list.files(path = gsm_dir, pattern = "*.jpg", full.names = TRUE)
      destination_file_dir = gsub(path_to_search, destination_path, gsm_dir)
      destination_file_path <- file.path(destination_file_dir, "tissue_lowres_image.png")
      file.copy(jpg_files, destination_file_path)
    }
    
  }else{
    
    for(gsm_dir in gsm_dirs) {
      # 在每个 GSM 文件夹下找到 spatial 文件夹
      spatial_dir <- list.dirs(path = gsm_dir, full.names = TRUE)
      spatial_dir <- spatial_dir[grepl("spatial", basename(spatial_dir))]
      
      # 找到 tissue_lowres_image.png 文件并复制到目标路径
      png_file <- list.files(path = spatial_dir, pattern = "tissue_lowres_image.png", full.names = TRUE)
      
      if (length(png_file) > 0) {
        
        destination_file_dir = gsub(path_to_search, destination_path, gsm_dir)
        destination_file_path <- file.path(destination_file_dir, "tissue_lowres_image.png")
        file.copy(png_file, destination_file_path)
        
      }
      
    }
    
  }
  
}




