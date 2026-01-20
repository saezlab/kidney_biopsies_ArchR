library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)
library(hdf5r)

# Helper functions
ArchR_write_h5 <- function(matrix = NULL, file = NULL, matrix.type = c("PeakMatrix", "GeneScoreMatrix")){
  if(is.null(file)){
    stop('No such file or directory')
  }
  if(matrix.type != "PeakMatrix"){
    stop("Currently only conversion from PeakMatrix is supported.")
  }
  if(class(matrix) != "RangedSummarizedExperiment"){
    stop("Class of", substitute(matrix), "is not RangedSummarizedExperiment")
  }
  h5 <- hdf5r::H5File$new(filename = file, mode = 'w')
  tryCatch({
    ArchR_to_h5(mat = matrix, h5 = h5, matrix.type = matrix.type)
    hdf5r::h5attr(h5, 'matrix_type') <- matrix.type
  },
  error = function(e) print(e),
  finally = {
    h5$close_all()
  })
}

ArchR_to_h5 <- function(matrix, h5, matrix.type){
  h5attr(h5, 'matrix_type') <- matrix.type
  # --- get peak names
  peak_names <- stringr::str_c(as.character(seqnames(rowRanges(matrix))),
                               as.character(ranges(rowRanges(matrix))),
                               sep = "_") %>%
    stringr::str_replace_all("-", "_")
  # --- save the matrix
  mat <- assay(matrix)
  stopifnot("Class of assay is not dgCMatrix" = 'dgCMatrix' %in% class(mat))

  h5mat <- h5$create_group("RawCounts")
  h5mat[['values']] <- slot(object = mat, name = 'x')
  h5mat[['indices']] <- slot(object = mat, name = 'i')
  h5mat[['indptr']] <- slot(object = mat, name = 'p')
  h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
  if (!is.null(slot(object = mat, name = 'Dimnames')[[1]])) {
    h5mat[['var_names']] <- slot(object = mat, name = 'Dimnames')[[1]]
  } else {
    h5mat[['var_names']] <- peak_names
  }
  h5mat[['obs_names']] <- slot(object = mat, name = 'Dimnames')[[2]]
  h5attr(h5mat, 'datatype') <- 'SparseMatrix'

  # --- save the cell annotations
  df_to_h5(df = colData(matrix), h5 = h5, gr_name = 'obs')
  # --- save peak information
  if(!is.null(rowData(matrix)@rownames)) {
    df_to_h5(df = rowData(matrix), h5 = h5, gr_name = 'var')
  }
}

#' Data frame to h5
#'
#' Data frame is converted to and saved within the h5 file
#' https://rdrr.io/github/JiekaiLab/RIOH5/src/R/baseIO.R
#' @param df Data frame of cell annotation or gene/peak annotation
#' @param h5 The h5 file
#' @param gr_name The group name represents the property of the data frame.
#'
df_to_h5 <- function(df, h5, gr_name=NULL){
  h5df <- h5$create_group(gr_name)
  h5df[['index']] = rownames(df)
  if(ncol(df)>0){
    h5df[['colnames']] = colnames(df)
  }
  # factor to levels,character to levels,logical to levels
  for(k in names(df)){
    if(is.factor(df[[k]])){
      h5df[[k]]<- as.integer(df[[k]]) - 1L # for 0 begin
      h5df[[paste0(k,'_levels')]]<- levels(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'category'
    }
    if(is.character(df[[k]])){
      str_to_lvl <- factor(df[[k]])
      h5df[[k]]<- as.integer(str_to_lvl) - 1L
      h5df[[paste0(k,'_levels')]]<- levels(str_to_lvl)
      h5attr(h5df[[k]], 'origin_dtype') = 'string'
    }
    if(is.logical(df[[k]])){
      h5df[[k]] <- as.integer(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'bool'
    }
    if(any(is.numeric(df[[k]]),is.integer(df[[k]]))){
      h5df[[k]] <- df[[k]]
      h5attr(h5df[[k]], 'origin_dtype') = 'number'
    }
  }
}

if (exists("snakemake")) {
  cat("Parsing snakemake object...\n")
  threads <- 24
  project_path <- snakemake@input$project_path
  output_file <- snakemake@output$output_file
} else {
  threads <- 6
  project_path <- "out/archr/Atlas/ATAC/proj_peaks"
  output_file <- "out/archr/Atlas/ATAC/peak_mat.h5"
}

cat("threads: ", threads, "\n")
cat("project_path: ", project_path, "\n")
cat("output_file:", output_file, "\n")

print("Set reference genome and load ArchR project")
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

# Extract peak matrix from ArchR project
peak_matrix <- getMatrixFromProject(ArchRProj = proj,
                                    useMatrix = "PeakMatrix",
                                    binarize = FALSE)

# ArchR matrix to h5 file
ArchR_write_h5(matrix = peak_matrix,
               file = output_file,
               matrix.type = "PeakMatrix")
