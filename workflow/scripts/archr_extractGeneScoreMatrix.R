library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)
library(hdf5r)

print("Parsing snakemake object...")
threads <- snakemake@params$threads
project_path <- snakemake@input$project_path
support_functions <- snakemake@input$support_functions
flag_file <- snakemake@output$flag_file
output_file <- snakemake@output$output_file

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

# Load ArchR project
proj <- loadArchRProject(path = project_path)

# Load support functions
source(support_functions)

# Extract gene score matrix
gene_score_matrix <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads()
)

# Write to h5
ArchR_write_h5(matrix = gene_score_matrix,
               file = output_file,
               matrix.type = "GeneScoreMatrix")