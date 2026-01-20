library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

print("Parsing snakemake object...")
threads <- 8 # To prevent memory problems
project_path <- snakemake@input$project_path
qc_dir <- snakemake@output$qc_dir
output_dir <- snakemake@output$output_dir
flag_file <- snakemake@output$flag_file
lsi_method <- as.integer(snakemake@params$lsi_method)
proj_qc_dir <- snakemake@params$proj_qc_dir

cat("project_path: ", project_path, "\n")
cat("qc_dir: ", qc_dir, "\n")
cat("output_dir: ", output_dir, "\n")
cat("flag_file: ", flag_file, "\n")
cat("lsi_method: ", lsi_method, "\n")
cat("proj_qc_dir: ", proj_qc_dir, "\n")

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

proj <- loadArchRProject(path = project_path)

# Create QC directory for each sample
samples <- list.dirs(proj_qc_dir,
                      full.names = FALSE,
                      recursive = FALSE)

sample_qc_dirs <- stringr::str_c(qc_dir, samples, sep = "/")
lapply(sample_qc_dirs, dir.create, recursive = TRUE)

proj <- addDoubletScores(
  input = proj,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = lsi_method,
  UMAPParams = list(n_neighbors = 10,
                    min_dist = 0.01,
                    metric = "euclidean",
                    verbose = FALSE),
  outDir = qc_dir
)

proj <- filterDoublets(ArchRProj = proj)

# Save ArchR project
saveArchRProject(ArchRProj =  proj, outputDirectory = output_dir)

print("Creating flag file")
writeLines("Complete", flag_file)
