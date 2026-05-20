library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

print("Parsing snakemake object...")
threads <- snakemake@params$threads
project_path <- snakemake@input$project_path
harmony_lsi_csv <- snakemake@output$harmony_lsi_csv

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

print("Loading ArchR project...")
proj <- loadArchRProject(path = project_path)

print("Retrieving Harmony LSI...")
# Retrieve dimensionality reduction
harmony_lsi <- getReducedDims(proj,
 reducedDims = "HarmonyLeiden",
 returnMatrix = TRUE)

print(head(harmony_lsi))
write.csv(harmony_lsi, harmony_lsi_csv)
