library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

print("Parsing snakemake object...")
threads <- snakemake@params$threads
project_path <- snakemake@input$project_path
seurat_obj <- snakemake@input$seurat_obj
flag_file <- snakemake@output$flag_file
label_transfer_csv <- snakemake@output$label_transfer_csv

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

proj <- loadArchRProject(path = project_path)
seRNA <- LoadH5Seurat(seurat_obj)

print("Adding ArchR gene integration matrix...")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "HarmonyLeiden",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "annotation",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

# Write label transfer annotations to csv
print("Writing label transfer annotations to csv...")
pred_annotations <- proj$predictedGroup_Un
write.csv(pred_annotations, label_transfer_csv, row.names = FALSE)

# Save ArchR project
print("Saving ArchR project...")
saveArchRProject(ArchRProj = proj,
                 outputDirectory = project_path,
                 overwrite = TRUE)

print("Creating flag file")
writeLines("Complete", flag_file)