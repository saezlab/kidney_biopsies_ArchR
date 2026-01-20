library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

cat("Parsing snakemake object...", "\n")
threads <- 24
project_path <- snakemake@input$project_path
group_by <- snakemake@params$group_by
pseudobulk_rds <- snakemake@output$pseudobulk_rds
flag_file <- snakemake@output$flag_file
# celltype_motif_enrichment_csv <- snakemake@output$celltype_motif_enrichment_csv
# motif_top_50_csv <- snakemake@output$motif_top_50_csv

cat("project_path: ", project_path, "\n")
cat("group_by: ", group_by, "\n")
cat("pseudobulk_rds: ",pseudobulk_rds, "\n")

# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

# Create sample x cell type column
proj <- addCellColData(ArchRProj = proj,
                       data = paste0(proj@cellColData$Sample,"_",proj@cellColData[[group_by]]),
                       name = "celltypexsample",
                       cells = getCellNames(proj))

# Save pseudobulks as a summarised experiment object
groupSE <- getGroupSE(ArchRProj = proj,
                      useMatrix = "PeakMatrix",
                      groupBy = "celltypexsample")

saveRDS(groupSE, file = pseudobulk_rds)

print("Creating flag file")
writeLines("Complete", flag_file)
