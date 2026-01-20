library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)
library(BSgenome) # Necessary for footprinting
library(BSgenome.Hsapiens.UCSC.hg38)

cat("Parsing snakemake object...", "\n")
threads <- 24
project_path <- snakemake@input$project_path
group_by <- snakemake@params$group_by
flag_file <- snakemake@output$flag_file

cat("project_path: ", project_path, "\n")
cat("group_by: ", group_by, "\n")

# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

# Add group coverages
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = group_by)

saveArchRProject(
    ArchRProj = proj,
    outputDirectory = project_path,
    overwrite = TRUE
    )

print("Creating flag file")
writeLines("Complete", flag_file)