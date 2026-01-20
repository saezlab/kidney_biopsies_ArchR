library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

cat("Parsing snakemake object...")
threads <- 24
project_path <- snakemake@input$project_path
annotations_path <- snakemake@input$annotations_path
flag_file <- snakemake@output$flag_file

cat("project_path: ", project_path, "\n")
cat("annotations_path: ", annotations_path, "\n")

# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
cat("Loading ArchR project...", "\n")
proj <- loadArchRProject(path = project_path)

# Retrieve the annotations and add them to the ArchR project
cat("Retrieving annotations...", "\n")
cell_annotations <- read.csv(annotations_path, row.names = 1)

cat("Updating annotations...", "\n")
head(cell_annotations)

for (c in colnames(cell_annotations)){
    proj@cellColData[[c]] <- cell_annotations[proj$cellNames, c]
}

# Save ArchR project
print("Saving ArchR project...")
saveArchRProject(ArchRProj =  proj,
                outputDirectory = project_path,
                overwrite = TRUE)

print("Creating flag file")
writeLines("Complete", flag_file)
