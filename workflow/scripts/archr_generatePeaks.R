library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)
library(BSgenome.Hsapiens.UCSC.hg38)

cat("Parsing snakemake object...")
threads <- 24
if (exists("snakemake")){
    project_path <- snakemake@input$project_path
    name_annotations <- snakemake@params$name_annotations
    name_clusters <- snakemake@params$name_clusters
    output_dir <- snakemake@output$output_dir
    pathToMacs2 <- findMacs2()
} else {
    if (Sys.info()["sysname"] == "Darwin"){
        project_path <- "out/archr/Atlas/ATAC/proj_2"
        name_annotations <- "level1"
        name_clusters <- "cluster"
        output_dir <- "out/archr/Atlas/ATAC/proj_peaks"
        pathToMacs2 <- "/Users/charlotteboys/SOFTWARE/mambaforge/bin/macs2"
    }
}

cat("project_path: ", project_path)
cat("name_annotations: ", name_annotations)
cat("name_clusters: ", name_clusters)
cat("output_dir: ", output_dir)
cat("path to MACS2: ", pathToMacs2)
# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

cat("Removing cells marked for removal...", "\n")
proj <- proj[!proj@cellColData[[name_annotations]] %in% c("remove"), ]

proj <- addGroupCoverages(proj,
                          groupBy = name_clusters) #, threads = 1)

proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = name_clusters,
    pathToMacs2 = pathToMacs2)

proj <- addPeakMatrix(proj)

saveArchRProject(ArchRProj = proj,
                 outputDirectory = output_dir,
                 )