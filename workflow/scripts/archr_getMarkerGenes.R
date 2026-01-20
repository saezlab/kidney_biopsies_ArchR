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

cat("project_path: ", project_path, "\n")
cat("group_by: ", group_by, "\n")

# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

# Retrieve marker genes
markers <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# Save marker genes
marker_list <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")
saveRDS(marker_list, file = snakemake@output$marker_rds)