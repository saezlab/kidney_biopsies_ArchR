library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)
library(BSgenome) # Necessary for footprinting
library(BSgenome.Hsapiens.UCSC.hg38)

cat("Parsing snakemake object...", "\n")
threads <- 24
project_path <- snakemake@input$project_path
pseudobulk_rds <- snakemake@output$pseudobulk_rds
motif_footprints_rds <- snakemake@output$motif_footprints_rds
motif_footprints_plots_rds <- snakemake@output$motif_footprints_plots_rds

cat("project_path: ", project_path, "\n")

# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

# Check group coverages
proj@projectMetadata$GroupCoverages$cluster_label

# Obtain motif positions (uses peakAnnotation slot by default)
motifPositions <- getPositions(proj)
motifPositions

# TFs of interest
motifs <- c("HNF4G", "HNF4A", "RXRG", "RXRB", "TEAD4", "TEAD1", "TEAD3", "TFCP2L1", "RELA", "SCX", "MEF2C")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
# We may need to remove some motifs which are found by grep but not 
# what we actually want, as in the ArchR documentation example
# markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "cluster_label"
)

plots <- plotFootprints(seFoot, plot = FALSE)

saveRDS(seFoot, file = motif_footprints_rds)
saveRDS(plots, file = motif_footprints_plots_rds)