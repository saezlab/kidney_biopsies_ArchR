library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

cat("Parsing snakemake object...", "\n")
threads <- 24
project_path <- snakemake@input$project_path
group_by <- snakemake@params$group_by
enrich_motifs_rds <- snakemake@output$enrich_motifs_rds
flag_file <- snakemake@output$flag_file

cat("project_path: ", project_path, "\n")
cat("group_by: ", group_by, "\n")
cat("enrich_motifs_rds: ", enrich_motifs_rds, "\n")

# Set reference genome and load ArchR project
addArchRGenome("hg38")
addArchRThreads(threads = threads)
proj <- loadArchRProject(path = project_path)

print("Cell counts:")
table(proj$cluster_label)
print("NA counts:")
sum(is.na(proj$cluster_label))

# 06.25 No longer necessary since label transfer v3:
# Remove muscle_mac, PL_B (too few cells, causes "Cloud has no points" error in getMarkerFeatures
# idxPass <- which(!proj$clean_label %in% c(NA, "muscle_mac", "PL_B"))
# cellsPass <- proj$cellNames[idxPass]
# proj[cellsPass, ]
# unique(proj[cellsPass,]$clean_label)

# Marker peaks
marker_peaks <- getMarkerFeatures(
    ArchRProj = proj,#proj[cellsPass,],
    useMatrix = "PeakMatrix",
    groupBy = group_by,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

proj <- addMotifAnnotations(
    ArchRProj = proj,
    motifSet = "cisbp",
    name = "Motif",
    force = TRUE
)

enrich_motifs <- peakAnnoEnrichment(
    seMarker = marker_peaks,
    ArchRProj = proj, #proj[cellsPass,],
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

saveArchRProject(
    ArchRProj = proj,
    outputDirectory = project_path,
    overwrite = TRUE
    )

saveRDS(enrich_motifs, file = enrich_motifs_rds)

print("Creating flag file")
writeLines("Complete", flag_file)
