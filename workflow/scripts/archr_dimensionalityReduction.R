library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

print("Parsing snakemake object...")
threads <- snakemake@params$threads
project_path <- snakemake@input$project_path
flag_file <- snakemake@output$flag_file
coordinates_csv <- snakemake@output$coordinates_csv
clusters_csv <- snakemake@output$clusters_csv

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

proj <- loadArchRProject(path = project_path)

# Dimensionality reduction
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSILeiden",
    iterations = 3,
    clusterParams = list( #See Seurat::FindClusters
        resolution = 0.4,
        sampleCells = 10000,
        n.start = 10,
        algorithm = 4 # Leiden algorithm on Pau's advice
    ),
    varFeatures = 200000,
    dimsToUse = 1:30
)

# Harmony batch effect correction
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSILeiden",
  name = "HarmonyLeiden",
  groupBy = "Sample"
)

# Clustering
proj <- addClusters(
  input = proj,
  reducedDims = "HarmonyLeiden",
  method = "Seurat",
  name = "HarmonyLeidenClusters",
  resolution = 0.2
)

# UMAP
proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = "HarmonyLeiden",
    name = "HarmonyLeidenUMAP",
    nNeighbors = 30,
    minDist = 0.1,
    metric = "cosine"
)

# Write embedding coordinates to csv
getCoordinates <- function(embedding) {
  embedding$df %>% setNames(c("X", "Y")) %>% tibble::rownames_to_column("cell")
}

coordinates <- lapply(proj@embeddings, getCoordinates) %>%
  setNames(names(proj@embeddings)) %>%
  dplyr::bind_rows(.id = "embedding")

write.csv(coordinates, coordinates_csv)

# Write clusters to csv
cellColData <- proj@cellColData
cluster_list <- cellColData[stringr::str_detect(names(cellColData),
                                                "Clusters")]
clusters <- as.data.frame(cluster_list)

write.csv(clusters, clusters_csv)

# UMAP plots
p1 <- plotEmbedding(ArchRProj = proj,
                    colorBy = "cellColData",
                    name = "Sample",
                    embedding = "HarmonyLeidenUMAP")

p2 <- plotEmbedding(ArchRProj = proj,
                    colorBy = "cellColData",
                    name = "diseaseNames",
                    embedding = "HarmonyLeidenUMAP")

p3 <- plotEmbedding(ArchRProj = proj,
                    colorBy = "cellColData",
                    name = "HarmonyLeidenClusters",
                    embedding = "HarmonyLeidenUMAP")

plotPDF(p1, p2, p3,
        name = "Plot-UMAP-HarmonyLeiden.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 8,
        height = 8)

# Save ArchR project
print("Saving ArchR project...")
saveArchRProject(ArchRProj =  proj,
                outputDirectory = project_path,
                overwrite = TRUE)

print("Creating flag file")
writeLines("Complete", flag_file)
