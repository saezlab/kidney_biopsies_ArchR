library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

print("Parsing snakemake object...")
threads <- 24
files <- snakemake@input$files
output_dir <- snakemake@output$output_dir
arrow_files_dir <- snakemake@params$arrow_files_dir
qc_dir <- snakemake@output$qc_dir

cat("files: ", files, "\n")
cat("output_dir: ", output_dir, "\n")
cat("arrow_files_dir: ", arrow_files_dir, "\n")
cat("qc_dir: ", qc_dir, "\n")

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

print("Samples:")
samplesWithATAC <- files %>%
    setNames(stringr::str_extract(string = files,
                                  pattern = "K[:digit:]+"))

print(samplesWithATAC)

print("Desired Arrow file destinations:")
outputNames <- stringr::str_c(
    arrow_files_dir,
    names(samplesWithATAC))

print(outputNames)

# Attempt to solve error
# <simpleError in H5Fcreate(file): HDF5.
#             File accessibility. Unable to open file.>
# when creating arrow files on the cluster, by using dev version of ArchR and
# adding ArchR locking
# See https://github.com/GreenleafLab/ArchR/issues/248
addArchRLocking(locking = TRUE)

ArrowFiles <- createArrowFiles(
  inputFiles = samplesWithATAC,
  sampleNames = names(samplesWithATAC),
  outputNames = stringr::str_c(
    arrow_files_dir,
    names(samplesWithATAC)),
  QCDir = qc_dir,
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  subThreading = FALSE,
  threads = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = output_dir,
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# Save ArchR project
saveArchRProject(ArchRProj =  proj, outputDirectory = output_dir)
