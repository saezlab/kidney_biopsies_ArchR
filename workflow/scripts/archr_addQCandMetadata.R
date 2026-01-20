library(SeuratDisk)
library(parallel)
library(ArchR)
library(magrittr)

print("Parsing snakemake object...")
threads <- 24
project_path <- snakemake@input$project_path
metadata <- snakemake@input$metadata
flag_file <- snakemake@output$flag_file

# Set reference genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = threads)
set.seed(123)

proj <- loadArchRProject(path = project_path)

print("Adding disease information...")
patient_metadata <- read.csv(metadata, row.names = 1)
diseases_list <- patient_metadata %>%
  dplyr::transmute(sample_ID = sample_ID,
                   disease = ifelse(is.na(disease), "unknown", disease)) %>%
  tibble::column_to_rownames("sample_ID")
cell_disease_labels <- diseases_list[proj$Sample, ]
proj$diseaseNames <- cell_disease_labels

print("Add information from qualitiative assessment of CellRanger QC plots...")
print("Fragment Size Distribution")
fsd_ideal <- c("K107", "K108", "K113", "K34", "K37", "K69", "K92")
fsd_medium <- c("K11", "K121", "K30", "K35", "K36", "K50", "K96", "K97")
fsd_compromised <- stringr::str_c("K", c(10, 109, 110:112, 114:120,
                                         122:123, 15, 29, 31:33, 45,
                                         51, 53, 9, 93:95, 98:99))
fsd <- data.frame(fsd = c(rep("ideal", length(fsd_ideal)),
                          rep("medium", length(fsd_medium)),
                          rep("compromised", length(fsd_compromised))),
                  sample = c(fsd_ideal, fsd_medium, fsd_compromised))
fsd_dict <- fsd %>% tibble::column_to_rownames("sample")
fsd_qc <- fsd_dict[proj$Sample, ]
proj$FragmentSizeDistribution <- fsd_qc

print("TSS by Unique Fragments")
tss_ideal <- stringr::str_c("K", c(108, 113, 30, 32:37, 69, 92))
tss_medium <- stringr::str_c("K", c(107, 109, 110, 114, 29, 31, 50, 53, 94))
tss_compromised <- stringr::str_c("K", c(10, 11, 111:112, 115:123,
                                         15, 45, 51, 9, 93, 95:99))
tss <- data.frame(tss = c(rep("ideal", length(tss_ideal)),
                          rep("medium", length(tss_medium)),
                          rep("compromised", length(tss_compromised))),
                  sample = c(tss_ideal, tss_medium, tss_compromised))

tss_dict <- tss %>% tibble::column_to_rownames("sample")
tss_qc <- tss_dict[proj$Sample, ]
proj$TSSbyUniqueFragments <- tss_qc

print("Adding information from qualitiative assessment of ArchR QC plots...")
print("Doublet projection R^2")
doublet_ideal <- stringr::str_c("K", c(113:115, 29, 35, 37, 94:95, 98))
doublet_compromised <- stringr::str_c("K", c(10, 107:109, 11, 110:112, 
                                          116:123, 15, 30:34, 36, 45, 
                                          50:51, 53, 69, 9, 92:93, 96:97, 99))
doublet <- data.frame(doublet = c(rep("ideal", length(doublet_ideal)),
                          rep("compromised", length(doublet_compromised))),
                  sample = c(doublet_ideal, doublet_compromised))

doublet_dict <- doublet %>% tibble::column_to_rownames("sample")
doublet_qc <- doublet_dict[proj$Sample, ]
proj$DoubletProjection <- doublet_qc

print("Median Fragments")
median_fragments_ideal <- stringr::str_c("K", c(10, 107, 11, 113, 116:119,
                                               122:123, 15, 37, 45, 9))
median_fragments_compromised <- stringr::str_c("K", c(108:109, 110:112,
                                                     114:115, 120:121, 29:36,
                                                     50, 51, 53, 69, 92:99))
median_fragments <- data.frame(
        median_fragments = c(rep("ideal",
                             length(median_fragments_ideal)),
                             rep("compromised",
                             length(median_fragments_compromised))),
        sample = c(median_fragments_ideal, median_fragments_compromised))

median_fragments_dict <- median_fragments %>%
                                tibble::column_to_rownames("sample")
mf_qc <- median_fragments_dict[proj$Sample, ]
proj$median_fragmentsQC <- mf_qc

# Save ArchR project
print("Saving ArchR project...")
saveArchRProject(ArchRProj =  proj,
                outputDirectory = project_path,
                overwrite = TRUE)

print("Creating flag file")
writeLines("Complete", flag_file)
