options(timeout=180) #prevent errors due to internet connection

remotes::install_github("GreenleafLab/ArchR", ref="dev", repos = BiocManager::repositories(), upgrade = "never")

library(ArchR)
ArchR::installExtraPackages()

# Requires hdf5r, which requires hdf5  v1.8.13
remotes::install_github("mojaveazure/seurat-disk")