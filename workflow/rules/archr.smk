
import glob
import re
# Find sample IDS for which we have ATAC data
SAMPLES = glob.glob("data/Samples/*/ATAC/fragments.tsv.gz")
IDS = [re.search(r"(?<=Samples/).*?(?=/ATAC)", k).group(0) for k in SAMPLES]

# This rule creates arrow files, the ArchR project,
# performs doublet filtering
rule ATAC_createArchRproject:
    input:
        files = expand("data/Samples/{id}/ATAC/fragments.tsv.gz",
            id = IDS)
    output:
        output_dir = directory('out/archr/Atlas/ATAC/proj_unfiltered'),
        qc_dir = directory('out/archr/Atlas/ATAC/QualityControl')
    params:
        arrow_files_dir = directory('out/archr/Atlas/ATAC/ArrowFiles')
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '5:00:00'
    script:
        '../scripts/archr_createArchRproject.R'
        
rule ATAC_doubletFiltering:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_unfiltered'
    output:
        output_dir = directory('out/archr/Atlas/ATAC/proj_{lsi_method}'),
        qc_dir = directory('out/archr/Atlas/ATAC/QualityControl_{lsi_method}'),
        flag_file = 'out/archr/Atlas/ATAC/flags/doubletFiltering_{lsi_method}.txt'
    params:
        lsi_method = '{lsi_method}',
        proj_qc_dir = directory('out/archr/Atlas/ATAC/QualityControl')
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '5:00:00'
    script:
        '../scripts/archr_doubletFiltering.R'

# Adds patient and QC metadata to ArchR project
rule ATAC_addQCandMetadata:
    input:
        project_path = ancient('out/archr/Atlas/ATAC/proj_2'),
        metadata = "data/clinical_data/patient_metadata.csv",
        flag_file = 'out/archr/Atlas/ATAC/flags/doubletFiltering_2.txt'
    output:
        # Would like to overwrite the input ArchR project, but 
        # can't have the same input and output.
        # So we introduce a 'flag file'
        # Note that this rule also overwrites the input project
        flag_file = 'out/archr/Atlas/ATAC/flags/addQCandMetadata.txt'
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '1:00:00'
    script:
        '../scripts/archr_addQCandMetadata.R'

rule ATAC_dimensionalityReduction:
    input:
        project_path = ancient('out/archr/Atlas/ATAC/proj_2'),
        flag_file = 'out/archr/Atlas/ATAC/flags/addQCandMetadata.txt'
    output:
        # Note that this rule also overwrites the input project
        coordinates_csv = 'out/archr/Atlas/ATAC/coordinates.csv',
        clusters_csv = 'out/archr/Atlas/ATAC/clusters.csv',
        flag_file = 'out/archr/Atlas/ATAC/flags/dimensionalityReduction.txt'
    params:
        threads = 24
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '2:00:00'
    script:
        '../scripts/archr_dimensionalityReduction.R'

rule ATAC_labelTransfer:
    input:
        project_path = ancient('out/archr/Atlas/ATAC/proj_2'),
        flag_file = 'out/archr/Atlas/ATAC/flags/dimensionalityReduction.txt',
        seurat_obj = 'out/exp/Atlas/RNA/seuratobject_minimal.h5Seurat'
    output:
        # Note that this rule also overwrites the input project
        label_transfer_csv = 'out/archr/Atlas/ATAC/label_transfer_annotations.csv',
        flag_file = 'out/archr/Atlas/ATAC/flags/labelTransfer.txt'
    params:
        threads = 24
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '4:00:00'
    script:
        '../scripts/archr_labelTransfer.R'

rule ATAC_extractGeneScoreMatrix:
    input:
        project_path = ancient('out/archr/Atlas/ATAC/proj_2'),
        flag_file = 'out/archr/Atlas/ATAC/flags/labelTransfer.txt',
        support_functions = 'workflow/src/archr_support_functions.R'
    output:
        output_file = 'out/archr/Atlas/ATAC/gene_score_matrix.h5'
    params:
        threads = 24
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '2:00:00'
    script:
        '../scripts/archr_extractGeneScoreMatrix.R'

rule ATAC_geneScoreMatrixToH5ad:
    input:
        mat_h5 = 'out/archr/Atlas/ATAC/gene_score_matrix.h5',
        embeddings_path = 'out/archr/Atlas/ATAC/coordinates.csv'
    output:
        mat_h5ad = 'out/archr/Atlas/ATAC/gene_score_matrix.h5ad'
    params:
        group = 'GeneScores'
    conda:
        '../envs/scanpy.yaml'
    threads: 24
    resources:
        mem_mb = 50000,
        time = '2:00:00'
    script:
        '../scripts/archr_h5_to_h5ad.py'

rule ATAC_addAnnotations:
    input:
        project_path = ancient('out/archr/Atlas/ATAC/proj_2'),
        annotations_path = 'out/archr/Atlas/ATAC/annotations.csv'
    output:
        flag_file = 'out/archr/Atlas/ATAC/flags/addAnnotations.txt'
    conda:
        '../envs/ArchR.yaml'
    threads: 24
    resources:
        mem_mb = 30000,
        time = '2:00:00'
    script:
        '../scripts/archr_addAnnotations.R'

rule ATAC_generatePeaks:
    input:
        project_path = ancient('out/archr/Atlas/ATAC/proj_2'),
        flag_file = 'out/archr/Atlas/ATAC/flags/addAnnotations.txt'
    output:
        output_dir = directory('out/archr/Atlas/ATAC/proj_peaks2')
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    params:
       name_annotations = 'level1',
       name_clusters = 'cluster'
    threads: 24
    resources:
        mem_mb = 90000,
        time = '3:00:00'
    script:
        '../scripts/archr_generatePeaks.R'

rule ATAC_addLabelTransfer:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2',
        annotations_path = 'data/annotations/label_transfer_RNA_ATAC_V3_clean.csv'
    output:
        flag_file = 'out/archr/Atlas/ATAC/flags/addLabelTransfer.txt'
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    threads: 24
    resources:
        mem_mb = 30000,
        time = '2:00:00'
    script:
        '../scripts/archr_addAnnotations.R'

rule ATAC_motifEnrichment:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2',
        flag_file = 'out/archr/Atlas/ATAC/flags/addLabelTransfer.txt'
    output:
        # Note that this rule also overwrites the input project
        flag_file = 'out/archr/Atlas/ATAC/flags/motifEnrichment.txt',
        enrich_motifs_rds = 'out/archr/Atlas/ATAC/enrich_motifs.rds'
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    params:
        group_by = 'cluster_label'
    threads: 24
    resources:
        mem_mb = 30000,
        time = '4:00:00'
    script:
        '../scripts/archr_motifEnrichment.R'

rule ATAC_getMarkerGenes:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2',
        flag_file = 'out/archr/Atlas/ATAC/flags/addLabelTransfer.txt'
    output:
        marker_rds = 'out/archr/Atlas/ATAC/marker_genes.rds'
    params:
        group_by = 'cluster_label'
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    threads: 24
    resources:
        mem_mb = 30000,
        time = '2:00:00'
    script:
        '../scripts/archr_getMarkerGenes.R'

rule ATAC_groupCoverages:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2',
        flag_file = 'out/archr/Atlas/ATAC/flags/addLabelTransfer.txt'
    output:
        # Note that this rule also overwrites the input project
        flag_file = 'out/archr/Atlas/ATAC/flags/groupCoverages.txt',
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    params:
        group_by = 'cluster_label'
    threads: 24
    resources:
        mem_mb = 100000,
        runtime=480
    script:
        '../scripts/archr_groupCoverages.R'

rule ATAC_motifFootprinting:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2',
        flag_file = 'out/archr/Atlas/ATAC/flags/groupCoverages.txt'
    output:
        motif_footprints_rds = 'out/archr/Atlas/ATAC/motiffootprint_cluster.rds',
        motif_footprints_plots_rds = 'out/archr/Atlas/ATAC/motiffootprint_cluster_plots.rds'
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    threads: 24
    resources:
        mem_mb = 40000,
        runtime=30
    script:
        '../scripts/archr_motifFootprinting_celltype.R'

rule ATAC_pseudobulk:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2',
        flag_file = 'out/archr/Atlas/ATAC/flags/addLabelTransfer.txt'
    output:
        pseudobulk_rds = 'out/archr/Atlas/ATAC/pseudobulk.rds',
        flag_file = 'out/archr/Atlas/ATAC/flags/pseudobulk.txt'
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    params:
        group_by = 'clean_label'
    threads: 24
    resources:
        mem_mb = 90000,
        time = '4:00:00'
    script:
        '../scripts/archr_pseudobulk.R'

rule ATAC_retrievePeakMatrix:
    input:
        project_path = 'out/archr/Atlas/ATAC/proj_peaks2'
    output:
        output_file = 'out/archr/Atlas/ATAC/peak_mat2.rds' # .h5'
    singularity:
        'workflow/envs/archr_1.0.3-base-r4.1_0.0.1.sif'
    threads: 24
    resources:
        mem_mb = 90000,
        time = '2:00:00'
    script:
        '../scripts/archr_retrievePeakMatrix.R'
        
rule ATAC_peak_to_h5ad:
    input:
        peak_mat_h5 = 'out/archr/Atlas/ATAC/peak_mat2.h5'
    output:
        peak_mat_h5ad = 'out/archr/Atlas/ATAC/peak_mat2.h5ad'
    singularity:
        '../singularities/Robin_sings/scanpy.0.0.9.sif'
    threads: 8
    resources:
        mem_mb = 10000,
        time = '1:00:00'
    script:
        '../scripts/archr_peak_to_h5ad.py'