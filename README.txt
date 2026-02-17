================================================================================
EPILEPSY RNA-SEQ ANALYSIS PIPELINE - SCRIPT DOCUMENTATION
================================================================================

OVERVIEW
--------
This repository contains analysis scripts for bulk RNA-seq and single-cell RNA-seq 
data from epilepsy research, including human and animal (mouse/rat) datasets. The 
pipeline covers the complete workflow from raw data preprocessing to advanced 
downstream analyses including differential expression, enrichment, deconvolution, 
and cross-species comparisons.

PROJECT STRUCTURE
-----------------
The scripts are organized into four main directories:
1. RNA_Human/     - Bulk RNA-seq analysis for human samples
2. RNA_Animal/   - Bulk RNA-seq analysis for mouse and rat samples
3. scRNA_Human/  - Single-cell RNA-seq analysis for human samples
4. scRNA_Animal/ - Single-cell RNA-seq analysis for mouse sample

SETUP
-----
Before running any scripts, source the configuration file:
  source Scratch_settings.sh

This file defines:
- Project directories and paths
- Reference genome paths (GRCh38, GRCm39, mRatBN7.2)
- Sample lists and metadata paths
- Output directory structures

GENE SET PREPARATION
--------------------
00.Geneset.R
  Purpose: Prepare gene sets for enrichment analysis and cross-species comparison
  Functions:
    - Ortholog mapping (Human-Mouse, Human-Rat, Human-Mouse-Rat)
    - BioMart GO term extraction and organization
    - Gene set creation for 9 categories: IC, NT, ND, NI, TNF, MF, NG, GG, SP
    - Upper category grouping and union gene sets (W set)
  Usage: Source this script before running enrichment analyses


================================================================================
RNA_HUMAN ANALYSIS PIPELINE
================================================================================

00. DATA PREPARATION
--------------------
00.Matrix_H.R
  Purpose: Import and prepare count/TPM matrices and metadata

00.Subgrouping_H.R
  Purpose: Define sample subgroups and filter datasets
  Functions:
    - Filter by brain location (HIPPO)
    - Create diagnosis groups (MTLEALL, MTLEHS, MTLE, NL)
    - Remove outlier samples
    - Generate filtered datasets for specific comparisons

01. PRE-PROCESSING
------------------
01.Pre-processing/03.fastp_paired.sh
  Purpose: Quality control and adapter trimming for paired-end reads

01.Pre-processing/03.fastp_single.sh
  Purpose: Quality control and adapter trimming for single-end reads

01.Pre-processing/04.Indexing.sh
  Purpose: Build STAR index for alignment

01.Pre-processing/04.aligned_paired_*.sh
  Purpose: Align paired-end reads to reference genome
  Scripts for different read lengths: 100bp, 125bp, 150bp

01.Pre-processing/04.aligned_single.sh
  Purpose: Align single-end reads to reference genome

02. QUANTIFICATION
------------------
02.Quantification/01.HTseq.sh
  Purpose: Count reads per gene using HTseq-count

02.Quantification/02.gtf_table.R
  Purpose: Extract gene information from GTF file

02.Quantification/03.Matrix.R
  Purpose: Merge individual count files into count matrix
  Functions:
    - Read all HTseq count files
    - Merge by gene name
    - Filter out non-gene features
    - Align with sample metadata

02.Quantification/04.combat.R
  Purpose: Batch effect correction using ComBat

03. NORMALIZATION
-----------------
03.Normalization/01.TPM_normalization.R
  Purpose: Calculate TPM (Transcripts Per Million) normalization

04. CLUSTERING
--------------
04.Clustering/01.Count_PCA.R
  Purpose: Principal Component Analysis on count data

04.Clustering/02.TPM_PCA.R
  Purpose: Principal Component Analysis on TPM data

05. COMPARISON
--------------
05.Comparison/01.TPM_Comparison.R
  Purpose: Compare expression levels across sample groups

06. SELECTION
------------
06.Selection/01.TPM_Outlier_Selection.R
  Purpose: Identify outlier samples based on TPM expression

06.Selection/02.TPM_Outlier_plot.R
  Purpose: Visualize TPM-based outlier samples

06.Selection/03.PCA_Outlier_Selection.R
  Purpose: Identify outlier samples based on PCA

06.Selection/04.PCA_Outlier_plot.R
  Purpose: Visualize PCA-based outliers

06.Selection/05.TPM_PCA_Outlier_Selection.R
  Purpose: Combined TPM and PCA outlier detection

06.Selection/06.TPM_PCA_Outlier_plot.R
  Purpose: Visualize combined outliers

07. DIFFERENTIAL EXPRESSION
----------------------------
07.DESeq/01.DESeq.R
  Purpose: Differential expression analysis using DESeq2
  Functions:
    - Create DESeqDataSet
    - Pre-filtering (rowSums >= 200)
    - DESeq analysis

07.DESeq/02.DESeq_ALL.R
  Purpose: Run DESeq2 for all comparison groups

08. ENRICHMENT ANALYSIS
-----------------------
08.Enrichment/ORA_01.Matrix_DESeq.R
  Purpose: Prepare DESeq results for Over-Representation Analysis (ORA)

08.Enrichment/ORA_02.ToppFun.R
  Purpose: Perform ORA using ToppFun

08.Enrichment/GSEA_01.Matrix_GTF.R
  Purpose: Prepare GTF-based expression matrix for GSEA
  Functions:
    - Extract gene ID and gene name from GTF
    - Convert gene name to gene ID in expression matrix
    - Aggregate by gene ID

08.Enrichment/GSEA_02.Matrix_GMT.R
  Purpose: Prepare GMT format gene sets for GSEA

08.Enrichment/GSEA_03.Matrix_TPM_hippo.R
  Purpose: Prepare TPM matrix for GSEA (hippocampus-specific)

08.Enrichment/GSEA_04.Output_GOBP.R
  Purpose: Extract and format GSEA GO Biological Process results

08.Enrichment/GSEA_05.Output_GOBP_combined.R
  Purpose: Combine multiple GSEA results

09. DECONVOLUTION
-----------------
09.Deconvolution/01.IMPORT_MIXTURE.R
  Purpose: Prepare TPM matrix for CIBERSORTx deconvolution
  Functions:
    - Add GeneSymbol header if missing
    - Format matrix for CIBERSORTx input

09.Deconvolution/02.IMPORT_GO.R
  Purpose: Import gene sets for deconvolution signature

09.Deconvolution/03.CIBER_01234Docker.sh
  Purpose: Run CIBERSORTx deconvolution analysis
  Functions:
    - High-resolution cell type deconvolution
    - Multiple CIBERSORTx modes

09.Deconvolution/04.CIBER_04HiRes_Merge.R
  Purpose: Merge high-resolution CIBERSORTx results

09.Deconvolution/05.DEG_Matrix_Integ.R
  Purpose: Prepare differential expression matrices for deconvolved cell types

09.Deconvolution/06.DEG_DESeq_Integ.R
  Purpose: DESeq2 analysis for deconvolved cell types

09.Deconvolution/07.DEG_DESeq_ALL_Integ.R
  Purpose: Run DESeq2 for all cell types and comparisons

09.Deconvolution/08.Correlation_Celltype.R
  Purpose: Correlation analysis between cell type proportions

10. CORRELATION ANALYSIS
------------------------
10.Correlation/01.Correlation_WholeGene.R
  Purpose: Correlation analysis using all genes
  Functions:
    - Cross-species correlation (Human vs Mouse/Rat)
    - Ortholog mapping
    - Correlation calculation and visualization

10.Correlation/01.Correlation_Epileptogenesis.R
  Purpose: Correlation analysis for epileptogenesis-related genes

10.Correlation/02.Dotplot_correlation_WholeGene.R
  Purpose: Create dot plots for correlation results (all genes)

10.Correlation/02.Dotplot_correlation_Epileptogenesis.R
  Purpose: Create dot plots for epileptogenesis correlations

10.Correlation/02.Dotplot_correlation_WholeGene_phase.R
  Purpose: Phase-specific correlation dot plots

10.Correlation/02.Dotplot_correlation_Epileptogenesis_phase.R
  Purpose: Phase-specific epileptogenesis correlation plots

10.Correlation/03.Dotplot_shared_genes_WholeGene.R
  Purpose: Visualize shared significant genes across species

10.Correlation/03.Dotplot_shared_genes_Epileptogenesis.R
  Purpose: Shared epileptogenesis genes visualization

10.Correlation/04.Celltype_Corr_Matrix_WholeGene.R
  Purpose: Cell type correlation matrix (all genes)

10.Correlation/04.Celltype_Corr_Matrix_Epileptogenesis.R
  Purpose: Cell type correlation matrix (epileptogenesis)

10.Correlation/04.Dotplot_Celltype_WholeGene.R
  Purpose: Cell type-specific dot plots (all genes)

10.Correlation/04.Dotplot_Celltype_Epileptogenesis.R
  Purpose: Cell type-specific dot plots (epileptogenesis)

10.Correlation/05.GO_Corr_Matrix.R
  Purpose: GO term correlation matrix

10.Correlation/05.GO_Dotplot.R
  Purpose: GO term correlation dot plots

11. SAMPLE VARIABLES
--------------------
11.Sample_variables/00.sankey_plot.R
  Purpose: Create Sankey diagrams for sample metadata

11.Sample_variables/01.matrix.R
  Purpose: Prepare matrices for sample variable analysis

11.Sample_variables/02.deseq.R
  Purpose: DESeq2 analysis for sample variables

11.Sample_variables/02.deseq_ALL.R
  Purpose: DESeq2 for all sample variables

11.Sample_variables/03.correlation.R
  Purpose: Correlation analysis with sample variables

11.Sample_variables/04.correlation_additional.R
  Purpose: Additional correlation analyses


================================================================================
RNA_ANIMAL ANALYSIS PIPELINE
================================================================================

00. DATA PREPARATION
--------------------
00.Matrix_MR.R
  Purpose: Import and prepare matrices for Mouse and Rat
  Functions:
    - Load count/TPM matrices for Mouse (M) and Rat (R)
    - Import DESeq results for both species
    - Load deconvolution results

00.Subgrouping_MR.R
  Purpose: Define subgroups for animal datasets
  Functions:
    - Filter by treatment type (KAI, PILO, TBI, AMG, PPS)
    - Group by time point and brain region
    - Create comparison groups

01. PRE-PROCESSING
------------------
01.Pre-processing/03.fastp_paired.sh
  Purpose: Quality control for paired-end animal reads

01.Pre-processing/03.fastp_single.sh
  Purpose: Quality control for single-end animal reads

01.Pre-processing/04.Indexing_Mouse.sh
  Purpose: Build STAR index for mouse (GRCm39)

01.Pre-processing/04.Indexing_Rat.sh
  Purpose: Build STAR index for rat (mRatBN7.2)

01.Pre-processing/04.aligned_M_Paired_*.sh
  Purpose: Align paired-end mouse reads
  Scripts for different read lengths: 35bp, 50bp, 100bp, 125bp, 150bp

01.Pre-processing/04.aligned_M_Single_*.sh
  Purpose: Align single-end mouse reads
  Scripts for 50bp and 75bp reads

01.Pre-processing/04.aligned_R_Paired_*.sh
  Purpose: Align paired-end rat reads
  Scripts for 75bp, 125bp, 150bp reads

01.Pre-processing/04.aligned_R_Single_35.sh
  Purpose: Align single-end rat reads (35bp)

02. QUANTIFICATION
------------------
02.Quantification/01.HTseq_M.sh
  Purpose: Count reads per gene for mouse samples

02.Quantification/01.HTseq_R.sh
  Purpose: Count reads per gene for rat samples

02.Quantification/02.gtf_table.R
  Purpose: Extract gene information from GTF (Mouse/Rat)

02.Quantification/03.Matrix.R
  Purpose: Merge count files into matrices for Mouse and Rat

02.Quantification/04.combat.R
  Purpose: Batch effect correction for animal data

03. NORMALIZATION
-----------------
03.Normalization/01.TPM_normalization.R
  Purpose: TPM normalization for Mouse and Rat

04. CLUSTERING
--------------
04.Clustering/01.Count_PCA_M.R
  Purpose: PCA on mouse count data

04.Clustering/01.Count_PCA_R.R
  Purpose: PCA on rat count data

04.Clustering/01.TPM_PCA_M.R
  Purpose: PCA on mouse TPM data

04.Clustering/01.TPM_PCA_R.R
  Purpose: PCA on rat TPM data

05. DIFFERENTIAL EXPRESSION
----------------------------
05.DESeq/02.DESeq_M_*.R
  Purpose: DESeq2 analysis for mouse datasets
  Scripts for different comparisons:
    - DESeq_M_ALL.R: All mouse comparisons
    - DESeq_M_PILO.R: Pilocarpine model
    - DESeq_M_KAI_IH.R: Kainic acid intrahippocampal
    - DESeq_M_KAI_IA.R: Kainic acid intraamygdala
    - DESeq_M_KAI_IP.R: Kainic acid intraperitoneal

05.DESeq/02.DESeq_R_*.R
  Purpose: DESeq2 analysis for rat datasets
  Scripts for different comparisons:
    - DESeq_R_ALL.R: All rat comparisons
    - DESeq_R_PILO.R: Pilocarpine model
    - DESeq_R_KAI_SUB.R: Kainic acid subcutaneous
    - DESeq_R_KAI_IP.R: Kainic acid intraperitoneal
    - DESeq_R_PPS.R: Perforant path stimulation
    - DESeq_R_AMG.R: Amygdala kindling
    - DESeq_R_TBI.R: Traumatic brain injury

06. DECONVOLUTION
-----------------
06.Deconvolution/01.IMPORT_MIXTURE.R
  Purpose: Prepare TPM matrices for CIBERSORTx (Mouse/Rat)

06.Deconvolution/02.IMPORT_GO.R
  Purpose: Import gene sets for deconvolution

06.Deconvolution/03.CIBER_01234Docker.sh
  Purpose: Run CIBERSORTx for animal data

06.Deconvolution/04.CIBER_04HiRes_Merge.R
  Purpose: Merge high-resolution results

06.Deconvolution/05.DEG_Matrix.R
  Purpose: Prepare DEG matrices for cell types

06.Deconvolution/06.DEG_DESeq.R
  Purpose: DESeq2 for specific cell types

06.Deconvolution/07.DEG_DESeq_ALL.R
  Purpose: DESeq2 for all cell types and comparisons

07. ENRICHMENT ANALYSIS
-----------------------
07.Enrichment/ORA_01.Matrix.R
  Purpose: Prepare matrices for ORA

07.Enrichment/ORA_02.ToppFun.R
  Purpose: Perform ORA using ToppFun

07.Enrichment/ORA_03.ToppFun_category.R
  Purpose: ORA by gene set categories

07.Enrichment/GSEA_01.Matrix_GTF_H.R
  Purpose: Prepare GTF-based matrices for GSEA

07.Enrichment/GSEA_02.Matrix_GMT_H.R
  Purpose: Prepare GMT files for GSEA

07.Enrichment/GSEA_03.Matrix_TPM_M_H.R
  Purpose: Prepare mouse TPM matrix for GSEA

07.Enrichment/GSEA_03.Matrix_TPM_R_H.R
  Purpose: Prepare rat TPM matrix for GSEA

07.Enrichment/GSEA_04.NES_NG.R
  Purpose: Extract NES (Normalized Enrichment Score) for Neurogenesis

07.Enrichment/GSEA_04.NES_NT_IC.R
  Purpose: Extract NES for Neurotransmission and Ion Channels

07.Enrichment/GSEA_04.NES_NI.R
  Purpose: Extract NES for Neuroinflammation

08. CROSS-SPECIES GENE COMPARISON
----------------------------------
08.Cross_species_gene_comparison/01.Gene_level_comparison_final.R
  Purpose: Compare gene expression across species
  Functions:
    - Map orthologs between Human-Mouse-Rat
    - Extract log2FC and padj for each species
    - Identify shared and species-specific DEGs
    - Create comparison matrices

08.Cross_species_gene_comparison/02.Dotplot_DEGs.R
  Purpose: Visualize cross-species DEG comparisons
  Functions:
    - Create dot plots for shared genes
    - Visualize correlation patterns
    - Highlight species-specific changes


================================================================================
scRNA_HUMAN ANALYSIS PIPELINE
================================================================================

01. SINGLE-CELL DATA PROCESSING
--------------------------------
01.H_GSE160189_Seurat.R
  Purpose: Process GSE160189 single-cell dataset
  Functions:
    - Create Seurat object from count matrix
    - Quality control (already done)
    - Cell type annotation
    - Prepare for integration

01.H_GSE186538_Seurat.R
  Purpose: Process GSE186538 single-cell dataset
  Functions:
    - Create Seurat object
    - Quality control
    - Cell type annotation
    - Prepare for integration

02. INTEGRATION
---------------
02.Integ_2_RPCA.R
  Purpose: Integrate multiple scRNA-seq datasets using RPCA
  Functions:
    - Reciprocal PCA integration
    - Batch correction
    - Unified cell type annotation

03. DEconvolution SIGNATURE EXTRACTION
---------------------------------------
03.Integ_2_Deconv_Sig_DEG.R
  Purpose: Extract cell type-specific signatures for deconvolution
  Functions:
    - Identify marker genes for each cell type
    - Create signature matrices
    - Prepare for CIBERSORTx

04. VALIDATION
--------------
04.Integ_2_Validation_AvgEpr.R
  Purpose: Validate deconvolution using average expression
  Functions:
    - Compare bulk and single-cell expression
    - Calculate correlation
    - Assess deconvolution accuracy

05.Integ_2_Validation_Pseudobulk.R
  Purpose: Validate using pseudobulk approach
  Functions:
    - Create pseudobulk samples
    - Compare with actual bulk data
    - Validate cell type proportions


================================================================================
scRNA_ANIMAL ANALYSIS PIPELINE
================================================================================

01. SINGLE-CELL DECONVOLUTION
------------------------------
01.GSE185862_10X_Deconv.ipynb
  Purpose: Deconvolution analysis for GSE185862 10X dataset
  Functions:
    - Process 10X Genomics data
    - Cell type identification
    - Deconvolution signature extraction

02.GSE185862_10X_Deconv_OPC.ipynb
  Purpose: OPC (Oligodendrocyte Precursor Cell) specific analysis
  Functions:
    - Focus on OPC cell type
    - Extract OPC signatures
    - Deconvolution for OPC

03.GSE185862_10X_Deconv_plot.ipynb
  Purpose: Visualization of deconvolution results
  Functions:
    - Create plots for cell type proportions
    - Visualize spatial/temporal patterns

================================================================================
END OF DOCUMENTATION
================================================================================
