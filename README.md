# Wei_et_al_CRC_2024
The scripts for Wei et. al paper on cancer cell calling in colorectal cancer 

Processed data is available from [zenodo](https://dx.doi.org/10.5281/zenodo.10692019) 
and should be unpacked into a directory structure like so (named ./data/):

```
.
├── cellbender
│   ├── p007n_cellbender_counts.h5
│   ├── p007t_cellbender_counts.h5
│   ├── p008n_cellbender_counts.h5
│   └── ...
├── cellranger
│   ├── p007n_raw_feature_bc_matrix.h5
│   ├── p007t_raw_feature_bc_matrix.h5
│   ├── p008n_raw_feature_bc_matrix.h5
│   └── ...
├── cellSNP
│   ├── downsampled
│   │   ├── p007t_05
│   │   │   ├── cellSNP.base.vcf
│   │   │   ├── cellSNP.samples.tsv
│   │   │   ├── cellSNP.tag.AD.mtx
│   │   │   ├── cellSNP.tag.DP.mtx
│   │   │   └── cellSNP.tag.OTH.mtx
│   │   ├── ...
│   ├── p007n
│   │   ├── cellSNP.base.vcf
│   │   ├── cellSNP.cells.vcf
│   │   ├── cellSNP.samples.tsv
│   │   ├── cellSNP.tag.AD.mtx
│   │   ├── cellSNP.tag.DP.mtx
│   │   └── cellSNP.tag.OTH.mtx
│   ├── ...
├── numbat
│   ├── p007_clone_post_2.tsv
│   ├── p008_clone_post_2.tsv
│   └── ...
└── WGS
    ├── p007_filtered.vcf.gz
    ├── p008_filtered.vcf.gz
    └── ...
```


## Analysis pipeline

### Ambient RNA background removal (Snakefile)
- Ambient RNA background removal by CellBender
    - input: {sample}_raw_feature_bc_matrix.h5
    - output: {sample}_cellbender_counts.h5

### Preprocessing
#### Integrate all samples from CellBender output and filter for QC parameters
#### Cell type annotation (coarse, epithelial cells vs immune cells vs stromal cells)
- run `1_preprocessing_h5.ipynb`
    - input: all cellbender_counts.h5 object
    - output: 
        - CB_all_cells.h5 and CB_epi_cells.h5 that contains high-quality all/epithelial cells with calculated PCA, UMAP, diffusion map, and louvain embeddings, as well as coarse cell type annotation
        - anno/{sample}_{cell_type}.txt cell barcode list by sample and cell type

### Identification of cancer cells
#### inferCNV: infer copy number alteration status from gene expression
- run `inferCNV/inferCNV.ipynb` to execute inferCNV 
    - input: CB_all_cells.h5 
    - output: all_cell_CB_counts.rds, all_epi_cell_CB_counts.rds, all_cell_anno.txt, all_epi_cell_anno.txt, inferCNV_expression_data.rds and other inferCNV intermediate objects
    
- run `inferCNV/inferCNV_result.ipynb` to collect inferCNV result
    - collect inferCNV result from the previous run
    - output: `infercnv_clone_scores.tsv`
    
#### Numbat: infer copy number alterations from phased gene expression profiles
- run `script/Snakefile` for all the samples except p009
- run `script//Numbat/combine_p009t1t2_and_run_numbat.ipynb` since p009 have two normal samples and two tumour samples, they were run separately with modified scripts
- run `script/Numbat/collect_numbat_result.ipynb` to collect Numbat result of all the samples
    - output: `numbat_all_output_clone_post_combined9.csv`

#### CCISM: identify cancer cells using SNVs
- rub `script/run_CCISM.sh` to run CCISM on all the sample

#### iCMS: annotates cancer cell phenotypes (iCMS2/iCMS3)
- input
    - download the h5 object from Joanito et al.: 
    - `CB_epi_cells.h5`
- run `iCMS/preprocessing_Joanito.ipynb` to filter the count matrix with the same criteria as our h5 objects
    - output `adata_concat_with_joanito.h5`
- run `iCMS/run_scvi_snakemake.sh` for scVI model training 
- run `iCMS/scvi_model_result_iCMS.ipynb` to inspect scVI modeling result, train and inspect scANVI models
    - output `CB_epi_cells_iCMS.h5

### Integrate results from inferCNV, Numbat, CCISM, and iCMS
- run `2_integrate_tools_result.ipynb`
- input
    - `CB_epi_cells_iCMS.h5` which contains iCMS result
    - CCISM result
    - Numbat: `numbat_all_output_clone_post_combined9.csv`
    - inferCNV: `infercnv_clone_scores.tsv`
- output: `CB_epi_Numbat_CCISM_inferCNV_iCMS.h5`

### Cell type annotation at finer resolution
- run `3_cell_type_annotation_finer.ipynb`
    - output `adata_all_full_cell_type_annotation.h5`

### Consensus cancer calls
- run `4_consensus_calls.ipynb`
    - output `CB_epi_Numbat_CCISM_inferCNV_icms_Uhlitz_resolved_identity.h5`

### Pseudotime and cell type enrichment
- run `5_cellrank.ipynb`

### Ligand and receptor expression
- run `6a_DEG_epi.ipynb` for expression and differential expression analysis of epithelial cells; 
`6b_DEG_imm_str.ipynb` for immune and stromal cells

### CRC pathways and signatures expression levels
- run `7_CRC_pathway_signature.ipynb`

### Cell-cell interaction
- run `8_cellchat.ipynb` to infer cell-cell interactions in the normal versus tumour samples



## Figures in the publication (Figures.ipynb)
- run `script/Figures.ipynb`
    - input: `data/data_consolidated.h5ad`
    - output: `figures/*`
- Figure 2 contains data simulation result, run `script/Fig2_simulations.Rmd`
    - output `Fig2_simulation_results.csv`



