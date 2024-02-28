# Wei_et_al_CRC_2024
The scripts for Wei et. al paper on cancer cell calling in colorectal cancer 

Processed data is available from [zenodo](https://dx.doi.org/10.5281/zenodo.10692019) 
and should be unpacked into a directory structure like so (under ./data/):

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

#### Preprocessing (Snakefile)
- Ambient RNA background removal by CellBender
    - input: {sample}_raw_feature_bc_matrix.h5
    - output: {sample}_cellbender_counts.h5

#### QC and UMAP (preprocessing_UMAP.ipynb)
- preprocessing_UMAP.ipynb
    - input: all cellbender_matrix.h5
    - output: CB_all_cells.h5 and CB_epi_cells_umap.h5
   
#### inferCNV (inferCNV.ipynb)
- 20230522_inferCNV.ipynb
    - input: CB_all_cells.h5 
    - output: all_cell_CB_counts.rds, all_epi_cell_CB_counts.rds, all_cell_anno.txt, all_epi_cell_anno.txt, inferCNV_expression_data.rds and other inferCNV intermediate objects
- 20230607_inferCNV_result.ipynb
    - collect inferCNV result
    - output: csv files
    
#### Copy number inference by Numbat (Snakefile)
- Since p009 have two normal samples and two tumour samples, they were run separately with modified scripts
    - combine_p009t1t2_and_run_numbat.ipynb
- collect Numbat result (collect_numbat_result.ipynb)

#### CCISM

- use `script/run_CCISM.sh`

#### iCMS 
- preprocessing of the count matrix in 20230530_preprocessing_Joanito.ipynb
- scVI and scANVI model training 

#### cell type annotation



#### Figures in the publication (Figures.ipynb)
- input: data_consolidated.h5ad
- output: figures/



Solving list
For Benedikt
- the h5ad object is missing from the file structure tree
- .git/ Permission denied

For cindy
- check what the zenodo repo actually looks like and adapt to that
- all hard-coded path
- name the script e.g. 1_, 2_xxx
    

