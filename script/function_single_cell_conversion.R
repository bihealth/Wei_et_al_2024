library(Seurat)
library(rhdf5)

# read scanpy generated h5 object

read_scanpy_h5 = function(path, if_pca=TRUE, if_umap=TRUE, if_raw_counts=TRUE, if_CB_counts=TRUE){
    
    #read scanpy generated h5 file and convert to Seurat object 

    # read in file
    adata = H5Fopen(path)    
    
    # X for Seurat
    feature_matrix <- Matrix::sparseMatrix(
        i = c(adata$X$indices),
        p = c(adata$X$indptr),
        x = c(adata$X$data),
        dimnames = list(
            adata$var$`_index`, #feature name
            adata$obs$`_index` #barcode
        ),
        index1 = FALSE
    )
    
    # convert to Seurat object
    sobj = CreateSeuratObject(counts = feature_matrix, project='sccrc')
    DefaultAssay(sobj) <- "RNA"
    
    if(if_raw_counts == TRUE){
    # raw count matrix for Seurat
        raw_matrix <- Matrix::sparseMatrix(
            i = c(adata$layers$counts$indices),
            p = c(adata$layers$counts$indptr),
            x = c(adata$layers$counts$data),
            dimnames = list(
                adata$var$`_index`, #feature name
                adata$obs$`_index` #barcode
            ),
            index1 = FALSE)
        
        raw_count_assay <- CreateAssayObject(counts = raw_matrix)
        sobj[["raw_counts"]] <- raw_count_assay
    
    } else {
    }
    
    
    if(if_CB_counts == TRUE){
    # raw or CellBender count matrix for Seurat
        raw_matrix <- Matrix::sparseMatrix(
            i = c(adata$layers$CB_counts$indices),
            p = c(adata$layers$CB_counts$indptr),
            x = c(adata$layers$CB_counts$data),
            dimnames = list(
                adata$var$`_index`, #feature name
                adata$obs$`_index` #barcode
            ),
            index1 = FALSE)
        
        raw_count_assay <- CreateAssayObject(counts = raw_matrix)
        sobj[["raw_counts"]] <- raw_count_assay
    
    } else {
    }
    
    if(if_pca == TRUE){
        
        # embedding
        pca_mx = as.matrix(t(adata$obsm$X_pca))
        colnames(pca_mx)=paste0('PC_', 1:dim(pca_mx)[2])
        rownames(pca_mx)=adata$obs$`_index`
        
        # loadings
        pca_loadings = as.matrix(t(adata$varm$PCs))
        colnames(pca_loadings)=paste0('PC_', 1:dim(pca_mx)[2])
        rownames(pca_loadings)=adata$var$`_index`
        
        # add to Seurat object
        sobj[['pca']] = CreateDimReducObject(embeddings = pca_mx,
                                             loadings = pca_loadings,
                                             key = "PC_", 
                                             assay = DefaultAssay(sobj))
    
    } else {
    }
    
    if(if_umap == TRUE){
        
        # embedding # assume dim=2
        umap_mx = as.matrix(t(adata$obsm$X_umap))
        rownames(umap_mx)=adata$obs$`_index`
        colnames(umap_mx)=paste0('UMAP_', 1:2)
        
        sobj[['umap']] = CreateDimReducObject(embeddings = umap_mx,
                                              key = "UMAP_", 
                                              assay = DefaultAssay(sobj))
    
    } else {
    }
    
    return(sobj)
    h5closeAll()
}



add_categorical_obs_to_sobj = function(h5fopen_obj, cat_name, sobj){
    
    adata = h5fopen_obj
 
    cat = c(adata$obs$cat_name$categories) # the dollar cannot take the variable

    # fix 0-index in py vs 1-index in R
    cat_obs = data.frame(cat_name = cat[apply(adata$obs$cat_name$codes, 1, function(x)x+1)],
                                              row.names = adata$obs$`_index`)

    sobj = AddMetaData(object = sobj,
                       metadata = cat_obs, col.name = cat_name)
    
    return(sobj)
                 

}
                                        

get_categorical_obs = function(h5fopen_obj, cat_name){
    
    adata = h5fopen_obj
 
    cat = c(adata$obs$cat_name$categories)

    # fix 0-index in py vs 1-index in R
    cat_obs = data.frame(cat_name = cat[apply(adata$obs$cat_name$codes, 1, function(x)x+1)],
                                              row.names = adata$obs$`_index`)
    return(cat_obs)                  

}





