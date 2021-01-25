#import rpy2.rinterface_lib.callbacks
#import logging
#from rpy2.robjects import pandas2ri
#import anndata2ri


def sctransform(adata,
                genes=2000,
                min_genes_per_cell=5,
                method='poisson',
                latent=None,
                batch=None,
                cores=1,
                memory=10,
                verbose=True):
    """
    Function to use scTransform. It needs at least the adata.obj['total_counts'] number of UMIs calculated in the data.
    """
    import numpy as np
    import rpy2.robjects as ro
    import anndata2ri
    import scanpy as sc
    from rpy2.robjects import pandas2ri
    from scipy.sparse import issparse
    import rpy2.rinterface_lib.callbacks
    import logging
    if not verbose:
        rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    
    ro.r('library(scater)')
    ro.r('library(sctransform)')
    ro.r('library(future)')
    pandas2ri.activate()
    anndata2ri.activate()

    print('Filtering genes')
    sc.pp.filter_genes(adata, min_cells=min_genes_per_cell)
    
    if issparse(adata.X):
        ro.globalenv['rawMatrix'] = adata.X.T.todense()
    else:
        ro.globalenv['rawMatrix'] = adata.X.T

    latent_var = []

    if latent is None:
        ro.r('cells_info = as.data.frame( colSums(rawMatrix) )')
        ro.globalenv['cellnames']=np.asarray(adata.obs_names)
        ro.r('rownames(cells_info) = cellnames')
    else:
        latent_var=latent
        ro.globalenv['cells_info'] = adata.obs[ latent_var ]
        latent_var = ['"data.'+i+'"' for i in latent_var]
    ro.globalenv['genes_name'] = adata.var_names
    
    ro.r('cell_df <- DataFrame(data = cells_info)')
    #ro.r('print(head(cell_df))')
    #ro.r('print(rownames(cell_df)[1:10])')
    #ro.r('rawMatrix=as.data.frame(rawMatrix)')
    ro.r('colnames(rawMatrix) <- rownames(cell_df)')
    ro.r('rownames(rawMatrix) <- genes_name')
    print('Configure future multithreading')
    ro.globalenv['cores']=cores
    ro.globalenv['memory']=memory
    ro.r('future::plan(strategy = \'multicore\', workers = cores)')
    ro.r('options(future.globals.maxSize = memory * 1024 ^ 3)')
    print('Run scTransform')
    ro.globalenv['genes']=int(genes)
    ro.globalenv['min_genes_per_cell']=int(min_genes_per_cell)
    ro.globalenv['method']=method
    stringCommand = 'vst_out=vst( as.matrix(rawMatrix), cell_attr=cell_df, n_genes=genes, method=method, show_progress=TRUE, min_cells=min_genes_per_cell, return_corrected_umi=TRUE'
    #latent_var = ['"data.'+i+'"' for i in latent_var]
    if batch is not None:
        batch = '"data.'+batch+'"'
        stringCommand = stringCommand + ', batch_var=' + batch
        if latent is not None:
            latent_var.remove(batch)
    if ((len(latent_var)>1) and (batch is not None))|((len(latent_var)>=1) and (batch is None)):
        #print(latent_var)
        stringCommand = stringCommand + ', latent_var=c(' + ','.join(latent_var) + ')'
    stringCommand += ')'
    print("Running the command:",stringCommand)
    ro.r(stringCommand)
    print('Extract results')
    new_matrix= ro.r('vst_out$y')
    sct_genes = ro.r('rownames(vst_out$model_pars)')
    all_genes = ro.r('rownames(vst_out$y)')
    umi_corrected = ro.r('vst_out$umi_corrected')

    adata = adata[:,all_genes].copy()
    adata.var['highly_variable'] = [i in sct_genes for i in adata.var_names]
    adata.layers['norm_sct'] = np.transpose( new_matrix )
    adata.layers['umi_corr'] = umi_corrected.T.copy()
    
    return adata
