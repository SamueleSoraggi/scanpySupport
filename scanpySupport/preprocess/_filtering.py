def calculateQC(adata, show=True):
    """
    Function for calculating and plotting QC values. 
    It can be followed by filteringQC to apply filters.
    """
    from scipy.sparse import issparse
    import scanpy as sc
    import sklearn
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import pandas as pd
    import numpy as np
    
    sc.preprocessing.calculate_qc_metrics(adata, inplace=True)
        
    if issparse(adata.X):
        if 'chromosome' in adata.var_keys():
            adata.obs['percent_mito'] = np.sum( adata[:, adata.var['chromosome']=='MT'].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
        else:
            print("Percentage of MT genes not calculated (adata.var['chromosome'] not found)")
        adata.obs['scaled_CDR'] = sklearn.preprocessing.scale( np.sum(adata.X > 0, 1).A1 )
    else:
        if 'chromosome' in adata.var_keys():
            adata.obs['percent_mito'] = np.sum( adata[:, adata.var['chromosome']=='MT'].X, axis=1) / np.sum(adata.X, axis=1)
        else:
            print("Percentage of MT genes not calculated (adata.var['chromosome'] not found)")
        adata.obs['scaled_CDR'] = sklearn.preprocessing.scale( np.sum(adata.X > 0, 1) )    

    layers = adata.layers.keys()
    if ('spliced' in layers)and('unspliced' in layers):
        print("Calculating spliced and unspliced proportions")
        if (issparse(adata.layers['spliced']))and(issparse(adata.layers['unspliced'])):
            denom = np.sum(adata.layers['spliced']+adata.layers['unspliced'], 1).A1
            adata.obs['prop_spl'] = np.sum(adata.layers['spliced'], 1).A1 / denom
            adata.obs['prop_unspl'] = 1 - adata.obs['prop_spl']
        else:
            denom = np.sum(adata.layers['spliced']+adata.layers['unspliced'], 1)
            adata.obs['prop_spl'] = np.sum(adata.layers['spliced'], 1) / denom
            adata.obs['prop_unspl'] = 1 - adata.obs['prop_spl']
    else:
        print("Spliced and unspliced proportions not calculated ('spliced' and 'unspliced' layers not found)")
   
    return adata

def pca_outliers(adata, min_genes_per_cell=5, verbose=True):
    """
    Function to filter outliers using scater PCA on quality measures
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

    pandas2ri.activate()
    anndata2ri.activate()

    print("Loading objects into R")
    if issparse(adata.X):
        ro.globalenv['rawMatrix'] = adata.X.T.todense()
    else:
        ro.globalenv['rawMatrix'] = adata.X.T
    ro.globalenv['variables'] = adata.var_names.copy()
    ro.globalenv['observations'] = adata.obs[ ['total_counts'] ]
    
    print('Calculate PCA outliers')
    
    ro.r('')
    ro.r('pd <- DataFrame(data = observations)')
    ro.r('colnames(rawMatrix) <- rownames(pd)')
    ro.r('rownames(rawMatrix) <- variables')
    ro.r('sce <- SingleCellExperiment(assays = list(counts = as.matrix(rawMatrix) ), colData = pd)')
    ro.r('sce <- calculateQCMetrics(sce)')
    ro.r('sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE)')
    ro.r('cat("Nr of outliers detected:", sum(sce$outlier), sep=" ")')
    ro.r('outlier2 = sce@colData@rownames[sce$outlier]')
    ro.r('plotReducedDim(sce, use_dimred="PCA", shape_by = "outlier", size_by = "total_counts", colour_by = "total_features_by_counts")')

    outlier2 = ro.r('outlier2')
    adata = adata[np.invert(np.in1d(adata.obs_names, outlier2))].copy()
    sc.pp.filter_genes(adata, min_cells=min_genes_per_cell)

    return adata


def pca_covariates(adata, covariates=['total_counts'], verbose=False):
    """
    Function to output R^2 of covariates against PCA projection
    """
    import numpy as np
    import pandas as pd
    import rpy2.robjects as ro
    import anndata2ri
    import scanpy as sc
    from rpy2.robjects import pandas2ri
    from scipy.sparse import issparse
    import rpy2.rinterface_lib.callbacks
    import logging
    if not verbose:
        rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    import seaborn as sns
    import matplotlib.pyplot as plt

    ro.r('library(scater)')

    pandas2ri.activate()
    anndata2ri.activate()

    print("Loading objects into R")
    if issparse(adata.X):
        ro.globalenv['rawMatrix'] = np.log1p(adata.X.T.todense())
    else:
        ro.globalenv['rawMatrix'] = np.log1p(adata.X.T)
    ro.globalenv['observations'] = adata.obs[ covariates ]
    
    print('Calculate PCA covariates')
    
    ro.r('pd <- DataFrame(data = observations)')
    #ro.r('print(pd[1:5,])')
    ro.r('colnames(rawMatrix) <- rownames(pd)')
    ro.r('sce <- SingleCellExperiment(assays = list(counts = as.matrix(rawMatrix) ), colData = pd)')
    commandString = 'getVarianceExplained(sce, exprs_values = "counts", variables = c('
    variables  = ['"data.'+i+'"' for i in covariates]
    commandString = commandString + ','.join(variables) + ') )'
    print("using the R command")
    print(commandString)
    vals = ro.r(commandString)
    medians = np.argsort( -np.median(vals, 0) )
    medianVals = -np.sort( -np.median(vals, 0) )
    vals = pd.DataFrame( vals[:,medians] )
    #print(covariates)
    #print(medians)
    vals.columns = np.asarray(covariates)[medians]
    plt.rcParams['figure.figsize']=(8,8)
    f, ax = plt.subplots(1)
    for nn,mm in zip(vals.columns,medianVals):
        sns.kdeplot(vals[nn], ax=ax, label=nn, clip=(mm,97), gridsize=100)
    ax.set_xscale("symlog")
    #plt.xlim(0,100)
    ax.legend(title="Covariates",loc='best')

    adata.uns['pca_covariates'] = vals

    return adata
