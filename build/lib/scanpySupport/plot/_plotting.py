def plotQC(adata):
    """
    Function for plotting QC values. 
    Needs to run preprocess.calculateQC first.
    """
    from scipy.sparse import issparse
    import scanpy as sc
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import pandas as pd
    import numpy as np
    import seaborn as sns
    
    
    df=pd.DataFrame(index=adata.obs_names)
    df['total_counts'] = adata.obs['total_counts']
    df['n_genes_per_cell'] = adata.obs['n_genes_by_counts']
    if 'percent_mito' in adata.obs_keys():
        df['percent_mito'] = adata.obs['percent_mito']
    else:
        df['percent_mito'] = np.zeros(adata.shape[0])
        
    plt.rcParams['figure.figsize']=(20,20)
    f, ax = plt.subplots(3,2)

    sns.scatterplot(x='total_counts', y='n_genes_per_cell', hue='percent_mito', data=df, ax=ax[0,0])
    ax[0,0].set_title('UMI vs GENES plot - percent mito genes')

    if 'prop_spl' in adata.obs_keys():
        df['rate_spliced'] = adata.obs['prop_spl']
    else:
        df['rate_spliced'] = np.zeros(adata.shape[0])        
    sns.scatterplot(x='total_counts', y='n_genes_per_cell', hue='rate_spliced', data=df, ax=ax[0,1])
    ax[0,1].set_title('UMI vs GENES plot - spliced proportions')

    sns.distplot(df['total_counts'], ax=ax[1,0])
    ax[1,0].set_title('UMI counts per cell')

    sns.distplot(df['n_genes_per_cell'], ax=ax[1,1])
    ax[1,1].set_title('Genes per cell')

    df['counts_over_genes'] = df['total_counts']/df['n_genes_per_cell']
    sns.scatterplot(x='total_counts', y='counts_over_genes', hue='percent_mito', data=df, ax=ax[2,0])
    ax[2,0].set_title('UMI vs UMI/GENES ratio plot - percent MT genes') 
    sns.scatterplot(x='total_counts', y='counts_over_genes', hue='rate_spliced', data=df, ax=ax[2,1])
    ax[2,1].set_title('UMI vs UMI/GENES ratio plot - spliced proportions')

    plt.rcParams['figure.figsize']=(6,6)


def plotClusterOverlap(adata, marker_file="./Markers.txt", cluster_key='leiden', gene_rank_key='rank_genes_leiden', remove_somatic=True, normalize='None', rename=True, method='overlap_count'):

    import scanpy as sc
    import string
    import numpy as np
    import pandas as pd 
    import matplotlib.pyplot as plt
    import seaborn as sns

    with open(marker_file,'r') as f: #read names both lower and upper case
        names=[x.strip().split('\t')[0] for x in f]
    with open(marker_file,'r') as f:    
        types=[x.strip().split('\t')[1] for x in f]
    with open(marker_file,'r') as f:
        names2=[x.strip().split('\t')[0].upper() for x in f]    

    del names[0] #remove headers
    del names2[0]
    del types[0]

    typesNew = np.concatenate( (np.array(types), np.array(types)) ) #put vectors together and find which genes are in our data
    namesNew = np.concatenate( (np.array(names),np.array(names2)) )
    idx = np.in1d(namesNew, adata.var_names)
    typesNew = typesNew[idx]
    namesNew = namesNew[idx]

    forTypes = np.unique(typesNew)
    clusterDict = dict()
    for tt in forTypes:
        idx = np.in1d( typesNew, tt )
        clusterDict[str(tt)] = np.array(namesNew)[idx]

    if(remove_somatic): #remove somatic cells if you want to
        del clusterDict['Somatic']
        
    cell_annotation = sc.tl.marker_gene_overlap(adata, clusterDict, key=gene_rank_key, normalize=normalize, method=method)

    nameCluster = []
    for i in range(cell_annotation.max().size): #find cluster names from overlapping matrix
        nameCluster.append( np.array( cell_annotation.index[ cell_annotation[str(i)] == cell_annotation.max()[i]  ] ) )  

    for i in range(cell_annotation.max().size): #resolve duplicated names
        nameCluster[i] = '+'.join(nameCluster[i])

    nameCluster=np.array(nameCluster)
    newCluster=nameCluster
    for i in np.unique(nameCluster):
        idx = [j==i for j in nameCluster]
        if sum(idx)>1:
            numbers=range(1, (sum(idx))+1 )
            newCluster[idx]= ["{}-{}".format(a,b) for a,b in zip( nameCluster[idx], numbers )]

    plt.rcParams['figure.figsize']=(15,15)
    sns.heatmap(cell_annotation, cbar=False, annot=True )
    plt.rcParams['figure.figsize']=(6,6)
