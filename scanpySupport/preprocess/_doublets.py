import scrublet as scr

def scrublet(adata, expected_rate=0.06, doublet_score=None, show=True):
    #import scrublet as scr
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_rate)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    if(show):
        scrub.plot_histogram()
    
    adata.obs['doublet_scores']=doublet_scores
    adata.obs['predicted_doublets']=predicted_doublets
    
    adata.obs['predicted_doublets'] = adata.obs['doublet_scores']>doublet_score  
