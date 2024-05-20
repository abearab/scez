import scanpy as sc


def normalization(adata, target_sum=1e4, max_value=10, final_layer='scaled', keep_initial_layer=True):
    if keep_initial_layer == True:
        adata.layers['raw_counts'] = adata.X.copy()
    elif type(keep_initial_layer) == str:
        adata.layers[keep_initial_layer] = adata.X.copy()
    
    # normalize counts to target_sum (default 1e4)
    counts = sc.pp.normalize_total(adata, target_sum=target_sum, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(counts["X"], copy=True)
    # scale counts
    adata.layers['scaled'] = sc.pp.scale(adata, max_value=max_value, copy=True).X
    # set the final layer
    adata.X = adata.layers[final_layer]


def clustering(
        adata
        ):
    pass
    # , n_pcs=50, n_neighbors=30, use_highly_variable='Yes',
    #     use_rep=None, resolution=None
    
    # if use_highly_variable == 'Yes':
    #     sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #     sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
    # else:
    #     sc.pp.pca(adata, n_comps=n_pcs)
    # sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)#, n_pcs=n_pcs)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, resolution=resolution)
