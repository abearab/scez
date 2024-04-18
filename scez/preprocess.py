import scanpy as sc


def normalization(adata,lognorm=True,final_layer='log1p_norm'):
    # keep raw counts as a layer
    adata.layers['raw_counts'] = adata.X.copy()
    scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
    if lognorm:
        # log1p transform
        adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    adata.X = adata.layers[final_layer]


def clustering(
        adata, n_pcs=50, n_neighbors=30, use_highly_variable=True,
        use_rep=None, resolution=None
        ):
    if use_highly_variable:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
    else:
        sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)#, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)
