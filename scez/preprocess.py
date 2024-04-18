import scanpy as sc


def normalization(adata,lognorm=True):
    # keep raw counts as a layer
    adata.layers['raw_counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    if lognorm:
        sc.pp.log1p(adata)


def clustering(
        adata, n_pcs=50, n_neighbor=30, use_highly_variable=True,
        use_rep=None, resolution=None
        ):
    if use_highly_variable:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    else:
        sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbor, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)
