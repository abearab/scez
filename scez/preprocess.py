import scanpy as sc


def normalization(adata, lognorm=False):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    if lognorm:
        sc.pp.log1p(adata)


def clustering(adata, lognorm=False, rep=None, n_neighbor=30,n_comps=50,res=None):
    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=n_neighbor)
    # why can't we just work with the default neighbors?
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=res)
