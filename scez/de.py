import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import anndata as ad

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from adpbulk import ADPBulk

inference = DefaultInference(n_cpus=8)


def pseudobulk_by_clusters(adt, condition, cluster_col='leiden', method="mean"):
    # initialize the object
    adpb = ADPBulk(adt, [cluster_col, condition], method=method)

    # perform the pseudobulking
    pseudobulk_matrix = adpb.fit_transform()

    # retrieve the sample metadata (useful for easy incorporation with edgeR)
    sample_meta = adpb.get_meta()

    out = ad.AnnData(
        X=pseudobulk_matrix,
        obs=sample_meta.set_index('SampleName')
    )

    return out


def run_deseq(adata, design):
    dds = DeseqDataSet(
        counts=adata.to_df().astype(int),
        metadata=adata.obs,
        design_factors=design,  # compare samples based on the "condition"
        # column ("B" vs "A")
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    stat_res = DeseqStats(dds, inference=inference)
    stat_res.summary()

    df = stat_res.results_df

    return df


def plot_volcano(df, title=None, font_scale=1):
    df['name'] = df.index.to_list()

    df['-log10(pvalue)'] = - np.log10(df.pvalue)

    fig, ax = plt.subplots(figsize=(3, 3))

    # Scatter plot
    ax.scatter(
        df['log2FoldChange'],
        df['-log10(pvalue)'],
        alpha=0.9, s=5, c='#fcae91'
    )

    # Set background color to transparent
    ax.set_facecolor('none')

    # Set smaller font size
    ax.tick_params(axis='both', which='both', labelsize=8 * font_scale)

    # Set labels
    ax.set_xlabel('log2FoldChange', fontsize=9 * font_scale)
    ax.set_ylabel('-log10(pvalue)', fontsize=9 * font_scale)

    # Set plot title
    if not title:
        ax.set_title('Volcano Plot', fontsize=10 * font_scale)
    else:
        ax.set_title(title, fontsize=10 * font_scale)

    ax.grid(False)

    top_genes = df.nlargest(4, '-log10(pvalue)')  # Adjust the number as needed
    for index, row in top_genes.iterrows():
        ax.annotate(row['name'], (row['log2FoldChange'], row['-log10(pvalue)']), fontsize=5 * font_scale, ha='right',
                    va='bottom')

    plt.tight_layout()
    plt.show()


def plot_top_DEG_violinplot(adata, df):
    top_genes = df.nlargest(20, '-log10(pvalue)')  # Adjust the number as needed

    # Filter the single-cell dataset for the top genes
    subset_adata = adata[:, top_genes.index]
    subset_adata.var.index = subset_adata.var.index.str.split('_').str[0]

    # Convert the subset of adata to a DataFrame
    subset_df = subset_adata.to_df()

    # Merge the DataFrame with .obs to include the 'sample' information
    merged_df = pd.merge(subset_df, adata.obs[['sample']], left_index=True, right_index=True)

    # Melt the DataFrame to prepare for violin plot
    melted_df = pd.melt(merged_df, id_vars='sample', var_name='Gene', value_name='Counts')

    # Create a violin plot
    plt.figure(figsize=(10, 4))
    sns.violinplot(x='Gene', y='Counts', hue='sample', data=melted_df, split=True, inner='quartile', palette='Set2')
    sns.stripplot(x='Gene', y='Counts', hue='sample', data=melted_df, dodge=True, jitter=True, color='black', size=1,
                  alpha=0.3)

    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.title('Top Differentially Expressed Genes')
    plt.show()


def write_top_DEGs(df, sample_id, n_hits=200):
    df['-log10(pvalue)'] = - np.log10(df.pvalue)
    df.nlargest(n_hits, '-log10(pvalue)').to_csv(f'{sample_id}_top_{n_hits}.csv')  # Adjust the number as needed
