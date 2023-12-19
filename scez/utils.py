import pandas as pd


def add_marker_feature(adata, marker, marker_name, clusters_name, thr = 0, figsize=(10, 4)):

    adata.obs[marker_name] = ''
    adata.obs.loc[adata.to_df().loc[:,marker] <= thr, marker_name] = f'{marker}-'
    adata.obs.loc[adata.to_df().loc[:,marker] > thr, marker_name] = f'{marker}+'

    df = pd.concat([
        adata.obs.groupby([marker_name,clusters_name]).size()[f'{marker}+'],
        adata.obs.groupby([marker_name,clusters_name]).size()[f'{marker}-']
    ],axis=1).rename(columns={0:f'{marker}+',1:f'{marker}-'})

    # Make some labels.
    labels = df[f'{marker}+'] / df.sum(axis=1) * 100
    labels = labels.round(decimals=1)
    labels.sort_values(ascending=False,inplace=True)
    df = df.loc[labels.index,]

    ax = df.plot.bar(stacked=True,rot=0,figsize=figsize)

    rects = ax.patches

    for rect, label in zip(rects, labels):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2, height + 5, str(label) + "%",
            ha="center", va="bottom", fontsize=8
        )

    ax.set_yscale('log')
    ax.set_ylabel('# of cells')
    return ax