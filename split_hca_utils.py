import anndata as ad
from pathlib import Path


def adata_filter_normal_cells(adata: ad.AnnData) -> ad.AnnData:
    """
    Removes all cells, which have for obs["disease"] != normal. Information is a requirement for the function to work
    obviously...

    Parameters
    ----------
    adata: ad.AnnData
        AnnData object from cellxgene.

    Returns
    -------
    ad.AnnData
        Sliced AnnData object.
    """

    list_normal_cells = [cell_idx for cell_idx, row in adata.obs.iterrows() if row["disease"] == "normal"]
    return adata[list_normal_cells, :]

def adata_split_by_tissue(adata: ad.AnnData, data_tag: str, out_path: (str, Path)):
    """
    Splits the AnnData object by the obs["tissue"] column and saves the data as .h5 file.
    Easier computation downstream with tissue agnostic metacell generation

    Parameters
    ----------
    adata: ad.AnnData
        AnnData object from cellxgene.
    data_tag: str
        Name for the export data.
    out_path: str|Path
        Directory for exporting the data.
    """

    Path(out_path).mkdir(exist_ok=True, parents=True)
    list_available_tissues = list(set(adata.obs["tissue"]))
    dict_sort = {}
    for tissue in list_available_tissues:
        tissue_df_obs =  adata.obs.where(adata.obs["tissue"] == tissue).dropna()
        dict_tissue = {}
        for key, row in tissue_df_obs.iterrows():
            # just use the standard key to prevent errors
            new_key = f"{key}__{row["cell_type"]}__{row["disease"]}_{row["sex"]}_{row["development_stage"].split("-")[0]}y/o"
            dict_tissue.update({key:new_key})
        dict_sort.update({tissue:dict_tissue})
    print(f"The data can be split in the following tissues: {" ".join(list(dict_sort.keys()))}")
    for key, dict_tissue in dict_sort.items():
        slice_adata = adata[adata.obs.tissue == key]
        slice_adata.obs.rename(index=dict_tissue, inplace=True)
        slice_adata.write(out_path  / f"{key}__{data_tag}.h5")
        print(f"Adata slice is saved as {key}__{data_tag}.h5 in {out_path}")



if __name__ == "__main__":
    pass