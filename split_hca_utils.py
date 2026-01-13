import anndata as ad
from pathlib import Path
from scipy import shuffle

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

        # performance issues with SEACells now force the splitting if datasets too large
        # threshold of splitting is set to 70,000 cells --> split by 2
        n_max_obs_fit = int(len(slice_adata.obs) / 70000)
        if n_max_obs_fit == 0:
            list_split_slices = [slice_adata]
        else:
            shuffle(slice_adata.obs, inplace=True)
            cell_tag_list_obs = list(slice_adata.obs.index)
            list_split = list(range(0, len(slice_adata.obs), int(len(slice_adata.obs)/n_max_obs_fit)))
            list_split.append(len(slice_adata.obs))
            i = 1
            list_split_slices = []
            while i < len(list_split):
                list_split_slices.append(slice_adata[cell_tag_list_obs[i-1:i], :])

        for i, anndata_slice in enumerate(list_split_slices):
            if len(list_split_slices) > 1:
                tag_save = f"{key}__{data_tag}__split_{i}.h5"
            else:
                tag_save = f"{key}__{data_tag}.h5"
            slice_adata.write(out_path  / tag_save)
            print(f"Adata slice is saved as {tag_save} in {out_path}")



if __name__ == "__main__":
    pass