import anndata as ad
from pathlib import Path


def adata_split_by_tissue(adata: ad.AnnData, data_tag: str, out_path: (str, Path)):
    out_path.mkdir(exist_ok=True, parents=True)
    list_available_tissues = list(set(adata.obs["tissue"]))
    dict_sort = {}
    for tissue in list_available_tissues:
        tissue_df_obs =  adata.obs.where(adata.obs["tissue"] == tissue).dropna()
        dict_tissue = {}
        for key, row in tissue_df_obs.iterrows():
            new_key = f"{key.split(":")[1]}__{row["cell_type"]}__{row["disease"]}_{row["sex"]}_{row["development_stage"].split("-")[0]}y/o"
            dict_tissue.update({key:new_key})
        dict_sort.update({tissue:dict_tissue})
    print(f"The data can be split in the following tissues: {" ".join(list(dict_sort.keys()))}")
    for key, dict_tissue in dict_sort.items():
        slice_adata = adata[adata.obs.tissue == key]
        slice_adata.obs.rename(index=dict_tissue, inplace=True)
        slice_adata.write_h5ad(out_path  / f"{key}__{data_tag}.h5")
        print()



if __name__ == "__main__":
    pass
