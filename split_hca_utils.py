import anndata as ad
from pathlib import Path

import scanpy


def adata_filter_normal_cells(adata: ad.AnnData) -> ad.AnnData:
    list_normal_cells = [cell_idx for cell_idx, row in adata.obs.iterrows() if row["disease"] == "normal"]
    return adata[list_normal_cells, :]

def adata_split_by_tissue(adata: ad.AnnData, data_tag: str, out_path: (str, Path)):
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
    #adata_split_by_tissue(scanpy.read_h5ad("/home/feilerwe/coding_playground/test_adata_cellxgene/data.h5ad"), "test_data",Path.cwd())
    pass