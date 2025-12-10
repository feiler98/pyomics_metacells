# imports
# ----------------------------------------------------------------------------------------------------------------------
from pyomics import utils as ut
from pathlib import Path
import pandas as pd
import scanpy as sc
from logging_utils import log_out
from split_hca_utils import adata_split_by_tissue, adata_filter_normal_cells
import urllib.request
# ----------------------------------------------------------------------------------------------------------------------

# pathlib.Path raw/in/out
path_raw = Path(__file__).parent / "raw_data"
print(path_raw)
path_in = Path(__file__).parent / "data"
path_out = Path(__file__).parent / "metacells_out"

# create raw_data and data if not exists
if not path_raw.exists():
    raise ValueError(f"Path {str(path_raw)} must exist to initiate the pipeline!")
if not path_out.exists():
    path_out.mkdir(exist_ok=True, parents=True)
if not path_in.exists():
    path_in.mkdir(exist_ok=True, parents=True)

# logging of stdin and stderr
# log_out(path_log=path_out, log_name=f"{Path(__file__).stem}__log")

# retrieve CellxGene Data as h5ad in raw_data folder
# ----------------------------------------------------------------------------------------------------------------------
print("""
###############################
# Download CellxGene datasets #
###############################
""")

"""
df_data_url = pd.read_csv(path_raw / "data_download_partial.txt", index_col="tag", sep=", ", engine="python")
for tags in list(df_data_url.index):
    url = df_data_url.loc[tags, "link"]
    retrive_as_path = path_raw / f"{tags}.h5ad"
    urllib.request.urlretrieve(url, retrive_as_path)
    print(f"Dataset < {url} > has been retrieved as {retrive_as_path}")
"""

# split data and save as .h5 in data
# ----------------------------------------------------------------------------------------------------------------------
print("""
################################
# HCA file splitting by tissue #
################################
""")

list_paths_h5ad = [p for p in path_raw.glob("*.h5ad")]
for p in list_paths_h5ad:
    import_adata = sc.read_h5ad(p, chunk_size=1000)
    adata_split_by_tissue(import_adata, p.stem, path_in)

# metacell generation using SeaCells
# ----------------------------------------------------------------------------------------------------------------------
print("""
####################################
# SEACells metacell generation log #
####################################
""")

# processing the preprepared files
list_paths_h5ad = [p for p in path_in.glob("*.h5")]
for p_h5 in list_paths_h5ad:
    header_string = f"Generation of metacells for < HCA non-neuronal cells | {p_h5.stem} >"
    print(header_string)
    print("#"*len(header_string))
    path_folder_seacells_out = path_out / f"{p_h5.stem}__seacells"
    path_folder_seacells_out.mkdir(exist_ok=True, parents=True)
    adata = sc.read_h5ad(p_h5)
    adata = adata_filter_normal_cells(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=30000, subset=True, flavor="seurat_v3_paper")
    ut.scrna_dim_red_by_metacells(adata=adata,
                                  data_tag=p_h5.stem,
                                  ad_obs_cell_type_assignment="cell_type",
                                  path_out=path_folder_seacells_out)


if __name__ == "__main__":
    pass
