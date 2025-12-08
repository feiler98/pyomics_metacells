# imports
# ----------------------------------------------------------------------------------------------------------------------
from pyomics import utils as ut
from pathlib import Path
import scanpy as sc
from logging_utils import log_out
from split_hca_utils import adata_split_by_tissue
# ----------------------------------------------------------------------------------------------------------------------

# pathlib.Path raw/in/out
path_raw = Path(__file__) / "raw_data"
path_in = Path(__file__) / "data"
path_out = Path(__file__) / "metacells_out"

# create raw_data and data if not exists
if not path_raw.exists():
    raise ValueError(f"Path {str(path_raw)} must exist to initiate the pipeline!")
if not path_out.exists():
    path_out.mkdir(exist_ok=True, parents=True)
if not path_in.exists():
    path_in.mkdir(exist_ok=True, parents=True)

# logging of stdin and stderr
log_out(path_log=path_out, log_name=f"{Path(__file__).stem}__log")

print("""
################################
# HCA file splitting by tissue #
################################
""")

list_paths_h5ad = [p for p in path_raw.glob("*.h5ad")]
for p in list_paths_h5ad:
    import_adata = sc.read_h5ad(p, chunk_size=1000)
    adata_split_by_tissue(import_adata, p.stem, path_in)

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
    path_folder_seacells_out = path_out / f"non_neuro_{p_h5.stem}"
    path_folder_seacells_out.mkdir(exist_ok=True, parents=True)
    adata = sc.read_h5ad(p_h5)
    sc.pp.highly_variable_genes(adata, n_top_genes=20000, subset=True, flavor="seurat_v3_paper")
    ut.scrna_dim_red_by_metacells(adata=adata,
                                  data_tag=f"non_neuro_{p_h5.stem}",
                                  ad_obs_cell_type_assignment="cell_type",
                                  path_out=path_folder_seacells_out)

if __name__ == "__main__":
    pass
