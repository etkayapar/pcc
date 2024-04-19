import pandas as pd

gene_table = pd.read_table("gene_table.tsv", index_col=False, dtype=str)
genes = gene_table["busco_id"]

