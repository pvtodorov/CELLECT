import pandas as pd
import numpy as np

output_dir = snakemake.params['out_dir']
specificity_matrix_name = snakemake.params['specificity_matrix_name']
specificity_matrix_file = snakemake.params['specificity_matrix_file']


specificity_df = pd.read_csv(specificity_matrix_file, index_col=None)

all_genes = pd.DataFrame(index = spec_matrix_df.index)
all_genes['all_genes_in_dataset'] = 1

all_genes_path = '{out_dir}/all_genes.multi_geneset.{sm_name}.txt'.format(
					  out_dir = output_dir,
					  sm_name = specificity_matrix_name)
all_genes.to_csv(mgs_all_genes_path, sep='\t')