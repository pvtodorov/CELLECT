import pandas as pd
import numpy as np
import os
from time import gmtime, strftime

#TODO: missing file file

def read_SNP_gene_file(SNPs_to_genes_dir, chromosome, input_bim_path):
	'''
	Reads in the SNPs and genes file for a given chromosome
	
	Args:
		SNPs_to_genes_dir (str): The path to the directory containing the SNPs_with_genes files
		chromosome (int): A chromosome number
		input_bim_path (str): The path to the file with the input SNPs
	Returns:
		The SNPs with genes in a long format for the specified chromosome

	'''
	input_bim_filename = os.path.basename(input_bim_path)
	SNPs_to_genes_path = os.path.join(SNPs_to_genes_dir,
		'SNPs_with_genes.{bim_name}.{chromosome}.txt'.format(bim_name = input_bim_filename,
															chromosome = chromosome))
	
	SNPs_to_genes = pd.read_csv(SNPs_to_genes_path, sep='\t', index_col=False)
	# Splitting the rows on the gene column so that each row contains a single SNP and a single gene
	# From https://stackoverflow.com/questions/12680754/ UPDATE: generic vectorized approach
	lst_col = 'ID_genes_in_matched_locus'
	SNP_to_gene_tmp = SNPs_to_genes.assign(**{lst_col:SNPs_to_genes[lst_col].str.split(';')})
	# SNP_to_gene_tmp.dropna(subset=[lst_col], inplace=True)
	SNP_to_gene_tmp.loc[SNP_to_gene_tmp.loc[:,lst_col].isnull(),lst_col] = SNP_to_gene_tmp.loc[SNP_to_gene_tmp.loc[:,lst_col].isnull(),lst_col].apply(lambda x: [np.nan])
	SNP_to_gene_long = pd.DataFrame({col:np.repeat(SNP_to_gene_tmp[col].values,
									 SNP_to_gene_tmp[lst_col].str.len())
									 for col in SNP_to_gene_tmp.columns.difference([lst_col])}).assign(
	**{lst_col:np.concatenate(SNP_to_gene_tmp[lst_col].values)})[SNP_to_gene_tmp.columns.tolist()]
	
	SNP_to_gene_long.rename(index=str, columns={lst_col: "gene"}, inplace = True)

	return(SNP_to_gene_long)



def create_annot_file_per_chromosome(chromosome, run_prefix, precomp_dir, bim_path, all_genes):
	'''
	Merges the dataframe of SNPs associated to genes with the dataframe of expression specificity
	for a given gene to all cell types. 
	Then calculates a given function over each SNP and cell type in the merged 
	dataframe to get a final annotation matrix to be used as input to LD score regression.
	Args:
		chromosome (int): A chromosome number
	Returns:
		Nothing
	'''
	input_dir = os.path.join(precomp_dir, 'SNPsnap')
	snps_and_genes = read_SNP_gene_file(input_dir, chromosome, bim_path)

	# Joining the potential annotation values to the SNPs
	print('Combining the two input files for chromosome' + str(chromosome))
	combined_ES2SNP_Gene2ES_df = snps_and_genes.join(genes_and_ES,
		on='gene', lsuffix='_in', rsuffix='_snap', how='left')
	# Taking the max (for now, will integrate sum, min, etc later) over each
	# column, for each SNP
	print('Calculating annotation file for chromosome' + str(chromosome))
	max_all_annots = combined_ES2SNP_Gene2ES_df.groupby('snpID').max()
	max_all_annots.sort_values(by='BP coordinate', inplace=True, ascending=True)
	annot_names = genes_and_ES.columns
	if all_genes == True:
		out_chr_filename = 'all_genes_in_{prefix}.{chromosome}.annot.gz'.format(prefix=run_prefix,
																				chromosome=chromosome)
		run_prefix = "control.all_genes_in_dataset"
	else:
		out_chr_filename = '{prefix}.COMBINED_ANNOT.{chromosome}.annot.gz'.format(prefix=run_prefix,
																				  chromosome=chromosome)
	max_all_annots[annot_names].to_csv(os.path.join(precomp_dir,
													run_prefix,
													out_chr_filename),
									sep='\t', compression ='gzip', na_rep='0', index=False)

	print('Finished making SNP to gene file for chromosome' + str(chromosome))


###################################### MAIN ######################################

genes_and_ES = pd.read_csv((snakemake.input[0], index_col='gene')
chromosome = snakemake.params['chromosome']
run_prefix = snakemake.params['run_prefix']
precomp_dir = snakemake.params['precomp_dir']
all_genes = snakemake.params['all_genes']
bfile_path = snakemake.params['bfile']

create_annot_file_per_chromosome(chromosome, run_prefix, precomp_dir, bfile_path, all_genes)
	
