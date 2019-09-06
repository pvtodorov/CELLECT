import pandas as pd
import os
import re

"""
@TODO:
    
    - Comment the code
    - implement h2 so that it fits the new script

Description:
    We directly use the declared files from listtargetfiles
    We find them and then we concatenate them.
    We also expand the csv with metadata on which specificity_id, 
    which gwas they're regressed on and etc 
"""

#base_output_dir = snakemake.params['BASE_OUTPUT_DIR']
results_output_dir = snakemake.params['results_out_dir'] 
#run_prefix = snakemake.params['run_prefix']
#analysis_type = snakemake.params['analysis_type']
#gwas = snakemake.params['gwas']
# file

targets = snakemake.input

prioritization_combined = pd.DataFrame()
conditional_combined = pd.DataFrame()
h2_combined = pd.DataFrame()
#analysis_type_dir = os.path.join(base_output_dir, analysis_type)

for file in targets:
    #print(file)
    
    if '/prioritization/' in file:
        prioritization = pd.read_csv(file, sep = '\t', header = 0)

        file = file.split('/')
        metadata = file[-1]
        metadata = re.sub(r'\.cell_type_results.txt', '', metadata)
        metadata = metadata.split('__')
        #dict_meta_data = {metadata[0]: metadata[1]}
        #metadata_df = pd.DataFrame(, columns = ['specificity_id', 'gwas'])
        #prioritization['Specificity_id'] = metadata[0]
        prioritization['GWAS'] = metadata[1]
        #display(priorization.head())
        prioritization_combined = pd.concat([prioritization_combined, prioritization])
        
        
    if '/conditional/' in file:
        conditional = pd.read_csv(file, sep = '\t', header = 0)
        file = file.split('/')
        metadata = file[-1]
        metadata = re.sub(r'\.cell_type_results.txt', '', metadata)
        metadata = metadata.split('__')
        conditional['GWAS'] = metadata[1]
        conditional['Conditional_annot'] = metadata[-1]
        conditional = conditional[['Name','GWAS', 'Conditional_annot', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
        conditional_combined = pd.concat([conditional_combined, conditional])

    if '/h2/' in file:
        h2 = pd.read_csv(file, sep = '\t', header = 0)
        file = file.split('/')
        metadata = file[-1].split('__')
        h2 = h2.iloc[-1:]
        h2['Specificity_id'], h2['GWAS'] = metadata[0], metadata[1]
        
        h2_combined = pd.concat([h2_combined, h2])


if not prioritization_combined.empty:
    tmp = prioritization_combined['Name'].str.split("__", 2, expand = True)
    prioritization_combined["Specificity_id"] = tmp[0]
    prioritization_combined["Annotation"] = tmp[1]
    prioritization_combined = prioritization_combined.drop(columns='Name')
    prioritization_combined = prioritization_combined[['GWAS', "Specificity_id", 'Annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
    prioritization_combined.to_csv(os.path.join(results_output_dir, 'prioritization.csv'), index = False, sep = '\t')

if not conditional_combined.empty:
    tmp = conditional_combined['Name'].str.split("__", 2, expand = True)
    conditional_combined["Specificity_id"] = tmp[0]
    conditional_combined["Annotation"] = tmp[1]
    conditional_combined = conditional_combined.drop(columns='Name')
    conditional_combined = conditional_combined[['GWAS', "Specificity_id", 'Conditional_annot', 'Annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]

    conditional_combined.to_csv(os.path.join(results_output_dir, 'conditional.csv'), index = False, sep = '\t')

if not h2_combined.empty:
    h2_combined = h2_combined[['GWAS', 'Specificity_id', 'Category', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error', 'Enrichment', 'Enrichment_std_error', 'Enrichment_p', 'Coefficient', 'Coefficient_std_error', 'Coefficient_z-score']]
    h2_combined.to_csv(os.path.join(results_output_dir, 'heritability.csv'), index = False, sep = '\t')

"""OLD - will be removed
if analysis_type == 'heritability':
"""
"""
    TODO: Throw this in a function
    
    args:
            analysis_type - if analysis_type is heriability, this function will be called
    
    output: CELLECT-ldsc/results/{prefix}__heriability.csv where all .results files will be concatenated in 
"""
"""
    h2_dir = snakemake.params['h2_dir']

    results_rows = []
    tmp_dict = {}

    print('Collecting all h2 results into one file now...')

    for h2_file in os.listdir(h2_dir):
        if h2_file.startswith(run_prefix) and h2_file.endswith('.results'):
            metadata = h2_file.split('__')
            results = [metadata[0], metadata[1]] #extract the dataset id and gwas
            with open(os.path.join(h2_dir, h2_file)) as f:
                results += list(f)[-1].split() #the last line of the .results has info on h2 of the given annotation. We take that line, strip it for delimiters sa \n and \t
                results_rows.append(results) 
                
                
    df = pd.DataFrame(results_rows, columns = ['Dataset_id',
                                                'Gwas', 
                                                'Category', 
                                                'Prop._SNPs', 
                                                'Prop._h2', 
                                                'Prop._h2_std_error', 
                                                'Enrichment', 
                                                'Enrichment_std_error',  
                                                'Enrichment_p',    
                                                'Coefficient',     
                                                'Coefficient_std_error',   
                                                'Coefficient_z-score'])
    
    path = os.path.join(results_output_dir, '{run_prefix}__heritability.csv'.format(run_prefix = run_prefix))
    df.to_csv(path, header = 0, index = False, sep = '\t')


else:

    for file in os.listdir(os.path.join(base_output_dir, analysis_type)):
        if file.endswith('cell_type_results.txt'):
            cell_type_result = pd.read_csv(file, sep = '\t', header = 0)
            
            metadata = re.sub('\.cell_type_results.txt$', '', file)
            metadata = re.split('__', metadata)

            df_metadata = pd.DataFrame([metadata * len(ccell_type_result), columns = 'Dataset_id', 'Gwas')

            cell_type_result = pd.merge(df_metadata, cell_type_result)

        combined_results_txt = pd.concat(cell_type_result)
"""



##
"""
To dos: Make a parser that takes in the list_target_files as input. 

"""
