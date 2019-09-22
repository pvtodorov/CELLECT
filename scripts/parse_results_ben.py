import pandas as pd
import os
import re
from glob import glob

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

base_output_dir = snakemake.params['BASE_OUTPUT_DIR']
results_output_dir = snakemake.params['results_out_dir'] 
#run_prefix = snakemake.params['run_prefix']
analysis_types = snakemake.params['analysis_types_performed']
#gwas = snakemake.params['gwas']
# file



for analysis_type in analysis_types:
    print('Compiling {} result files...'.format(analysis_type))

    if analysis_type == 'prioritization':
        
        prioritization_combined = pd.DataFrame()
        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', analysis_type, '*cell_type_results.txt'))
        result_files = [f for f in glob(result_path)]
    
        for f in result_files:
            prioritization = pd.read_csv(f, sep = '\t', header = 0)
            
            
            f = f.split('/')
            metadata = f[-1].replace('.cell_type_results.txt', '') # modules.mousebrain_bmicelltypes__T1D_Bradfield2011        
            metadata = metadata.split('__') # ['modules.mousebrain_bmicelltypes', 'T1D_Bradfield2011']
            prioritization['gwas'] = metadata[-1] # ['T1D_Bradfield2011']
            
            prioritization_combined = pd.concat([prioritization_combined, prioritization])
        
        specificity_id_annotation = prioritization_combined['Name'].str.split('__', expand = True) #https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
        prioritization_combined['specificity_id'] = specificity_id_annotation[0] #mousebrain
        prioritization_combined['annotation'] = specificity_id_annotation[1] #TEGLU32
        prioritization_combined = prioritization_combined[['gwas', "specificity_id", 'annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
        prioritization_combined.rename(columns = {'Coefficient': 'tau', 'Coefficient_std_error': 'se', 'Coefficient_P_value': 'pvalue'}, inplace = True)

        prioritization_combined.to_csv(os.path.join(results_output_dir, 'prioritization.csv'), index = False, sep = '\t')


    if analysis_type == 'conditional':
        
        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', analysis_type, '*cell_type_results.txt'))
        result_files = [f for f in glob(result_path)]
        
        conditional_combined = pd.DataFrame()

        for f in result_files:
            
            conditional = pd.read_csv(f, sep = '\t', header = 0)
            
            f = f.split('/')
            metadata = f[-1].replace('.cell_type_results.txt', '') # tabula_muris__EA3_Lee2018__CONDITIONAL__Brain_Non-Myeloid.oligodendrocyte      
            metadata = metadata.split('__') # ['tabula_muris', 'EA3_Lee2018', 'CONDITIONAL', 'Brain_Non-Myeloid.oligodendrocyte']
    
            conditional['gwas'] = metadata[1] # 'EA3_Lee2018'
            conditional['conditional_annotation'] = metadata[-1] # 'Brain_Non-Myeloid.oligodendrocyte'
            
            conditional_combined = pd.concat([conditional_combined, conditional])
        
        
        specificity_id_annotation = conditional_combined['Name'].str.split('__', expand = True) #https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
        conditional_combined['specificity_id'] = specificity_id_annotation[0] #mousebrain
        conditional_combined['annotation'] = specificity_id_annotation[1] #TEGLU32
        conditional_combined = conditional_combined[['gwas', "specificity_id", 'conditional_annotation', 'annotation', 'Coefficient', 'Coefficient_std_error',  'Coefficient_P_value']]
        conditional_combined.rename(columns = {'Coefficient': 'tau', 'Coefficient_std_error': 'se', 'Coefficient_P_value': 'pvalue'}, inplace = True)   

        conditional_combined.to_csv(os.path.join(results_output_dir, 'conditional.csv'), index = False, sep = '\t')


    if analysis_type == 'heritability': 

        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', 'h2', '*.results')) #'*results'))
        result_files = [f for f in glob(result_path)]
        
        h2_combined = pd.DataFrame() 

        for f in result_files:
            
            h2 = pd.read_csv(f, sep = '\t', header = 0)
            h2 = h2.tail(1) #.iloc[-1:] didn't exactly as intended, so we try this alternative method

            f = f.split('/')
            metadata = f[-1].replace('.results', '') # tabula_muris__EA3_Lee2018__CONDITIONAL__Brain_Non-Myeloid.oligodendrocyte      
            metadata = metadata.split('__') # ['tabula_muris', 'EA3_Lee2018', 'CONDITIONAL', 'Brain_Non-Myeloid.oligodendrocyte']
            h2['specificity_id'] = metadata[0]
            h2['gwas'] = metadata[1] # 'EA3_Lee2018'
            h2['annotation'] = metadata[-1] # 'Brain_Non-Myeloid.oligodendrocyte'
           
            h2_combined = pd.concat([h2_combined, h2])
           
        h2_combined = h2_combined[['gwas', "specificity_id", 'annotation', 'Prop._SNPs', 'Prop._h2', 'Prop._h2_std_error', 'Enrichment' , 'Enrichment_std_error', 'Enrichment_p']]
        h2_combined.rename(columns = {'Enrichment': 'h2_enrichment', 'Enrichment_std_error': 'h2_enrichment_se', 'Enrichment_p': 'h2_enrichment_pvalue'}, inplace = True)

        h2_combined.to_csv(os.path.join(results_output_dir, 'heritability.csv'), index = False, sep = '\t')


    if analysis_type == 'heritability_intervals':
    
        h2_int_combined = pd.DataFrame()

        result_path = os.path.abspath(os.path.join(base_output_dir, 'out', 'h2', '*results_intervals'))
        result_files = [f for f in glob(result_path )]
        
        for f in result_files:
            
            h2_int = pd.read_csv(f, sep = '\t', header = 0)
            h2_int['q'] = pd.Series(range(h2_int.shape[0]))
            
            f = f.split('/')
            metadata = f[-1].replace('.results_intervals', '') # tabula_muris__EA3_Lee2018__CONDITIONAL__Brain_Non-Myeloid.oligodendrocyte      
            metadata = metadata.split('__') #['tabula_muris', 'BMI_UKBB_Loh2018', 'h2_intervals', 'Liver.hepatocyte.qfixed']
            
            h2_int['specificity_id'] = metadata[0]
            h2_int['gwas'] = metadata[1] # 'BMI_UKBB_Loh2018'
            annotation = metadata[-1] # 'Liver.hepatocyte.qfixed'
            annotation = annotation.split('.')  # ['Liver', 'hepatocyte' , 'qfixed']
            annotation = '.'.join(annotation[:-1]) # 'Liver.hepatocyte'
            h2_int['annotation'] = annotation
            
            h2_int_combined = pd.concat([h2_int_combined, h2_int])
        
            
        h2_int_combined = h2_int_combined[['specificity_id', 'gwas', 'annotation', 'q', 'h2g', 'h2g_se', 'prop_h2g', 'prop_h2g_se', 'enr', 'enr_se', 'enr_pval']]

        h2_int_combined.to_csv(os.path.join(results_output_dir, 'heritability_intervals.csv'), index = False, sep = '\t')

    print('{} result files done compiling.'.format(analysis_type))



#analysis_type_dir = os.path.join(base_output_dir, analysis_type)
"""
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
"""

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
