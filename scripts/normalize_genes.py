import glob
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def count_16s_sequences():
    """
    Counts the number of 16s sequences in files matching '*.level_1.txt'.
    
    Returns:
        dict: A dictionary mapping sample IDs to the number of 16s sequences.
    """
    # Find all relevant MetaxaQR files
    files = glob.glob('*.level_1.txt')
    sequence_count_dict = {}
     
    # Load files and append sample_id and nr SSU to dict   
    for file in files:
        sample_id = file.split('.level_1.txt')[0]
        df = pd.read_csv(file, sep='\t', names=['kingdom', 'nr_ssu'])
        bacteria_count = list(df[df.kingdom == 'Bacteria'].nr_ssu)[0]
        # Dict with sample ID as key and nr bacterial SSU as value
        sequence_count_dict[sample_id] = bacteria_count
    
    return sequence_count_dict

def filter_group_for_best_alignment(group):
    """
    Filters a DataFrame group to find the best alignment based on criteria.

    Args:
        group (DataFrame): The DataFrame group to filter.

    Returns:
        DataFrame: A DataFrame containing the filtered results.
    """
    # Best his has lowest e-value, highest pident and length
    # If multiple hits has same value return a randomly selected best hit
    group = group[group['evalue'] == group['evalue'].min()]
    group = group[group['pident'] == group['pident'].max()]
    group = group[group['length'] == group['length'].max()]
    
    return group.sample(1, random_state=1)

def process_file(file, mode):
    """
    Processes a single file to group and filter data based on sequence identity.

    Args:
        file (str): The file path to process.
        mode (str): The mode of processing, affects file parsing.

    Returns:
        DataFrame: A DataFrame with processed and filtered data.
    """
    sample_id = file.split(f'_{mode}.tsv')[0]
    colnames = ['qseqid', 'sseqid', 'pident', 'evalue', 'length', 'slen']
    df = pd.read_csv(file, sep='\t', names=colnames)
    # Function for selecting best hit
    df_grouped = df.groupby('qseqid').apply(filter_group_for_best_alignment).reset_index(drop=True)
    df_grouped['sample_id'] = sample_id
    return df_grouped[['sseqid', 'slen', 'sample_id']]

def process_file_wrapper(file, mode):
    """
    Wrapper function for process_file to make it compatible with ProcessPoolExecutor.

    Args:
        file (str): The file path to process.
        mode (str): The mode of processing.

    Returns:
        DataFrame: A DataFrame with processed and filtered data.
    """
    return process_file(file, mode)

def get_files(mode):
    """
    Retrieves and processes files in parallel, concatenating the results.

    Args:
        mode (str): The mode of processing to determine file patterns.

    Returns:
        DataFrame: A concatenated DataFrame of all processed files.
    """
    files = glob.glob(f'*_{mode}.tsv')
    
    # Parallelizes processing files
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_file_wrapper, files, [mode]*len(files)))
    
    # Output here is the best hit for each read
    return pd.concat(results, ignore_index=True)

def normalize_gene_data(df_temp, nr_ssu):
    """
    Normalizes gene data based on sequence length and 16s sequence counts.

    Args:
        df_temp (DataFrame): The DataFrame with gene data.
        nr_ssu (dict): A dictionary mapping sample IDs to 16s sequence counts.

    Returns:
        DataFrame: A DataFrame with normalized gene data.
    """
    # Counts how many reads aligned with each ARG/MGE
    gene_counts = df_temp['sseqid'].value_counts().reset_index()
    gene_counts.columns = ["Genes", "Count"]
    
    # Add data like gene length etc
    df_temp = df_temp.merge(gene_counts, left_on='sseqid', right_on='Genes', how='right')
    
    # add nr of ssu for each sample
    df_temp['nr_ssu'] = df_temp['sample_id'].map(nr_ssu)
    
    # Normalize gene abundance  
    df_temp['Gene_normalized'] = (df_temp['Count'] / df_temp['slen']) / (df_temp['nr_ssu'] / 720)

    # Due to me writing the script in a strange way we will get duplicate rows with identical data
    # Output here is a non redundant df.
    df_temp.drop_duplicates(subset=['sseqid', 'sample_id'], inplace=True)
    
    return df_temp[['Genes', 'Count', 'nr_ssu', 'slen', 'Gene_normalized', 'sample_id']]

def read_and_normalize_files(mode, nr_ssu):
    """
    Reads and normalizes files based on the specified mode and 16s sequence counts.

    Args:
        mode (str): The mode of processing.
        nr_ssu (dict): A dictionary mapping sample IDs to 16s sequence counts.

    Returns:
        DataFrame: A DataFrame with normalized gene data.
    """
    df_temp = get_files(mode)
    return normalize_gene_data(df_temp, nr_ssu)

def assign_phenotypes_to_args(diamond_out):
    """
    Assigns phenotypes to ARGs based on gene accession numbers.

    Args:
        diamond_out (DataFrame): The DataFrame containing ARG data.
    """
    # Load phenotypic information from ResFinder
    cols = ['Gene_accession no.', 'Class']
    phenotype_data = pd.read_csv('phenotypes.txt', sep='\t', usecols=cols)
    
    # Handling multiple resistance phenotypes
    # looking for ',' since phenotypes are written "tetracycline, aminoglycosides" etc.
    # eg. If a gene confers resistance towards more than 1 antibiotic is it called "multiple"
    multiple_resistance_classes = [cls for cls in phenotype_data['Class'].unique() if ',' in cls]
    
    # Prodigal adds _ and a number to each translated gene name. 
    # Here I remove those characters from the gene name so it matches 
    # with the one in the ResFinder Phenotype file.
    diamond_out['Phenotype'] = diamond_out['Genes'].apply(lambda x: '_'.join(x.rsplit('_', 1)[0:-1]))
    
    # Replace phenotype to multiple if in multiple_resistance_classes list
    phenotype_data.loc[phenotype_data['Class'].isin(multiple_resistance_classes), 'Class'] = 'Multiple'
    
    # There are some genes which have multiple entries in the Phenotype list
    # Here I make sure that they are classified as Multiple
    duplicated_args = phenotype_data[phenotype_data.duplicated(['Gene_accession no.'])]
    phenotype_data.loc[phenotype_data['Gene_accession no.'].isin(duplicated_args["Gene_accession no."]), 'Class'] = 'Multiple'
    phenotype_data.drop_duplicates(subset=['Gene_accession no.'], inplace=True, keep='first')
    phenotype_data.set_index('Gene_accession no.', inplace=True)
    
    # Converting phenotype_data to dict and use it to translate gene manes to phenotype
    phenotype_dict = phenotype_data['Class'].to_dict()
    diamond_out.replace({"Phenotype": phenotype_dict}, inplace=True)
    diamond_out.to_csv('ARG_norm.tsv', sep='\t', index=False)

def main():
    """
    Main function to process files and infer phenotypes.
    """
    modes = ['ResFinder', 'mobileOG']
    for mode in modes: 
        SSU_counts = count_16s_sequences()
        normalized_files = read_and_normalize_files(mode, SSU_counts)  
        if mode == 'ResFinder':
            assign_phenotypes_to_args(normalized_files)
        elif mode == 'mobileOG':
            normalized_files.to_csv('MGE_norm.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()