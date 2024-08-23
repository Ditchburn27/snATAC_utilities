#### Functions for snATAC-seq preprocessing with snapATAC2
from pybedtools import BedTool
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

##############################################################
# Function for identifying genomic bins
# that overlap with features in provided
# bed file e.g. promoter, enhancers
def region_overlaps_bed(args):
    region, bed_file_bed = args
    chrom, start_end = region.split(':')
    start, end = map(int, start_end.split('-'))
    region_bed = BedTool(f"{chrom} {start} {end}", from_string=True)
    return bool(region_bed.intersect(bed_file_bed, u=True))

def select_bin_region_features(adata, region_dict):
    """
    Adds columns to adata.var indicating if each region overlaps with any regions in the provided BED files.
    
    Parameters:
    adata (AnnData): The AnnData object containing the data.
    region_dict (dict): A dictionary where keys are the output column names and values are the paths to BED files.
    region_column (str): The column in adata.var that contains the regions.
    """
    # Create regions column in adata.var
    adata.var['regions'] = adata.var_names
    # Loop through each item in the dictionary
    for output_column, bed_file_path in region_dict.items():
        # Read the BED file using pybedtools
        bed_file_bed = BedTool(bed_file_path)

        # Prepare arguments for multiprocessing
        args = [(region, bed_file_bed) for region in adata.var['regions']]

        # Use multiprocessing to process regions in parallel with progress bar
        with Pool(cpu_count()) as pool:
            results = list(tqdm(pool.imap(region_overlaps_bed, args), total=len(args), desc=f"Processing {output_column}"))

        # Add the results to the new column in adata.var
        adata.var[output_column] = results

##############################################################
# Function to print counts of regions that are identified to overlap with
# features from bed file. 
def count_region_features(adata, region_feature_label='selected'):
    true_count = adata.var[region_feature_label].sum()
    print(f"Number of genomic bins that are {region_feature_label}: {true_count}")

##############################################################
# Function for counting number of doublets and removing doublets
def count_remove_doublets(adata=adata, groupby='library_id'):
    """
    Count and remove doublets from an AnnData object.

    This function identifies and counts doublets in the `adata` object for each group specified
    by the `groupby` parameter, prints the number of doublets for each group, the total number
    of doublets, and removes these doublets from the dataset.

    Parameters:
    ----------
    adata : AnnData
        The AnnData object containing single-cell data. The `.obs` attribute should have a 
        boolean column 'predicted_doublet' that indicates whether each cell is a predicted doublet.
    groupby : str, optional (default: 'library_id')
        The column name in `adata.obs` by which to group the data and count doublets.

    Returns:
    -------
    None
        The function modifies the input `adata` in place by removing the doublets.
        
    Notes:
    -----
    - Ensure that `adata` contains a boolean column named 'predicted_doublet' in its `.obs` attribute.
    - The function prints the number of doublets for each group and the total number of doublets.
    - The function filters out the doublets from the dataset.
    """
    # Print number of doublets for each library
    doublet_counts = adata.obs.groupby(groupby)['predicted_doublet'].sum()
    print('Doublets identified for each sample:')
    for library_id, counts in doublet_counts.items():
        print(f'{groupby}: {library_id}, number of doublets = {counts}')
    total = doublet_counts.sum()
    print(f'Total number of doublets = {total}')
    # Remove doublets from dataset
    dub_mk = adata.obs['predicted_doublet']
    adata = adata[~dub_mk]
    return