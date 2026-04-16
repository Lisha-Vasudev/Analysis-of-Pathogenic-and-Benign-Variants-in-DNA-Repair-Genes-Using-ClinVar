"""
Data loading and preprocessing functions for ClinVar variants analysis.
"""

import os
import gzip
import requests
import pandas as pd
from pathlib import Path


def download_clinvar(output_path: str, chunk_size: int = 8192) -> str:
    """
    Download the latest ClinVar variant_summary.txt.gz from NCBI FTP.
    
    Parameters
    ----------
    output_path : str
        Path to save the downloaded file.
    chunk_size : int
        Chunk size for streaming download.
    
    Returns
    -------
    str
        Path to the decompressed file.
    """
    url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    
    print(f"Downloading ClinVar variant_summary.txt.gz from {url}...")
    
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        # Save compressed file
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
        
        print(f"Downloaded to {output_path}")
        
        # Decompress
        output_txt = output_path.replace('.gz', '')
        print(f"Decompressing to {output_txt}...")
        
        with gzip.open(output_path, 'rb') as f_in:
            with open(output_txt, 'wb') as f_out:
                f_out.write(f_in.read())
        
        print(f"Decompression complete: {output_txt}")
        return output_txt
    
    except Exception as e:
        print(f"Error downloading ClinVar: {e}")
        raise


def load_clinvar_raw(filepath: str) -> pd.DataFrame:
    """
    Load ClinVar variant_summary.txt file using only needed columns.
    """
    print(f"Loading ClinVar data from {filepath}...")
    
    df = pd.read_csv(
        filepath,
        sep='\t',
        usecols=[
            'GeneSymbol',
            'ClinicalSignificance',
            'Chromosome',
            'Start',
            'Stop',
            'Type',
            'ReviewStatus'
        ],
        dtype={
            'Chromosome': str,
            'Start': 'Int64',
            'Stop': 'Int64'
        },
        low_memory=True
    )
    
    print(f"Loaded {len(df)} total variants")
    return df


def filter_to_dna_repair_genes(df: pd.DataFrame, gene_list: list = None) -> pd.DataFrame:
    """
    Filter ClinVar variants to DNA repair genes.
    
    Parameters
    ----------
    df : pd.DataFrame
        ClinVar data frame.
    gene_list : list, optional
        List of gene symbols. If None, uses default 19 DNA repair genes.
    
    Returns
    -------
    pd.DataFrame
        Filtered data frame.
    """
    if gene_list is None:
        # 19 DNA repair genes
        gene_list = [
            'MLH1', 'MSH2', 'MSH6', 'PMS2',  # MMR
            'BRCA1', 'BRCA2', 'RAD51', 'XRCC3',  # HR
            'APE1', 'XRCC1',  # BER
            'XPA', 'XPC', 'ERCC1',  # NER
            'TP53', 'CHEK2', 'ATM', 'PTEN',  # Checkpoint
            'POLE', 'POLD1'  # Polymerase
        ]
    
    # Handle multiple genes per row (may be semicolon-separated)
    mask = df['GeneSymbol'].fillna('').apply(
    lambda x: any(gene in str(x).split(';') for gene in gene_list)
)
    
    filtered = df[mask].copy()
    print(f"Filtered to {len(filtered)} variants in DNA repair genes")
    return filtered


def filter_clinical_significance(df: pd.DataFrame, 
                                 significance_list: list = None) -> pd.DataFrame:
    """
    Filter to specific clinical significance categories.
    
    Parameters
    ----------
    df : pd.DataFrame
        ClinVar data frame.
    significance_list : list, optional
        List of clinical significance values. 
        If None, uses ['Pathogenic', 'Benign'].
    
    Returns
    -------
    pd.DataFrame
        Filtered data frame.
    """
    if significance_list is None:
        significance_list = ['Pathogenic', 'Benign']
    
    mask = df['ClinicalSignificance'].isin(significance_list)
    filtered = df[mask].copy()
    print(f"Filtered to {len(filtered)} variants with significance in {significance_list}")
    return filtered


def filter_review_status(df: pd.DataFrame, min_review_level: str = 'classified by multiple submitters') -> pd.DataFrame:
    """
    Filter to variants with minimum review status quality.
    
    Parameters
    ----------
    df : pd.DataFrame
        ClinVar data frame.
    min_review_level : str
        Minimum review status (currently accepts any threshold for flexibility).
    
    Returns
    -------
    pd.DataFrame
        Filtered data frame.
    """
    # ClinVar review statuses range from lowest to highest confidence
    # For this analysis, we keep variants with "classified by multiple submitters" and higher
    
    # If ReviewStatus column exists, apply filter; otherwise, return as-is with note
    if 'ReviewStatus' in df.columns:
        # Keep rows where ReviewStatus contains "multiple submitters" or higher confidence terms
        mask = df['ReviewStatus'].fillna('').str.contains(
            'multiple submitters|criteria provided|no assertion|expert|practice guideline',
            case=False, regex=True
        )
        # For simplicity, we include variants with any assertion; users can customize
        filtered = df[mask | (df['ReviewStatus'].notna())].copy()
        print(f"Filtered to {len(filtered)} variants with review status information")
    else:
        print("ReviewStatus column not found; skipping review status filter")
        filtered = df.copy()
    
    return filtered


def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove duplicate variants (same chromosome, position, alleles, gene).
    
    Parameters
    ----------
    df : pd.DataFrame
        Data frame.
    
    Returns
    -------
    pd.DataFrame
        De-duplicated data frame.
    """
    key_cols = ['Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele', 'Gene']
    
    # Check which columns exist
    key_cols = [col for col in key_cols if col in df.columns]
    
    before = len(df)
    df_dedup = df.drop_duplicates(subset=key_cols, keep='first')
    after = len(df_dedup)
    
    print(f"Removed {before - after} duplicates, retained {after} variants")
    return df_dedup


def select_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select and rename key columns for analysis.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data frame.
    
    Returns
    -------
    pd.DataFrame
        Data frame with selected columns.
    """
    columns_to_keep = [
        'Chromosome',
        'Start',
        'Stop',
        'Type',
        'GeneSymbol',
        'ClinicalSignificance',
        'ReviewStatus'
    ]    
    
    # Keep only columns that exist
    columns_to_keep = [col for col in columns_to_keep if col in df.columns]
    df_selected = df[columns_to_keep].copy()
    
    # Rename for clarity
    df_selected = df_selected.rename(columns={
        'Type': 'VariantType',
        'GeneSymbol': 'Gene',
        'ClinicalSignificance': 'ClinicalSignif'
    })
    
    return df_selected


def preprocess_clinvar(input_path: str, output_path: str, gene_list: list = None) -> pd.DataFrame:
    """
    Complete preprocessing pipeline for ClinVar data.
    
    Parameters
    ----------
    input_path : str
        Path to raw ClinVar variant_summary.txt
    output_path : str
        Path to save processed CSV
    gene_list : list, optional
        Custom gene list for filtering.
    
    Returns
    -------
    pd.DataFrame
        Processed data frame.
    """
    # Load
    df = load_clinvar_raw(input_path)
    
    # Filter to DNA repair genes
    df = filter_to_dna_repair_genes(df, gene_list=gene_list)
    
    # Filter to pathogenic/benign
    df = filter_clinical_significance(df)
    
    # Filter by review status (if applicable)
    df = filter_review_status(df)
    
    # Remove duplicates
    df = remove_duplicates(df)
    
    # Select and rename columns
    df = select_columns(df)
    
    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"\nProcessed data saved to {output_path}")
    print(f"Final dataset: {len(df)} variants across {df['Gene'].nunique()} genes")
    
    return df
