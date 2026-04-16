"""
Analysis functions for variant statistics and summaries.
"""

import pandas as pd
import numpy as np


def gene_summary_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate gene-level summary statistics.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    pd.DataFrame
        Gene summary with counts and statistics.
    """
    summary = df.groupby('Gene').agg(
        TotalVariants=('Gene', 'count'),
        PathogenicCount=('ClinicalSignif', lambda x: (x == 'Pathogenic').sum()),
        BenignCount=('ClinicalSignif', lambda x: (x == 'Benign').sum()),
        ChromosomesAffected=('Chromosome', 'nunique'),
        VariantTypesCount=('VariantType', 'nunique'),
    ).reset_index()
    
    # Calculate pathogenic-to-benign ratio
    summary['PathogenicToBenignRatio'] = summary['PathogenicCount'] / (summary['BenignCount'] + 1)
    
    # Calculate percentages
    summary['PercentPathogenic'] = (summary['PathogenicCount'] / summary['TotalVariants'] * 100).round(1)
    summary['PercentBenign'] = (summary['BenignCount'] / summary['TotalVariants'] * 100).round(1)
    
    # Sort by total variants descending
    summary = summary.sort_values('TotalVariants', ascending=False).reset_index(drop=True)
    
    return summary


def chromosome_summary_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate chromosome-level summary statistics.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    pd.DataFrame
        Chromosome summary.
    """
    summary = df.groupby('Chromosome').agg(
        TotalVariants=('Chromosome', 'count'),
        GenesAffected=('Gene', 'nunique'),
        PathogenicCount=('ClinicalSignif', lambda x: (x == 'Pathogenic').sum()),
        BenignCount=('ClinicalSignif', lambda x: (x == 'Benign').sum()),
    ).reset_index()
    
    summary['PercentPathogenic'] = (summary['PathogenicCount'] / summary['TotalVariants'] * 100).round(1)
    summary = summary.sort_values('TotalVariants', ascending=False).reset_index(drop=True)
    
    return summary


def variant_type_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate variant type distribution summary.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    pd.DataFrame
        Variant type summary split by clinical significance.
    """
    summary = df.groupby(['VariantType', 'ClinicalSignif']).size().reset_index(name='Count')
    summary_pivot = summary.pivot(index='VariantType', columns='ClinicalSignif', values='Count').fillna(0)
    
    # Calculate totals and percentages
    if 'Pathogenic' in summary_pivot.columns and 'Benign' in summary_pivot.columns:
        summary_pivot['Total'] = summary_pivot['Pathogenic'] + summary_pivot['Benign']
        summary_pivot['PercentPathogenic'] = (
            summary_pivot['Pathogenic'] / summary_pivot['Total'] * 100
        ).round(1)
    
    summary_pivot = summary_pivot.sort_values('Total', ascending=False)
    
    return summary_pivot


def functional_consequence_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate functional consequence distribution summary.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    pd.DataFrame
        Functional consequence summary.
    """
    # Check if Consequence column exists and is not empty
    if 'Consequence' not in df.columns or df['Consequence'].isna().all():
        print("Warning: Consequence column not available or empty. Skipping consequence analysis.")
        return pd.DataFrame()
    
    # Parse consequences (may be comma-separated)
    df_consequence = df.copy()
    df_consequence['Consequence'] = df_consequence['Consequence'].fillna('Unknown')
    
    # Expand consequences (each consequence gets its own row)
    df_consequence['Consequence'] = df_consequence['Consequence'].str.split(',')
    df_consequence = df_consequence.explode('Consequence')
    df_consequence['Consequence'] = df_consequence['Consequence'].str.strip()
    
    summary = df_consequence.groupby(['Consequence', 'ClinicalSignif']).size().reset_index(name='Count')
    summary_pivot = summary.pivot(index='Consequence', columns='ClinicalSignif', values='Count').fillna(0)
    
    if 'Pathogenic' in summary_pivot.columns and 'Benign' in summary_pivot.columns:
        summary_pivot['Total'] = summary_pivot['Pathogenic'] + summary_pivot['Benign']
        summary_pivot['PercentPathogenic'] = (
            summary_pivot['Pathogenic'] / summary_pivot['Total'] * 100
        ).round(1)
        summary_pivot = summary_pivot.sort_values('Total', ascending=False)
    else:
        summary_pivot = summary_pivot.sort_values(summary_pivot.columns[0], ascending=False)
    
    return summary_pivot


def clinical_significance_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate clinical significance distribution summary.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    pd.DataFrame
        Clinical significance summary.
    """
    summary = df['ClinicalSignif'].value_counts().reset_index()
    summary.columns = ['ClinicalSignificance', 'Count']
    summary['Percentage'] = (summary['Count'] / len(df) * 100).round(1)
    summary = summary.sort_values('Count', ascending=False).reset_index(drop=True)
    
    return summary


def allele_frequency_summary(df: pd.DataFrame) -> dict:
    """
    Generate allele frequency statistics if available.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    dict
        Summary statistics for allele frequency by clinical significance.
    """
    if 'AlleleFreq' not in df.columns:
        print("AlleleFreq column not available.")
        return None
    
    # Convert to numeric, coercing errors
    df_freq = df.copy()
    df_freq['AlleleFreq'] = pd.to_numeric(df_freq['AlleleFreq'], errors='coerce')
    
    # Separate pathogenic and benign
    pathogenic = df_freq[df_freq['ClinicalSignif'] == 'Pathogenic']['AlleleFreq'].dropna()
    benign = df_freq[df_freq['ClinicalSignif'] == 'Benign']['AlleleFreq'].dropna()
    
    summary = {
        'pathogenic': {
            'count': len(pathogenic),
            'mean': pathogenic.mean(),
            'median': pathogenic.median(),
            'min': pathogenic.min(),
            'max': pathogenic.max(),
        },
        'benign': {
            'count': len(benign),
            'mean': benign.mean(),
            'median': benign.median(),
            'min': benign.min(),
            'max': benign.max(),
        },
        'coverage_percent': {
            'pathogenic': (len(pathogenic) / len(df_freq[df_freq['ClinicalSignif'] == 'Pathogenic']) * 100).round(1),
            'benign': (len(benign) / len(df_freq[df_freq['ClinicalSignif'] == 'Benign']) * 100).round(1),
        }
    }
    
    return summary


def variant_count_summary(df: pd.DataFrame) -> dict:
    """
    Generate overall variant count summary.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame.
    
    Returns
    -------
    dict
        Summary counts.
    """
    summary = {
        'total_variants': len(df),
        'total_genes': df['Gene'].nunique(),
        'total_chromosomes': df['Chromosome'].nunique(),
        'pathogenic': len(df[df['ClinicalSignif'] == 'Pathogenic']),
        'benign': len(df[df['ClinicalSignif'] == 'Benign']),
        'variant_types': df['VariantType'].nunique(),
    }
    
    summary['percent_pathogenic'] = round(summary['pathogenic'] / summary['total_variants'] * 100, 1)
    summary['percent_benign'] = round(summary['benign'] / summary['total_variants'] * 100, 1)
    return summary
    
   