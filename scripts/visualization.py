"""
Visualization functions for variant analysis.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def set_style():
    """Set matplotlib and seaborn style for publication-quality figures."""
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    plt.rcParams['figure.figsize'] = (12, 6)
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.labelsize'] = 11
    plt.rcParams['axes.titlesize'] = 12


def plot_gene_distribution(gene_summary: pd.DataFrame, 
                          top_n: int = 15,
                          save_path: str = None):
    """
    Plot variant distribution across genes.
    
    Parameters
    ----------
    gene_summary : pd.DataFrame
        Gene summary statistics.
    top_n : int
        Number of top genes to display.
    save_path : str, optional
        Path to save figure.
    """
    set_style()
    
    gene_summary_top = gene_summary.head(top_n).copy()
    gene_summary_top = gene_summary_top.sort_values('TotalVariants')
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    x = np.arange(len(gene_summary_top))
    width = 0.35
    
    ax.barh(x - width/2, gene_summary_top['PathogenicCount'], width, 
            label='Pathogenic', color='#d62728', alpha=0.8)
    ax.barh(x + width/2, gene_summary_top['BenignCount'], width, 
            label='Benign', color='#2ca02c', alpha=0.8)
    
    ax.set_yticks(x)
    ax.set_yticklabels(gene_summary_top['Gene'])
    ax.set_xlabel('Number of Variants')
    ax.set_title(f'Variant Distribution Across Top {top_n} Genes\n(Pathogenic vs Benign)')
    ax.legend()
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax


def plot_variant_type_distribution(variant_type_summary: pd.DataFrame,
                                   save_path: str = None):
    """
    Plot distribution of variant types by clinical significance.
    
    Parameters
    ----------
    variant_type_summary : pd.DataFrame
        Variant type summary from analysis.
    save_path : str, optional
        Path to save figure.
    """
    set_style()
    
    df = variant_type_summary.copy()
    
    # Ensure we have the right columns
    if 'Pathogenic' in df.columns and 'Benign' in df.columns:
        df_plot = df[['Pathogenic', 'Benign']].sort_values('Pathogenic', ascending=False)
    else:
        df_plot = df.iloc[:, :2]
    
    fig, ax = plt.subplots(figsize=(11, 6))
    
    df_plot.plot(kind='bar', ax=ax, color=['#d62728', '#2ca02c'], alpha=0.8)
    
    ax.set_title('Variant Type Distribution (Pathogenic vs Benign)')
    ax.set_xlabel('Variant Type')
    ax.set_ylabel('Number of Variants')
    ax.legend(title='Clinical Significance')
    ax.grid(axis='y', alpha=0.3)
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax


def plot_clinical_significance_pie(clinical_sig_summary: pd.DataFrame,
                                  save_path: str = None):
    """
    Plot pie chart of clinical significance distribution.
    
    Parameters
    ----------
    clinical_sig_summary : pd.DataFrame
        Clinical significance summary.
    save_path : str, optional
        Path to save figure.
    """
    set_style()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    colors = ['#d62728', '#2ca02c', '#1f77b4']
    
    ax.pie(clinical_sig_summary['Count'], 
           labels=clinical_sig_summary['ClinicalSignificance'],
           autopct='%1.1f%%',
           colors=colors[:len(clinical_sig_summary)],
           startangle=90)
    
    ax.set_title('Distribution of Clinical Significance Categories')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax


def plot_pathogenic_to_benign_ratio(gene_summary: pd.DataFrame,
                                   save_path: str = None):
    """
    Plot pathogenic-to-benign ratio by gene.
    
    Parameters
    ----------
    gene_summary : pd.DataFrame
        Gene summary statistics.
    save_path : str, optional
        Path to save figure.
    """
    set_style()
    
    gene_summary_sorted = gene_summary.sort_values('PathogenicToBenignRatio', ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = ['#d62728' if x > 1 else '#2ca02c' for x in gene_summary_sorted['PathogenicToBenignRatio']]
    
    ax.barh(gene_summary_sorted['Gene'], gene_summary_sorted['PathogenicToBenignRatio'], 
            color=colors, alpha=0.8)
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1, label='Ratio = 1')
    
    ax.set_xlabel('Pathogenic-to-Benign Ratio')
    ax.set_title('Pathogenic-to-Benign Variant Ratio by Gene')
    ax.legend()
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax


def plot_allele_frequency_spectrum(df: pd.DataFrame,
                                   save_path: str = None):
    """
    Plot allele frequency spectrum comparing pathogenic and benign variants.
    
    Parameters
    ----------
    df : pd.DataFrame
        Processed variant data frame with AlleleFreq column.
    save_path : str, optional
        Path to save figure.
    """
    set_style()
    
    if 'AlleleFreq' not in df.columns:
        print("AlleleFreq column not available. Skipping allele frequency plot.")
        return None, None
    
    df_freq = df.copy()
    df_freq['AlleleFreq'] = pd.to_numeric(df_freq['AlleleFreq'], errors='coerce')
    
    pathogenic = df_freq[df_freq['ClinicalSignif'] == 'Pathogenic']['AlleleFreq'].dropna()
    benign = df_freq[df_freq['ClinicalSignif'] == 'Benign']['AlleleFreq'].dropna()
    
    if len(pathogenic) == 0 or len(benign) == 0:
        print("Insufficient allele frequency data for comparison. Skipping plot.")
        return None, None
    
    fig, ax = plt.subplots(figsize=(11, 6))
    
    ax.hist(pathogenic, bins=30, alpha=0.6, label=f'Pathogenic (n={len(pathogenic)})', 
            color='#d62728', edgecolor='black')
    ax.hist(benign, bins=30, alpha=0.6, label=f'Benign (n={len(benign)})', 
            color='#2ca02c', edgecolor='black')
    
    ax.set_xlabel('Allele Frequency')
    ax.set_ylabel('Number of Variants')
    ax.set_title('Allele Frequency Spectrum: Pathogenic vs Benign Variants')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax


def plot_consequence_summary(consequence_summary: pd.DataFrame,
                            top_n: int = 15,
                            save_path: str = None):
    """
    Plot functional consequence distribution.
    
    Parameters
    ----------
    consequence_summary : pd.DataFrame
        Consequence summary from analysis.
    top_n : int
        Number of top consequences to display.
    save_path : str, optional
        Path to save figure.
    """
    set_style()
    
    if consequence_summary.empty:
        print("Consequence summary is empty. Skipping plot.")
        return None, None
    
    df = consequence_summary.copy()
    
    # Get top consequences
    if 'Total' in df.columns:
        df_top = df.nlargest(top_n, 'Total')
    else:
        df_top = df.head(top_n)
    
    # Ensure we have Pathogenic and Benign columns
    if 'Pathogenic' in df_top.columns and 'Benign' in df_top.columns:
        df_top = df_top.sort_values('Pathogenic')
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        x = np.arange(len(df_top))
        width = 0.35
        
        ax.barh(x - width/2, df_top['Pathogenic'], width, 
                label='Pathogenic', color='#d62728', alpha=0.8)
        ax.barh(x + width/2, df_top['Benign'], width, 
                label='Benign', color='#2ca02c', alpha=0.8)
        
        ax.set_yticks(x)
        ax.set_yticklabels(df_top.index)
        ax.set_xlabel('Number of Variants')
        ax.set_title(f'Functional Consequence Distribution (Top {top_n})')
        ax.legend()
        ax.grid(axis='x', alpha=0.3)
    else:
        fig, ax = plt.subplots(figsize=(10, 8))
        df_top.plot(kind='barh', ax=ax, color='#1f77b4', alpha=0.8)
        ax.set_title(f'Functional Consequence Distribution (Top {top_n})')
        ax.set_xlabel('Number of Variants')
        ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax
