# Raw Data

## ClinVar Variant Summary

This directory should contain the ClinVar variant summary data used for analysis.

### How to obtain the data:

1. **Download ClinVar variant_summary.txt**
   - Source: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
   - File: `variant_summary.txt.gz` (latest release)
   - Size: ~200 MB (compressed), ~700 MB (uncompressed)

2. **Option A: Automatic download (recommended)**
   - The data preparation notebook (`notebooks/01_data_preparation.ipynb`) includes a function to automatically download and decompress this file.
   - Simply run the notebook; it will place the file here automatically.

3. **Option B: Manual download**
   - Download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/variant_summary.txt.gz
   - Decompress: `gunzip variant_summary.txt.gz`
   - Place `variant_summary.txt` in this directory

### Data structure:

The variant_summary.txt file is a tab-separated file with the following key columns:
- `#Chromosome`: Chromosome number or identifier
- `Start`: Variant start position (1-based)
- `Stop`: Variant stop position (1-based)
- `ReferenceAllele`: Reference allele sequence
- `AlternateAllele`: Alternate allele sequence
- `VariationType`: Type of variant (Insertion, Deletion, Substitution, etc.)
- `Gene` or `GeneSymbol`: Gene symbol(s) associated with the variant
- `ClinicalSignificance`: Reported clinical significance (Pathogenic, Benign, Uncertain significance, etc.)
- `ReviewStatus`: Assertion criteria review status
- `AlleleFrequency`: Allele frequency (if available; often sparse)
- Additional columns for reference information

### Processing in this project:

Notebook 01_data_preparation.ipynb will:
1. Download or load this file
2. Filter to DNA repair gene list (19 genes: MLH1, MSH2, MSH6, PMS2, BRCA1, BRCA2, RAD51, XRCC3, APE1, XRCC1, XPA, XPC, ERCC1, TP53, CHEK2, ATM, PTEN, POLE, POLD1)
3. Filter to ClinicalSignificance = "Pathogenic" or "Benign"
4. Filter to ReviewStatus = "classified by multiple submitters" or higher
5. Remove duplicates
6. Output processed CSV to `data/processed/`

### License and attribution:

ClinVar data is in the public domain. Attribution: "ClinVar at NCBI (https://www.ncbi.nlm.nih.gov/clinvar/)"
