# Methods and Interpretation

## Study Design

This project employs a descriptive bioinformatics analysis of curated genetic variant data to systematically characterize pathogenic and benign variants in clinically important DNA repair genes. The study is observational and comparative in nature, using public data without original experimental work.

---

## Data Source and Filtering

### ClinVar Data

**Source:** ClinVar (NCBI Variation Services)
- **URL:** https://www.ncbi.nlm.nih.gov/clinvar/
- **Data file:** variant_summary.txt (downloaded from FTP: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/)
- **Version:** Latest release at time of analysis
- **Size:** ~700 MB uncompressed; ~200 MB compressed

**Citation:**
Landrum, M. J., et al. (2023). ClinVar at NCBI: Improving access and support for clinical variant interpretation. _Human Mutation_, 44(5), 1108–1120. https://doi.org/10.1002/humu.24614

### Gene Selection

A curated set of 19 clinically and biologically important DNA repair genes were selected across four major DNA repair pathways:

| Pathway | Genes | Relevance |
|---------|-------|-----------|
| **Mismatch Repair (MMR)** | MLH1, MSH2, MSH6, PMS2 | Lynch syndrome (hereditary colorectal cancer); HNPCC |
| **Homologous Recombination (HR)** | BRCA1, BRCA2, RAD51, XRCC3 | Breast/ovarian cancer susceptibility; genome stability |
| **Base Excision Repair (BER)** | APE1, XRCC1 | Oxidative damage repair; cancer prevention |
| **Nucleotide Excision Repair (NER)** | XPA, XPC, ERCC1 | UV damage, xenobiotic repair; xeroderma pigmentosum |
| **Checkpoint & Tumor Suppressors** | TP53, CHEK2, ATM, PTEN | Cell cycle control; cancer susceptibility |
| **DNA Polymerase Fidelity** | POLE, POLD1 | Replication fidelity; polymerase-associated polyposis |

### Filtering Criteria

Variants were filtered using the following criteria (applied sequentially):

1. **Gene filter:** Restrict to 19 selected DNA repair genes
2. **Clinical significance filter:** Keep only "Pathogenic" or "Benign" classifications
   - Exclude "Uncertain significance" to enable clear binary comparison
   - Rationale: Uncertain variants introduce ambiguity; binary classification supports interpretable comparison
3. **Review status filter:** Retain variants with quality assertions
   - Include "classified by multiple submitters" or higher confidence levels
   - Rationale: Multi-submitter consensus indicates higher confidence
4. **Deduplication:** Remove exact duplicate entries (same chromosome, position, alleles, gene)

**Expected outcome after filtering:** 5,000–15,000 variants representing a balance of pathogenic and benign variants in these genes.

---

## Data Processing and Validation

### Column Standardization

ClinVar raw data columns were standardized for clarity:
- `#Chromosome` → `Chromosome`
- `ReferenceAllele` → `Ref`
- `AlternateAllele` → `Alt`
- `VariationType` → `VariantType`
- `ClinicalSignificance` → `ClinicalSignif`
- `AlleleFrequency` → `AlleleFreq`

### Data Quality Checks

The following checks were performed:
1. **Completeness:** Assess missing values in key columns (gene, chromosome, position, clinical significance)
2. **Coordinate validation:** Ensure Start ≤ Stop and both are positive integers
3. **Duplicate detection:** Identify and remove exact duplicate records
4. **Data type validation:** Confirm chromosome and position columns are numeric/categorical

### Allele Frequency Handling

Allele frequency (AF) data from ClinVar was evaluated for completeness:
- If ≥10% of variants have AF data: Conditional analysis performed with transparency
- If <10% AF coverage: AF analysis omitted and documented as limitation
- **No external frequency data merged:** gnomAD, 1000 Genomes, or other population databases were not integrated
- Rationale: Use only data originally present in ClinVar to avoid data integration complexity

---

## Statistical Analysis

### Descriptive Statistics

**Variant counts:**
- Total variants retained post-filtering
- Count by clinical significance (pathogenic, benign)
- Count by variant type (insertion, deletion, substitution, complex)
- Count per gene (rank order)
- Count per chromosome

**Proportions:**
- Percentage pathogenic: (# pathogenic) / (# total) × 100
- Percentage benign: (# benign) / (# total) × 100

### Gene-Level Analysis

**Pathogenic-to-benign ratio per gene:**
$$\text{P/B Ratio}_{\text{gene}} = \frac{\text{# pathogenic variants in gene}}{\text{# benign variants in gene} + 1}$$

The +1 denominator (Laplace smoothing) prevents division by zero for genes with zero benign variants.

**Interpretation:**
- Ratio > 1: More pathogenic than benign (higher disease burden)
- Ratio = 1: Equal representation
- Ratio < 1: More benign than pathogenic (lower disease burden)

### Variant Type Analysis

Distribution of variant types (e.g., Insertion, Deletion, Substitution) was tabulated and stratified by clinical significance. Percent pathogenic calculated for each type.

### Functional Consequence Analysis

Variants were categorized by functional consequence (missense, frameshift, splice-site, synonymous, etc.) if this annotation was present in ClinVar. Consequences were analyzed stratified by pathogenicity to identify which types associate with disease.

### Allele Frequency Comparison

For variants with AF data available:
- Mean, median, min, max AF for pathogenic variants
- Mean, median, min, max AF for benign variants
- Percent coverage calculated separately for each group

**Expected pattern:** Pathogenic variants typically have lower AF due to purifying selection; benign variants more tolerated at higher frequencies.

---

## Visualization

Publication-quality figures were generated using matplotlib and seaborn:

1. **Gene distribution plot:** Horizontal bar chart showing top N genes by variant count, stratified by pathogenic/benign
2. **Variant type plot:** Bar chart comparing type distributions across clinical significance categories
3. **Pathogenic-to-benign ratio plot:** Horizontal bar chart with reference line at ratio=1
4. **Clinical significance pie chart:** Proportional representation of pathogenic vs benign
5. **Allele frequency histogram:** Overlaid distributions (if AF data sufficient)
6. **Consequence plot:** Bar chart of functional consequence types (if available)

All figures saved at 300 DPI in PNG format for publication quality.

---

## Important Limitations

### 1. Curated Annotations Only
This analysis uses exclusively ClinVar curated community annotations. **It does not represent ground truth.** ClinVar classifications reflect expert consensus and submitter assertions, which may evolve as new evidence emerges.

### 2. No Population-Scale Frequency Data
Analysis does **not** incorporate gnomAD allele frequencies or other population databases. Any allele frequency presented is from ClinVar's limited AF field, which is often sparse. For allele frequency interpretation, users should consult gnomAD directly.

### 3. No Experimental Validation
This analysis does **not** include functional studies, protein structure predictions, molecular dynamics simulations, or in vitro/in vivo validation. Clinical significance is based on literature consensus, not experimental evidence in this study.

### 4. Binary Classification
"Uncertain significance" variants are excluded to enable clear pathogenic vs. benign comparison. However, many variants in ClinVar are genuinely uncertain; excluding them creates a selection bias toward well-characterized variants.

### 5. Gene-Level Mapping Only
Variants are mapped to genes, not specific isoforms, protein domains, or regulatory elements. Multiple isoforms of the same gene may be differentially affected. Domain-level analysis treated as extension.

### 6. Overrepresentation of Disease Associations
ClinVar contains a higher proportion of disease-associated variants than the general population. This is intentional—ClinVar is a clinical database. Frequency statistics should not be interpreted as population-representative.

### 7. Allele Frequency Data Sparserity
AF information in ClinVar is incomplete and often missing, particularly for rare variants and less-studied genes. Conditional AF analysis documents this limitation explicitly.

### 8. No Regulatory Annotation
Analysis does not include regulatory region classification (enhancers, silencers, etc.) as this requires specialized tools and data not integrated here.

---

## Interpretation Framework

### Biological Validity
Results are interpreted in the context of known biology:
- DNA repair pathways are essential for genome stability
- Pathogenic variants in these genes associate with cancer predisposition and other genetic disease
- Patterns observed (e.g., missense enrichment in pathogenic) are consistent with literature

### Clinical Relevance
Findings connect to medical genetics and cancer genomics:
- Genes with high P/B ratios may warrant heightened clinical scrutiny
- Variant type patterns inform variant classification algorithms
- Consequence mapping aligns with variant effect prediction tools (e.g., PolyPhen, SIFT)

### Caveats for Interpretation
- **Not causal:** Associations observed in this descriptive analysis do not imply causation
- **Population-specific:** Results reflect ClinVar, which has demographic and ascertainment bias (more representation of European ancestry, clinically affected individuals)
- **Not a diagnostic tool:** Results should not be used clinically without independent validation

---

## Reproducibility

All analyses are implemented in Python with pandas, matplotlib, and seaborn. Code is modular and documented. The complete workflow is reproducible on a standard laptop:
- **Runtime:** ~5–10 minutes (first run includes ClinVar download)
- **Storage:** <500 MB total
- **RAM:** <200 MB peak

Notebook 01 downloads and preprocesses data. Notebook 02 performs all analyses and generates figures/tables. Code is version-controlled and published to GitHub for transparency.

---

## Future Directions

**Extensions not included in current analysis:**
1. Integrate gnomAD population frequencies for allele frequency context
2. Transcript isoform and domain-level mapping
3. In silico prediction tool integration (PolyPhen, SIFT, Conservation scores)
4. Pathway enrichment analysis across repair mechanisms
5. Machine learning models for pathogenicity prediction
6. Comparative genomics across species orthologs
7. Interactive web visualization (Streamlit/Shiny)

---

## Software and Versions

- **Python:** 3.8+
- **pandas:** 2.0.3
- **numpy:** 1.24.3
- **matplotlib:** 3.7.2
- **seaborn:** 0.12.2
- **scipy:** 1.11.2 (optional, for future statistical tests)

---

## References

1. Landrum, M. J., et al. (2023). ClinVar at NCBI: Improving access and support for clinical variant interpretation. _Human Mutation_, 44(5), 1108–1120.

2. Karczewski, K. J., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. _Nature_, 581(7809), 434–443. (gnomAD)

3. The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. _Nature_, 526(7571), 68–74.

4. Wood, R. D., Mitchell, M., & Lindahl, T. (2005). Human DNA repair genes. _Mutation Research/Fundamental and Molecular Mechanisms of Mutagenesis_, 569(1–2), 15–128.

---

**Analysis date:** April 2026  
**Data version:** Latest ClinVar release at time of analysis  
**Status:** Graduate-level course submission
