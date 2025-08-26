# Somatic Variant ML Pipeline

A machine learning pipeline for predicting pathogenicity of somatic variants in cancer genomics.

## Features
- Automated database downloading (COSMIC, ClinVar, OncoKB)
- VAF-based somatic variant filtering
- Gene annotation and biological context
- ML feature extraction from genomic data
- Random Forest pathogenicity prediction

## Usage
```bash
# Run complete pipeline
snakemake -j 1

# Run individual steps
snakemake filter_variants -j 1
snakemake train_model -j 1
