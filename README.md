# Somatic Variant ML Pipeline

A machine learning pipeline for predicting pathogenicity of somatic variants in cancer genomics.

## Features
- Automated database downloading (COSMIC, ClinVar, OncoKB)
- VAF-based somatic variant filtering
- Gene annotation and biological context
- ML feature extraction from genomic data
- Random Forest pathogenicity prediction

## Installation


# Clone repository
- git clone https://github.com/Janespercy/somatic-variant-ml-pipeline.git
- cd somatic-variant-ml-pipeline

# Create conda environment

- conda env create -f envs/annotation.yaml
- conda activate annotation

# Install additional dependencies
pip install -r requirements.txt

## Usage

# Run complete pipeline
snakemake -j 1

# Run individual steps
- snakemake filter_variants -j 1
- snakemake train_model -j 1

## Pipeline Overview

1. Database Download - Retrieves COSMIC, ClinVar, OncoKB databases
2. Variant Filtering - Applies tumor VAF ≥5%, normal VAF ≤2% thresholds
3. Gene Annotation - Maps variants to cancer genes and functional roles
4. Feature Extraction - Generates ML features (VAF ratios, mutation types, gene context)
5. Model Training - Trains Random Forest classifier for pathogenicity prediction

## Biological Rationale

The pipeline addresses key challenges in somatic variant interpretation:

- VAF filtering distinguishes true somatic mutations from germline contamination
- Gene context prioritizes variants in known cancer genes
- Feature engineering captures mutation patterns relevant to pathogenicity

## Requirements

- Python 3.9+
- Snakemake ≥7.0
- pandas, scikit-learn, numpy EOF
