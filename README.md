# Somatic Variant ML Pipeline

A machine learning pipeline for predicting pathogenicity of somatic variants in cancer genomics. Implements VAF-based filtering, gene annotation, and Random Forest classification.

## Features

- **Automated database integration** (COSMIC, ClinVar, OncoKB)
- **Somatic variant filtering** using biologically-informed VAF thresholds
- **Gene annotation** with cancer gene classification
- **ML feature extraction** from genomic and clinical data
- **Pathogenicity prediction** using Random Forest

## Installation

```bash
# Clone repository
git clone https://github.com/Janespercy/somatic-variant-ml-pipeline.git
cd somatic-variant-ml-pipeline

# Create conda environment
conda env create -f envs/annotation.yaml
conda activate annotation

# Install additional dependencies
pip install -r requirements.txt
