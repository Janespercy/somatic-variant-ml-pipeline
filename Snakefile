# Somatic Variant ML Pipeline
import os
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Define all final outputs (this tells Snakemake what to build)
rule all:
    input:
        # Database downloads
        config["data"]["databases"]["cosmic"]["local_paths"]["census"],
        config["data"]["databases"]["clinvar"]["local_paths"]["vcf"],
        config["data"]["databases"]["oncokb"]["local_paths"]["actionable"],
        # Pipeline outputs
        "results/04_models/pathogenicity_model.pkl"

# Rule 1: Download COSMIC Cancer Gene Census
rule download_cosmic_census:
    output:
        config["data"]["databases"]["cosmic"]["local_paths"]["census"]
    params:
        url = config["data"]["databases"]["cosmic"]["urls"]["census"]
    shell:
        """
        mkdir -p $(dirname {output})
        curl -L -o {output} {params.url}
        """

# Rule 2: Download ClinVar VCF
rule download_clinvar:
    output:
        config["data"]["databases"]["clinvar"]["local_paths"]["vcf"]
    params:
        url = config["data"]["databases"]["clinvar"]["urls"]["vcf"]
    shell:
        """
        mkdir -p $(dirname {output})
        curl -L -o {output} {params.url}
        """

# Rule 3: Download OncoKB actionable variants
rule download_oncokb:
    output:
        config["data"]["databases"]["oncokb"]["local_paths"]["actionable"]
    params:
        url = config["data"]["databases"]["oncokb"]["urls"]["actionable"]
    shell:
        """
        mkdir -p $(dirname {output})
        curl -L -o {output} {params.url}
        """

# Rule 4: Filter variants by VAF and quality
rule filter_variants:
    input:
        vcf = "data/input/test_somatic.vcf"
    output:
        vcf = "results/01_filtered/filtered_variants.vcf"
    params:
        tumor_min_vaf = config["filtering"]["quality"]["tumor"]["min_vaf"],
        normal_max_vaf = config["filtering"]["quality"]["normal"]["max_vaf"],
        min_depth = config["filtering"]["quality"]["min_depth"]
    script:
        "scripts/filter_variants.py"

# Rule 5: Annotate variants with COSMIC data
rule annotate_cosmic:
    input:
        vcf = "results/01_filtered/filtered_variants.vcf",
        cosmic = "data/databases/cosmic/census.csv"
    output:
        vcf = "results/02_annotated/cosmic_annotated.vcf"
    conda:
        "envs/annotation.yaml"
    script:
        "scripts/annotate_cosmic.py"

# Rule 6: Extract features for ML
rule extract_features:
    input:
        vcf = "results/02_annotated/cosmic_annotated.vcf"
    output:
        features = "results/03_features/variant_features.csv"
    script:
        "scripts/extract_features.py"

# Rule 7: Train ML model
rule train_model:
    input:
        features = "results/03_features/variant_features.csv"
    output:
        model = "results/04_models/pathogenicity_model.pkl"
    params:
        algorithm = config["ml_models"]["algorithms"][0]  # Use first algorithm
    script:
        "scripts/train_model.py"