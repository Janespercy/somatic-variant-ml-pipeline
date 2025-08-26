#!/usr/bin/env python3

import pandas as pd
import numpy as np

def extract_variant_features(vcf_line):
    """Extract ML features from annotated VCF line"""
    fields = vcf_line.strip().split('\t')
    chrom, pos, ref, alt, qual, info = fields[0], fields[1], fields[3], fields[4], fields[5], fields[7]
    normal, tumor = fields[9:11]
    
    # Parse sample data
    format_fields = fields[8].split(':')
    tumor_data = dict(zip(format_fields, tumor.split(':')))
    normal_data = dict(zip(format_fields, normal.split(':')))
    
    # Calculate VAFs
    tumor_ad = list(map(int, tumor_data['AD'].split(',')))
    normal_ad = list(map(int, normal_data['AD'].split(',')))
    
    tumor_vaf = tumor_ad[1] / sum(tumor_ad)
    normal_vaf = normal_ad[1] / sum(normal_ad) if sum(normal_ad) > 0 else 0
    
    # Extract gene annotation
    gene = "UNKNOWN"
    gene_role = "UNKNOWN"
    for field in info.split(';'):
        if field.startswith('GENE='):
            gene = field.split('=')[1]
        elif field.startswith('GENE_ROLE='):
            gene_role = field.split('=')[1]
    
    return {
        'chrom': chrom,
        'pos': int(pos),
        'ref': ref,
        'alt': alt,
        'qual': float(qual),
        'tumor_vaf': tumor_vaf,
        'normal_vaf': normal_vaf,
        'tumor_depth': sum(tumor_ad),
        'normal_depth': sum(normal_ad),
        'gene': gene,
        'gene_role': gene_role,
        'is_transition': (ref + alt) in ['AG', 'GA', 'CT', 'TC'],
        'vaf_ratio': tumor_vaf / max(normal_vaf, 0.001)  # Tumor/normal VAF ratio
    }

def main():
    input_vcf = snakemake.input[0]
    output_csv = snakemake.output[0]
    
    features = []
    with open(input_vcf, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                features.append(extract_variant_features(line))
    
    # Create feature DataFrame
    df = pd.DataFrame(features)
    df.to_csv(output_csv, index=False)
    
    print(f"Extracted {len(df)} variant features")

if __name__ == "__main__":
    main()
