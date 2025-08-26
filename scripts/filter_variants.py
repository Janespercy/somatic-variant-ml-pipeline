#!/usr/bin/env python3

import sys

def parse_vcf_line(line):
    """Parse a VCF data line and calculate VAF"""
    fields = line.strip().split('\t')
    chrom, pos, id_col, ref, alt, qual, filter_col, info, format_col = fields[:9]
    normal, tumor = fields[9:11]
    
    # Parse tumor sample data
    tumor_data = dict(zip(format_col.split(':'), tumor.split(':')))
    normal_data = dict(zip(format_col.split(':'), normal.split(':')))
    
    # Calculate VAF for tumor and normal
    tumor_ad = list(map(int, tumor_data['AD'].split(',')))
    normal_ad = list(map(int, normal_data['AD'].split(',')))
    
    tumor_vaf = tumor_ad[1] / sum(tumor_ad) if sum(tumor_ad) > 0 else 0
    normal_vaf = normal_ad[1] / sum(normal_ad) if sum(normal_ad) > 0 else 0
    
    return {
        'line': line,
        'qual': float(qual),
        'tumor_vaf': tumor_vaf,
        'normal_vaf': normal_vaf,
        'tumor_depth': sum(tumor_ad),
        'normal_depth': sum(normal_ad)
    }

def main():
    input_vcf = snakemake.input[0]
    output_vcf = snakemake.output[0]
    
    # Get filtering parameters from config
    tumor_min_vaf = snakemake.params.tumor_min_vaf
    normal_max_vaf = snakemake.params.normal_max_vaf
    min_depth = snakemake.params.min_depth
    
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            # Copy header lines
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Parse and filter variant lines
            variant = parse_vcf_line(line)
            
            # Apply filters
            if (variant['tumor_vaf'] >= tumor_min_vaf and 
                variant['normal_vaf'] <= normal_max_vaf and
                variant['tumor_depth'] >= min_depth and
                variant['normal_depth'] >= min_depth):
                outfile.write(variant['line'])

if __name__ == "__main__":
    main()
