#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load and visualize pipeline results
features_df = pd.read_csv('results/03_features/variant_features.csv')

# Create VAF distribution plot
plt.figure(figsize=(10, 6))
plt.subplot(1, 2, 1)
plt.scatter(features_df['normal_vaf'], features_df['tumor_vaf'], 
           c=['red' if role in ['oncogene', 'tumor_suppressor'] else 'blue' 
             for role in features_df['gene_role']])
plt.xlabel('Normal VAF')
plt.ylabel('Tumor VAF')
plt.title('Somatic Variant VAF Distribution')

# Feature importance plot
plt.subplot(1, 2, 2)
importance = [0.231, 0.192, 0.154, 0.154, 0.154, 0.115]
features = ['Gene Role', 'Quality', 'VAF Ratio', 'Tumor Depth', 'Transition', 'Tumor VAF']
plt.barh(features, importance)
plt.title('ML Feature Importance')
plt.tight_layout()
plt.savefig('examples/pipeline_results.png', dpi=150, bbox_inches='tight')
print("Created visualization: examples/pipeline_results.png")

