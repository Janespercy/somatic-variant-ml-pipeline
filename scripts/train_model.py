#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
import pickle

def main():
    input_csv = snakemake.input[0]
    output_model = snakemake.output[0]
    
    # Load features
    df = pd.read_csv(input_csv)
    
    # Create synthetic pathogenicity labels for demo
    # In real analysis, these would come from clinical databases
    df['pathogenic'] = np.where(
        (df['gene_role'] == 'oncogene') | (df['gene_role'] == 'tumor_suppressor'), 1, 0
    )
    
    # Prepare features for ML
    le = LabelEncoder()
    feature_cols = ['qual', 'tumor_vaf', 'vaf_ratio', 'tumor_depth', 'is_transition']
    
    # Add encoded categorical features
    df['gene_role_encoded'] = le.fit_transform(df['gene_role'])
    feature_cols.append('gene_role_encoded')
    
    X = df[feature_cols]
    y = df['pathogenic']
    
    # Train model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    # Save model and feature info
    model_data = {
        'model': model,
        'features': feature_cols,
        'label_encoder': le
    }
    
    with open(output_model, 'wb') as f:
        pickle.dump(model_data, f)
    
    print(f"Trained model on {len(df)} variants")
    print(f"Feature importance: {dict(zip(feature_cols, model.feature_importances_))}")

if __name__ == "__main__":
    main()
