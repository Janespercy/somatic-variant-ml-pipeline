FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget curl git \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip first
RUN pip install --upgrade pip

# Install Python packages separately
RUN pip install snakemake
RUN pip install pyyaml requests pandas numpy scikit-learn matplotlib seaborn

# Copy pipeline files
COPY . .

# Create data directories
RUN mkdir -p data/databases/{cosmic,clinvar,oncokb} \
    data/reference \
    results/{01_filtered,02_annotated,03_features,04_models,05_predictions,06_reports}

# Default command
CMD ["snakemake", "--help"]