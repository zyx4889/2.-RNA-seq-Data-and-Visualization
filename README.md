# 2. RNA-seq Data and Visualization

---

---

## RNA-seq Analysis Pipeline

This comprehensive RNA-seq analysis pipeline details steps from raw sequencing data (FASTQ files) to differential expression analysis and pathway enrichment.

---

## Table of Contents
1. [Pre-processing: FASTQ Download and Quantification](#pre-processing)
2. [Visualization and Differential Expression Analysis](#visualization-and-deseq-analysis)
3. [Pathway Enrichment Analysis](#pathway-enrichment)

---

## 1. Pre-processing

### Overview
This section outlines how to efficiently download raw FASTQ files from public repositories using SRA Toolkit's `fasterq-dump`.

### Step 1: FASTQ Download

This step is optional and only required if you are analyzing publicly available data. If you already have your FASTQ files available, you can skip directly to the next step.

We utilize `fasterq-dump` from the `sratoolkit` to download FASTQ files based on accession numbers listed in a text file (`accession.txt`). The process is parallelized to expedite the download.

### Creating the Accession File

Create a simple text file (`accession.txt`) listing all accession numbers (one per line) you wish to download from the Sequence Read Archive (SRA). For example:

```
SRR1234567
SRR1234568
SRR1234569
```

Save this file in your working directory.

### Automated Bash Script (`0_accession.sh`)

This script is designed for an HPC environment (High-Performance Computing) using SLURM. The SBATCH directives (`#SBATCH`) specify the computational resources needed. If you are using another system, these directives may need to be modified accordingly.

```bash
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=08:00:00
#SBATCH --mem=50G
#SBATCH --job-name=genomics
#SBATCH --output=/projects/p31832/RNAseq/Commands/%j.out
#SBATCH --error=/projects/p31832/RNAseq/Commands/%j.err

# Load the sratoolkit module
module load sratoolkit/3.0.0

# User-defined directories (change these paths as needed)
SCRATCH_DIR=/path/to/your/scratch_directory
DATA_DIR=/path/to/your/data_directory
ACCESSIONS_FILE=/path/to/your/accession.txt

# Ensure the directories exist (script creates them if not)
mkdir -p $SCRATCH_DIR
mkdir -p $DATA_DIR

# Function to download SRA data
download_sra() {
    ACCESSION=$1
    fasterq-dump -p -t $SCRATCH_DIR -O $DATA_DIR $ACCESSION
}

# Export the function and variables for GNU Parallel
export -f download_sra
export SCRATCH_DIR
export DATA_DIR

# Download FASTQ files in parallel
cat $ACCESSIONS_FILE | parallel download_sra
```

### How to Customize the Script

- Replace `/path/to/your/scratch_directory` and `/path/to/your/data_directory` with paths specific to your HPC environment.
- Update `/path/to/your/accession.txt` to the correct path where your `accession.txt` file is saved.

### Execution Steps

1. **Make the script executable**:

```bash
chmod +x 0_accession.sh
```

2. **Submit the script to your HPC job scheduler (SLURM)**:

```bash
sbatch 0_accession.sh
```

3. **Verify the downloads**:

```bash
ls $DATA_DIR
```

### Expected Output

After execution, your data directory should contain downloaded FASTQ files:

```plaintext
📁 rawdata/
├── SRR1234567_1.fastq
├── SRR1234567_2.fastq
├── SRR1234568_1.fastq
├── SRR1234568_2.fastq
└── ... (other samples)
```

---

## 2. Visualization and Differential Expression Analysis

*Detailed steps and automated scripts will be provided in this section, including QC visualization and differential expression analysis using DESeq2.*

---

## 3. Pathway Enrichment Analysis

*Comprehensive instructions for pathway analysis using pathfindR or GSEA will follow, integrating DE results.*

---

