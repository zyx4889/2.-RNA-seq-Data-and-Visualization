# 2. RNA-seq Data and Visualization


## RNA-seq Analysis Pipeline

This comprehensive RNA-seq analysis pipeline details steps from raw sequencing data (FASTQ files) to differential expression analysis and pathway enrichment.

---

## Table of Contents
1. [Pre-processing: FASTQ Download and Quantification](#Pre-processing)
2. [Visualization and Differential Expression Analysis](#visualization-and-deseq-analysis)
3. [Pathway Enrichment Analysis](#pathway-enrichment)

---

## 1. Pre-processing: FASTQ Download and Quantification

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

### Automated Bash Script ([0_accession.sh](./script/0_accession.sh))

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

## Step 1.5: Merging Paired-end Reads from Multiple Runs

### Overview

This step is optional and only required if your sequencing samples are spread across multiple runs and need to be combined into single paired-end files per sample. If each sample was sequenced only once, you can skip this step.

### Why Merge Paired-end Reads?

Sequencing facilities often split samples across multiple lanes or runs. Merging these files ensures you have a single set of paired-end FASTQ files per sample for downstream analysis.

### Automated Bash Script ([1.5A_merge.sh](./script/1.5A_merge.sh)) 

Here's an easy-to-use interactive script that guides you step-by-step to merge paired-end FASTQ files:

```bash
#!/bin/bash

echo "Merging paired-end reads for each sample..."

# Prompt user for the directory containing raw FASTQ files
read -p "Enter the working directory (containing your raw FASTQ files): " WORK_DIR

# Prompt user for the desired output directory for merged files
read -p "Enter the output directory for merged files: " MERGE_DIR
mkdir -p ${MERGE_DIR}

# Prompt user for the total number of samples to merge
read -p "Enter the number of samples you want to merge: " NUM_SAMPLES

# Declare an array to store sample names and accession IDs
declare -A SAMPLES

# Get sample names and accession IDs from the user
for i in $(seq 1 $NUM_SAMPLES); do
    read -p "Enter a short, descriptive name for sample $i (e.g., Control_1): " SAMPLE_NAME
    read -p "Enter the accession IDs for sample $i (without _R1/_R2), separated by space (e.g., SRR1234567 SRR1234568): " ACCESSION_IDS
    SAMPLES+=( ["${SAMPLE_NAME}"]="${ACCESSION_IDS}" )
done

# Perform the merging for each sample
for SAMPLE in "${!SAMPLES[@]}"; do
    PAIRED_READS=(${SAMPLES[$SAMPLE]})
    # Merge R1 (forward reads)
    cat ${WORK_DIR}/${PAIRED_READS[0]}_R1.fastq.gz ${WORK_DIR}/${PAIRED_READS[1]}_R1.fastq.gz > ${MERGE_DIR}/${SAMPLE}_1.fastq.gz
    # Merge R2 (reverse reads)
    cat ${WORK_DIR}/${PAIRED_READS[0]}_R2.fastq.gz ${WORK_DIR}/${PAIRED_READS[1]}_R2.fastq.gz > ${MERGE_DIR}/${SAMPLE}_2.fastq.gz
done

echo "All paired-end reads have been successfully merged."
```

### Detailed Tutorial: How to Execute the Script

1. **Place the Script**: Save the script (`1.5A_merge.sh`) into your scripts directory.

2. **Make the Script Executable**:

```bash
chmod +x 1.5A_merge.sh
```

3. **Run the Script**:

Navigate to your script directory and run the script:

```bash
cd scripts
./1.5A_merge.sh
```

### Interactive Inputs During Execution

When executing the script, you'll be asked to provide:

- **Working directory**: The path to your raw FASTQ files (e.g., `/data/raw_fastq`).
- **Output directory**: The path where merged files will be saved (e.g., `/data/merged_fastq`).
- **Number of samples**: The number of unique samples to merge.
- **Sample details**:
  - A clear sample name (e.g., `Treatment_A`).
  - Accession IDs corresponding to that sample (e.g., `SRR1234567 SRR1234568`).

### Example Interaction

```
Enter the working directory (containing your raw FASTQ files): /projects/data/raw_fastq
Enter the output directory for merged files: /projects/data/merged_fastq
Enter the number of samples you want to merge: 2
Enter a short, descriptive name for sample 1 (e.g., Control_1): Control
Enter the accession IDs for sample 1 (without _R1/_R2), separated by space (e.g., SRR1234567 SRR1234568): SRR1111111 SRR1111112
Enter a short, descriptive name for sample 2 (e.g., Treatment_1): Treated
Enter the accession IDs for sample 2 (without _R1/_R2), separated by space (e.g., SRR1234567 SRR1234568): SRR2222222 SRR2222223
```

### Expected Output Structure

After running the script, your merged data directory should contain:

```
📁 merged_fastq/
├── Control_1.fastq.gz
├── Control_2.fastq.gz
├── Treated_1.fastq.gz
├── Treated_2.fastq.gz
└── ... (additional samples)
```

These files are now ready for quality control and further downstream analysis.

---

## Step 2: Quality Control with FastQC

### Overview
Quality control (QC) ensures your sequencing data is high-quality before analysis. FastQC provides visual summaries to evaluate FASTQ file quality.

### Why Perform Quality Control?
QC helps detect issues like low base quality, adapter contamination, and biased nucleotide composition, ensuring accurate downstream analysis.

### Automated Bash Script ([2A_fastqc.sh](./script/2A_fastqc.sh))

This script runs FastQC in an HPC environment using SLURM for resource allocation. Adjust SLURM settings for different computing environments.

```bash
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=08:00:00
#SBATCH --mem=30G
#SBATCH --job-name=genomics
#SBATCH --output=/projects/b1042/AbdulkadirLab/Will/22RV1_CHIP/%j.out
#SBATCH --error=/projects/b1042/AbdulkadirLab/Will/22RV1_CHIP/%j.err

# Capture command-line arguments
OUTPUT=$1
RAW_FASTQ_DIR=$2
KEYWORDS=$3

# Navigate to FASTQ directory
cd ${RAW_FASTQ_DIR}

# Load FastQC module
module load fastqc/0.12.0

# Run FastQC
fastqc -o ${OUTPUT} -t 6 ${KEYWORDS}*.fastq
```

### How to Execute the Script

Save the script ([2A_fastqc.sh](./script/2A_fastqc.sh)) and execute these commands:

1. **Make Script Executable**:
```bash
chmod +x 2A_fastqc.sh
```

2. **Submit to HPC Scheduler (SLURM)**:

Replace placeholders with actual paths:
```bash
sbatch 2A_fastqc.sh <output_directory> <raw_fastq_directory> <keywords>
```
- `<output_directory>`: Directory for saving FastQC reports.
- `<raw_fastq_directory>`: Directory containing FASTQ files.
- `<keywords>` (optional): Filename keywords for selective analysis (use `""` for all files).

### Example Use Case
Analyze all FASTQ files in a directory:
```bash
sbatch 2A_fastqc.sh /projects/b1042/AbdulkadirLab/Will/macrophage_RNAseq/2_fastqc /projects/b1042/AbdulkadirLab/Will/macrophage_RNAseq/Raw_data ""
```

### Expected Output
The output folder will include FastQC reports:
```
📁 2_fastqc/
├── sample1_fastqc.html
├── sample1_fastqc.zip
├── sample2_fastqc.html
├── sample2_fastqc.zip
└── ... (other samples)
```
### Example FastQC Report Interpretation ([Example FastQC HTML Report](./ouput/Example_fastqc.html))

| Metric                           | Expected Results                                        | Example Result | Interpretation and Next Steps                      |
|----------------------------------|---------------------------------------------------------|----------------|----------------------------------------------------|
| Basic Statistics                 | Confirm sequence counts, GC content (~40-60%), read length | ✅ PASS        | Data meets basic expectations                      |
| Per base sequence quality        | Scores consistently >Q30                                | ⚠️ WARNING    | Trim low-quality bases at read ends                |
| Per tile sequence quality        | Uniform quality across tiles                            | ❌ FAIL       | Usually acceptable; investigate if severe          |
| Per sequence quality scores      | High average quality scores (>Q30)                      | ✅ PASS        | Good overall read quality                          |
| Per base sequence content        | Even nucleotide distribution                            | ❌ FAIL       | Often acceptable; check for severe anomalies       |
| Per sequence GC content          | Matches organism's expected GC distribution (~40-60%)   | ✅ PASS        | Appropriate GC content                             |
| Per base N content               | Minimal N content                                       | ✅ PASS        | Good sequencing accuracy                           |
| Sequence Length Distribution     | Uniform read lengths                                    | ⚠️ WARNING    | Usually acceptable; verify trimming strategy       |
| Sequence Duplication Levels      | Low duplication levels (<20%)                           | ❌ FAIL       | Often acceptable for RNA-seq; review if extreme    |
| Overrepresented sequences        | Few or no overrepresented sequences                     | ⚠️ WARNING    | Usually acceptable; remove if adapters/contaminants|
| Adapter Content                  | Low or no detectable adapter sequences                  | ✅ PASS        | Adapter contamination minimal or none              |
| Kmer Content                     | No biased kmer enrichment                               | ❌ FAIL       | Often acceptable; check for contamination or biases|


## Step 3: Trimming and Filtering Reads (Trimmomatic)

### Overview
Trimming and filtering raw FASTQ files remove adapters, low-quality bases, and short sequences, ensuring only high-quality reads are analyzed downstream. This step uses **Trimmomatic**.

### Paired-end Reads Trimming ([3A_trim_filter.sh](./script/3A_trim_filter.sh.txt))
This script specifically trims paired-end sequencing reads, removing adapters and filtering based on quality parameters.

### Key Script Components
- **Adapter removal** (`ILLUMINACLIP`): Removes adapter sequences based on matching criteria.
- **Quality thresholds** (`LEADING:3`, `TRAILING:3`): Removes bases from the start and end of a read if below quality score 3.
- **Sliding window** (`SLIDINGWINDOW:4:15`): Cuts read when average quality within a 4-base sliding window drops below 15.
- **Minimum read length** (`MINLEN:25`): Discards reads shorter than 25 bases after trimming.

### Script (with brief explanations)
```bash
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=08:00:00
#SBATCH --mem=60G

module load trimmomatic/0.39

input_dir=$1
output_dir=$2
keywords=$3

java -jar /hpc/software/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 12 -phred33 \
     ${input_dir}/${keywords}*1.fastq.gz ${input_dir}/${keywords}*2.fastq.gz \      # Input paired reads
     ${output_dir}/forward_paired.fastq.gz ${output_dir}/forward_unpaired.fastq.gz \ # Output forward reads
     ${output_dir}/reverse_paired.fastq.gz ${output_dir}/reverse_unpaired.fastq.gz \ # Output reverse reads
     ILLUMINACLIP:/hpc/software/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \   # Adapter trimming
     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25                               # Quality filtering
```

### Tutorial: How to Execute (Paired-end)
Follow these steps for paired-end trimming:

1. **Make the script executable**:
```bash
chmod +x 3A_trim_filter.sh
```

2. **Submit to HPC Scheduler (SLURM)**:
Replace placeholders with your directories.
```bash
sbatch 3A_trim_filter.sh <input_directory> <output_directory> <keywords>
```
- `<input_directory>`: Directory with your raw FASTQ files.
- `<output_directory>`: Directory for trimmed output.
- `<keywords>`: Common filename prefix or leave blank (`""`) for all files.

### Example Use Case (Paired-end)
```bash
sbatch 3A_trim_filter.sh /projects/b1042/AbdulkadirLab/Will/macrophage_RNAseq/Raw_data /projects/b1042/AbdulkadirLab/Will/macrophage_RNAseq/3_trimmed ""
```

### Single-end Reads Trimming ([3B_trim_filter_single.sh](./script/3B_trim_filter_single.sh))
If your sequencing data is single-end, use the provided single-end trimming script instead.

### How to Execute (Single-end)
```bash
chmod +x 3B_trim_filter_single.sh
sbatch 3B_trim_filter_single.sh <input_directory> <output_directory> <keywords>
```

### Example Use Case (Single-end)
```bash
sbatch 3B_trim_filter_single.sh /projects/b1042/AbdulkadirLab/Will/Songhua_RNAseq/Raw_data /projects/b1042/AbdulkadirLab/Will/Songhua_RNAseq/3_trimmed ""
```

### Expected Output
Trimmed and quality-filtered FASTQ files ready for alignment:
```
📁 3_trimmed/
├── forward_paired.fastq.gz
├── reverse_paired.fastq.gz
├── forward_unpaired.fastq.gz
├── reverse_unpaired.fastq.gz
└── ... (additional samples)
```

These trimmed files will be used for downstream analysis, such as alignment or quantification.

## Step 4: Post-Trimming Quality Control (FastQC)

### Overview
After trimming, it's essential to re-evaluate the quality of your reads. This quick FastQC step confirms the effectiveness of the trimming step.

### Execution Instructions
1. **Make the script executable**:
```bash
chmod +x 4A_fastqctrim.sh
```

2. **Submit to HPC Scheduler (SLURM)**:
```bash
sbatch 4A_fastqctrim.sh <output_folder> <trimmed_fastq_directory> <keywords>
```

### Example Use Case
```bash
sbatch 4A_fastqctrim.sh /projects/p31832/RNAseq/MYCi975/MYCCAP_results/4_fastqc_trim /projects/p31832/RNAseq/MYCi975/MYCCAP_results/3_trimmed ""
```

### Expected Output
Quality control reports for trimmed FASTQ files:
```
📁 4_fastqc_trim/
├── sample1_fastqc.html
├── sample1_fastqc.zip
├── sample2_fastqc.html
├── sample2_fastqc.zip
└── ... (additional samples)
```
Review these HTML reports to confirm trimming success.



## 2. Visualization and Differential Expression Analysis

*Detailed steps and automated scripts will be provided in this section, including QC visualization and differential expression analysis using DESeq2.*

---

## 3. Pathway Enrichment Analysis

*Comprehensive instructions for pathway analysis using pathfindR or GSEA will follow, integrating DE results.*

---

