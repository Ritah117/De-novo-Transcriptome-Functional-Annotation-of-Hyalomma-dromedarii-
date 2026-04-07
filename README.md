# De-novo-Transcriptome-Functional-Annotation-of-Hyalomma-dromedarii-
This project utilizes TransDecoder and HMMER to identify and characterize key chemosensory gene families, including Odorant (ORs), Gustatory (GRs), and Ionotropic (IRs) Receptors, alongside Odorant Binding Proteins (OBPs) and Chemosensory Proteins (CSPs).
#  1. Pre-Assembly Quality Control & Data Integrity
Before functional annotation, raw sequencing reads were migrated to the High-Performance Computing (HPC) cluster and verified for structural integrity.

#Data Integrity Check

Because chemosensory genes (ORs, GRs, IRs) are often expressed at low levels, ensuring file integrity is critical. A single corrupted byte could cause an assembler to miss a low-abundance Odorant Receptor.

```bash
# Verify integrity of compressed fastq files
gzip -t /nfs/amukami/raw_data/*.fastq.gz
```
# 2. Sequence Cleaning with fastp
Objective: Removal of Illumina adapters and low-quality bases ($Q < 20$).

Structural Complexity: Chemosensory receptors are 7-transmembrane proteins; high-quality reads are essential to resolve these complex domains.

Chimera Prevention: Trimming prevents the formation of artificial hybrids between tick DNA and adapters.Read Smoothing: Ensures the assembler generates full-length receptors rather than fragmented sequences.
```bash
#!/bin/bash
#SBATCH --job-name=fastp_batch
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/nfs/amukami/Denovo/fastp_173.log
#SBATCH --error=/nfs/amukami/Denovo/fastp_173.err

# --- 1. Point to the NEW NFS Environment ---
FASTP_BIN="/nfs/amukami/tick_tools_nfs/bin/fastp"

# --- 2. Define the NFS Data Paths ---
BASE_DIR="/nfs/amukami/Denovo"
RAW_DIR="${BASE_DIR}/raw_data"
CLEAN_DIR="${BASE_DIR}/cleaned_data"
REPORT_DIR="${BASE_DIR}/fastp_reports"

# --- 3. Ensure Output Folders Exist ---
mkdir -p $CLEAN_DIR
mkdir -p $REPORT_DIR

# --- 4. Run the Loop ---
for SAMPLE in Hdrom1 Hdrom2 Hdrom3; do
    echo "Starting processing for: ${SAMPLE}"

    $FASTP_BIN \
      -i ${RAW_DIR}/${SAMPLE}_1.fastq.gz \
      -I ${RAW_DIR}/${SAMPLE}_2.fastq.gz \
      -o ${CLEAN_DIR}/${SAMPLE}_clean_R1.fastq.gz \
      -O ${CLEAN_DIR}/${SAMPLE}_clean_R2.fastq.gz \
      -h ${REPORT_DIR}/${SAMPLE}_report.html \
      -j ${REPORT_DIR}/${SAMPLE}_report.json \
      --thread 8

    echo "Finished processing for: ${SAMPLE}"
done
```
# 3. Transcriptome Assembly (Trinity)
The cleaned reads were assembled de novo to reconstruct the full H. dromedarii message.

Strategy: Combined multi-replicate reads to maximize sensitivity for rare sensory transcripts.

Result: Generation of a comprehensive transcriptome representing the tick's metabolic and sensory state.
```bash
#!/bin/bash
#SBATCH --job-name=Tick_Trinity_NFS
#SBATCH --partition=debug
#SBATCH --output=/nfs/amukami/Denovo/trinity_%j.log
#SBATCH --error=/nfs/amukami/Denovo/trinity_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=48:00:00

# --- 1. SETUP ENVIRONMENT ---
source /etc/profile.d/modules.sh
export LC_ALL=C
export LANG=C

module purge
module load samtools/1.14
module load trinity/v2.12.0

# --- 2. RUN TRINITY FROM NFS ---
# Note: We are using the /nfs/amukami/Tick_Reads/ path now!
Trinity --seqType fq \
        --max_memory 55G \
        --left /nfs/amukami/Tick_Reads/Hdrom1_clean_R1.fastq.gz,/nfs/amukami/Tick_$        --right /nfs/amukami/Tick_Reads/Hdrom1_clean_R2.fastq.gz,/nfs/amukami/Tick$        --CPU 32 \
        --normalize_reads \
        --output /nfs/amukami/trinity_out

```
# 4. Assembly Quality Control (BUSCO)
The completeness of the H. dromedarii assembly was evaluated against the Arachnida lineage to ensure we haven't missed critical conserved genes.
```bash
busco -i Hdrom_new.fasta \
      -m transcriptome \
      -l arachnida_odb10 \
      -o Hdrom_Busco_Report \
      --cpu 16
```
# 5. Protein Prediction & ORF Identification (TransDecoder)

To identify the actual protein-coding potential of the H. dromedarii transcriptome, we transitioned from nucleotide sequences to a functional proteome using TransDecoder v5.5.0.

The Prediction Workflow
Identification of Longest ORFs: Extracted potential Open Reading Frames with a minimum length of 100 amino acids.

Candidate Prediction: Predicted likely coding sequences based on hexameric frequencies and homology scores.

SLURM Execution Script (transdecoder_step1.sh)
```bash
#!/bin/bash
#SBATCH --job-name=Hdrom_TransDecoder
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --output=longorfs_%j.log

# Path to verified environment binaries
BIN="/home/amukami/.conda/envs/tick_tools/bin"
cd /nfs/amukami/Annotation

# Step 1: Extract Longest ORFs (>100aa)
$BIN/TransDecoder.LongOrfs -t Hdrom_new.fasta

# Step 2: Predict Coding Regions
$BIN/TransDecoder.Predict -t Hdrom_new.fasta
```





