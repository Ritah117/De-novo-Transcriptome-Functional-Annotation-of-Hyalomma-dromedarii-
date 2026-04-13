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
#SBATCH --job-name=Hdrom_ORFs
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/amukami/Denovo/trinity_out/trans_%j.log
#SBATCH --error=/nfs/amukami/Denovo/trinity_out/trans_%j.err

# 1. The path to your tools (Updated to your actual tick_tools location)
BIN="/nfs/amukami/tick_tools_nfs/bin"

# 2. Navigate to your current project directory
cd /nfs/amukami/Denovo/trinity_out

# 3. Define the input (Using the Trinity.fasta we just analyzed)
TRINITY_FILE="Trinity.fasta"

echo "--- Starting TransDecoder for H. dromedarii ---"
date

# 4. Run TransDecoder.LongOrfs (Finds potential proteins)
echo "Running TransDecoder.LongOrfs..."
$BIN/TransDecoder.LongOrfs -t $TRINITY_FILE

# 5. Run TransDecoder.Predict (Finalizes the protein sequences)
echo "Running TransDecoder.Predict..."
$BIN/TransDecoder.Predict -t $TRINITY_FILE

date
echo "--- Pipeline Finished ---"
```
# 6. Chemosensory Gene Discovery & Targeted Annotation
Following protein prediction, we implemented a high-stringency bioinformatics screen to isolate specific chemosensory gene families. This targeted approach is critical for understanding the molecular basis of host-seeking in H. dromedarii.

# 7. The "Gold Standard" Database Strategy

We utilized a curated 10-species arthropod sensory database  to perform a reciprocal BLAST search against the 53,646 predicted proteins.

Database Composition: Curated sequences of IRs, GRs, ORs, OBPs, CSPs, and SNMPs.
## Dataset: Arthropod Species of Medical, Veterinary, and Research Importance
The organisms are grouped into three strategic categories:

Medical and Veterinary Vectors – Blood-feeding species responsible for the transmission of human and animal pathogens, including the Brown Ear Tick (East Coast Fever) and a diverse Tsetse fly trio (Trypanosomiasis).

Livestock and Urban Pests – Species such as the Stable Fly and Bed Bug that cause significant economic losses in animal production or represent major urban nuisance challenges.

Biological Model Organisms – Highly characterized species like the Vinegar Fly and Silk Moth, which serve as the gold standard for understanding insect olfaction and pheromone detection. 

### Species and Their Significance

This table summarizes key arthropod species included in the study, their ecological category, and their biological or economic significance.

| Category | Species | Common Name | Ecological/Economic Significance |
|---|---|---|---|
| Acarine Vector | *Rhipicephalus appendiculatus* | Brown Ear Tick | Primary vector of East Coast Fever (Theileriosis) in African cattle. |
| Dipteran Vectors | *Glossina morsitans* | Savannah Tsetse Fly | Major vector of Nagana in cattle and Sleeping Sickness in humans. |
|  | *Glossina fuscipes* | Riverine Tsetse Fly | Responsible for over 90% of human Sleeping Sickness cases in Africa. |
|  | *Glossina brevipalpis* | Forest Tsetse Fly | Large-bodied vector contributing to animal trypanosomiasis transmission. |
|  | *Anopheles gambiae* | Malaria Mosquito | Primary vector of human malaria (*Plasmodium falciparum*) in sub-Saharan Africa. |
|  | *Aedes aegypti* | Yellow Fever Mosquito | Global vector of Zika, Dengue, and Yellow Fever viruses. |
| Parasites | *Stomoxys calcitrans* | Stable Fly | Significant livestock pest causing severe milk and meat production losses. |
|  | *Cimex lectularius* | Bed Bug | Obligate blood-feeder and significant urban nuisance pest. |
| Model Organisms | *Drosophila melanogaster* | Vinegar Fly | Benchmark model for genetics and insect olfaction research. |
|  | *Bombyx mori* | Silk Moth | Model for Lepidopteran pheromone detection and silk production. |



#  Data Sourcing & Proteome Assembly 

Proteomes were retrieved in FASTA format from UniProt (Reference Proteomes) and NCBI RefSeq. To ensure a standardized comparative analysis, these individual proteomes were concatenated into a single master database.

Master File Name: Master_Arthropod_Reference_FINAL.fasta

Total Sequence Count: 234,439 proteins

## Protein Sequence Dataset Summary

This table summarizes the number of protein sequences retrieved for each species and the primary biological database used as the source.

## Protein Sequence Sources by Species

This table summarizes the total number of protein sequences retrieved for each species and the primary biological database used.

| Species | Sequence Count | Primary Source |
|---|---|---|
| *Rhipicephalus appendiculatus* | 56,508 | NCBI RefSeq |
| *Drosophila melanogaster* | 30,802 | UniProt |
| *Stomoxys calcitrans* | 26,015 | NCBI RefSeq |
| *Cimex lectularius* | 24,194 | NCBI RefSeq |
| *Bombyx mori* | 22,510 | UniProt |
| *Aedes aegypti* | 20,643 | UniProt |
| *Anopheles gambiae* | 14,102 | UniProt |
| *Glossina morsitans* | 12,910 | UniProt |
| *Glossina fuscipes* | 13,500 | NCBI/UniProt |
| *Glossina brevipalpis* | 13,200 | NCBI/UniProt |

\*Approximate sequence counts based on combined database records.

#  HMMER Analysis & Gene Discovery
To identify the chemosensory repertoire, we performed a profile-based search using HMMER 3.3.2. 

#Targeted Gene Families

We searched for six distinct families essential for arthropod chemical sensing:

Odorant Receptors (OR)

Ionotropic Receptors (IR)

Gustatory Receptors (GR)

Odorant Binding Proteins (OBP)

Chemosensory Proteins (CSP)

Sensory Neuron Membrane Proteins (SNMP)

#The Search Command

We utilized the hmmsearch tool with a stringent E-value cutoff of 1e-5. The following loop was used to automate the search across all six families:

```bash
# 1. Create the results directory if it doesn't exist
mkdir -p results

# 2. Automated search loop across the Official 10-species Master Reference
for family in OR IR GR OBP CSP SNMP; do
    echo "Processing family: $family"
    
    # Using the full path to HMM profiles identified in our workspace
    hmmsearch --tblout results/${family}_hits_FINAL.tbl \
    --noali -E 1e-5 \
    insect_chemo_project/hmm_profiles/${family}.hmm \
    Master_Arthropod_Reference_FINAL.fasta
done
```

###  Summary Table: Identified Gene Counts

## Comparative Chemosensory Gene Distribution

This table compares the distribution of major chemosensory gene families across selected arthropod species, including ticks, mosquitoes, flies, and model organisms.

| Species | OR | IR | GR | OBP | CSP | SNMP |
|---|---|---|---|---|---|---|
| *R. appendiculatus* (Tick) | 0 | 84 | 46 | 18 | 0 | 6 |
| *G. morsitans* (Savannah Tsetse) | 46 | 20 | 15 | 5 | 5 | 15 |
| *G. fuscipes* (Riverine Tsetse) | 42 | 28 | 15 | 6 | 5 | 15 |
| *G. brevipalpis* (Forest Tsetse) | 41 | 24 | 13 | 6 | 4 | 13 |
| *C. lectularius* (Bedbug) | 65 | 65 | 28 | 12 | 19 | 29 |
| *A. aegypti* (Mosquito) | 85 | 53 | 54 | 6 | 48 | 25 |
| *A. gambiae* (Mosquito) | 80 | 31 | 92 | 7 | 7 | 17 |
| *S. calcitrans* (Stable Fly) | 75 | 51 | 69 | 6 | 13 | 34 |
| *D. melanogaster* (Model Organism) | 64 | 53 | 84 | 9 | 9 | 35 |
| *B. mori* (Silkworm) | 109 | 43 | 26 | 5 | 43 | 35 |


# Data Extraction and Curated Databases
To facilitate downstream analysis, full protein sequences were extracted from the master file using the hit IDs from the HMMER results. This process creates specific FASTA files for each gene family, which are essential for building phylogenetic trees or performing BLAST searches.

#The Automated HMMER Search

This code runs the search for OR, IR, GR, OBP, CSP, and SNMP simultaneously. It uses the curated .hmm profiles against the 10-species consolidated proteome.

```bash
# Define the families we are looking for
for family in OR IR GR OBP CSP SNMP; do
    echo "Starting search for: $family"
    
    # Run hmmsearch with the 1e-5 threshold
    # --tblout saves the results in an easy-to-read table format
    # We point to the specific path identified in the workspace
    hmmsearch --tblout results/${family}_hits_FINAL.tbl \
    --noali -E 1e-5 \
    insect_chemo_project/hmm_profiles/${family}.hmm \
    Master_Arthropod_Reference_FINAL.fasta
done
```
#Counting the Hits (To Build the Comparative Table)

After the search, we use this command to count how many high-confidence sequences were found for each family.

```bash
for file in results/*_FINAL.tbl; do
    printf "$file: "
    grep -v "^#" "$file" | wc -l
done
```
#Extracting the FASTA Sequences

Once the hit IDs were identified in the .tbl files, we pulled the actual protein sequences out of the master file to create curated family-specific databases.

```bash
# Ensure the output directory exists
mkdir -p family_databases

for family in OR IR GR OBP CSP SNMP; do
    echo "Extracting sequences for $family..."
    
    # 1. Get the IDs from the first column of the HMMER table
    awk '!/^#/ {print $1}' results/${family}_hits_FINAL.tbl > temp_ids.txt
    
    # 2. Use the ID list to pull sequences from the Master Reference
    # We use -A 1 to grab the header and the sequence line immediately following it
    grep -Ff temp_ids.txt -A 1 Master_Arthropod_Reference_FINAL.fasta | \
    grep -v "^--" > family_databases/${family}_arthropod.fasta
done

# Clean up temporary file
rm temp_ids.txt
```
# Summary of Curated Databases
After executing the extraction pipeline above, the following family-specific FASTA databases were generated.
##  Sequence Files Summary

This table summarizes the **protein sequence files used for comparative genomics**, including the number of sequences and a brief description of each gene family.

## Chemosensory Sequence Files

The following FASTA files contain curated protein sequences used for comparative analysis of major arthropod chemosensory gene families.

## Chemosensory Gene Family Sequence Files

| File Name | Sequence Count | Description |
|---|---|---|
| OR_arthropod.fasta | 607 | Odorant Receptors: Insect-specific proteins; key for volatile sensing. |
| IR_arthropod.fasta | 452 | Ionotropic Receptors: Evolutionarily ancestral; found in all 10 species. |
| GR_arthropod.fasta | 442 | Gustatory Receptors: Taste (sugar/bitter) and CO₂ sensors. |
| OBP_arthropod.fasta | 80 | Odorant Binding Proteins: Small carriers for hydrophobic odors. |
| CSP_arthropod.fasta | 153 | Chemosensory Proteins: Broadly expressed soluble carriers. |
| SNMP_arthropod.fasta | 224 | Sensory Neuron Membrane Proteins: Essential pheromone co-factors. |

# 8. Target Discovery: The H. dromedarii Chemosensory Repertoire

Using the curated arthropod databases generated, a reciprocal blastp search was performed against the 53,646 predicted proteins of H. dromedarii to identify high-confidence orthologs.

```bash
#!/bin/bash
#SBATCH --job-name=Hdrom_Extract
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=extract_%j.log

# 1. Setup Paths
QUERY_PEP="/nfs/amukami/Denovo/trinity_out/Trinity.fasta.transdecoder.pep"
RESULTS_DIR="/nfs/amukami/Denovo/trinity_out/scoreboard_results"
OUTPUT_DIR="/nfs/amukami/Denovo/trinity_out/extracted_sequences"

mkdir -p $OUTPUT_DIR

# 2. Families to extract
FAMILIES=("IR" "OBP" "SNMP" "GR")

echo "--- Starting Sequence Extraction (Python Version) ---"

for FAM in "${FAMILIES[@]}"; do
    TAB_FILE="$RESULTS_DIR/${FAM}_Hdrom_matches.tab"
    ID_FILE="$OUTPUT_DIR/${FAM}_ids.txt"
    FASTA_OUT="$OUTPUT_DIR/Hdrom_${FAM}_sequences.fasta"

    if [ -f "$TAB_FILE" ]; then
        # Get unique Trinity IDs
        cut -f2 "$TAB_FILE" | sort | uniq > "$ID_FILE"

        # Use Python to extract sequences safely
        python3 -c "
import sys
ids = set(line.strip() for line in open('$ID_FILE'))
with open('$QUERY_PEP') as f:
    keep = False
    for line in f:
        if line.startswith('>'):
            # Grab the ID before any spaces or dots
            header_id = line[1:].split()[0].split('.')[0]
            # Also check full ID in case versioning matters
            full_id = line[1:].split()[0]
            keep = header_id in ids or full_id in ids
        if keep:
            print(line, end='')
        " > "$FASTA_OUT"

        COUNT=$(grep -c ">" "$FASTA_OUT")
        echo "Family $FAM: Extracted $COUNT sequences."
    else
        echo "Family $FAM: No BLAST results found, skipping."
    fi
done

echo "Extraction complete."
```

Results: The H. dromedarii Sensory Scoreboard

| Gene Family                                   | H. dromedarii Sequence Count |
|----------------------------------------------|------------------------------|
| Odorant Binding Proteins (OBPs)              | 64                           |
| Ionotropic Receptors (IRs)                   | 42                           |
| Sensory Neuron Membrane Proteins (SNMPs)     | 17                           |
| Gustatory Receptors (GRs)                    | 1                            |
| Chemosensory Proteins (CSPs)                 | 0                            |
| Odorant Receptors (ORs)                      | 0                            |
# Phylogenetic Tree Construction





