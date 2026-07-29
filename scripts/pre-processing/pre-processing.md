# Different steps involved in pre-processing of 10X Chromium Single Cell RNA Seq Data

> Note all the bash commands/scripts should be executed from the directory this file is in originally

### Download the dataset (PRJNA815362)

```bash
# Note this command should be executed in the data/ directory
# Use aria2 for downloading the compressed fastq files
aria2c \
  -i fastq_urls.txt \
  -x 2 \
  -s 2 \
  -j 4 \
  -c \
  --max-tries=10 \
  --retry-wait=5

```

### Quality Control

```bash
# Run fastqc on the samples
mkdir -p ../../results/fastqc_results
fastqc -t 4 -o ../../results/fastqc_results/ ../../data/*.fastq.gz

# Run multiqc to summarize the fastqc results
mkdir -p ../../results/multiqc_results
multiqc -o ../../results/multiqc_results/ ../../results/fastqc_results/
```

### Download pre-made kallisto index

```bash

```

### Generate Kallisto index reference from scratch

```bash
# create a directory to store the index
mkdir -p ../../data/index

# download the fasta and gtf files for mus musculus
gget ref -w dna,gtf mus_musculus

wget -O ../../data/index/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz "http://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

wget -O ../../data/index/Mus_musculus.GRCm39.116.gtf.gz "http://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz"

gunzip ../../data/index/*.gz
```
### Create the reference index

```bash
# create a kallisto index for mouse genome
kb ref \
  --workflow=nucleus \
  -i ../../data/index/index.idx \
  -g ../../data/index/t2g.txt \
  -c1 ../../data/index/cdna_t2c.txt \
  -c2 ../../data/index/intron_t2c.txt \
  -f1 ../../data/index/cdna.fasta \
  -f2 ../../data/index/intron.fasta \
  ../../data/index/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  ../../data/index/Mus_musculus.GRCm39.116.gtf

```

### Pseudoalignment and Quantification using Kallisto

```bash
# create a directory to store the quantified counts
mkdir -p ../../data/processed-data

# use kb count to pseudoalign and quantify the counts from fastq files
kb count \
  --workflow=nucleus \
  -t 8 \
  -i ../../data/index/index.idx \
  -g ../../data/index/t2g.txt \
  -c1 ../../data/index/cdna_t2c.txt \
  -c2 ../../data/index/intron_t2c.txt \
  -x 10xv3 \
  -o ../../data/counts/GSM9119989 \
  ../../data/fastq/GSM9119989_R1.fastq.gz \
  ../../data/fastq/GSM9119989_R2.fastq.gz

```