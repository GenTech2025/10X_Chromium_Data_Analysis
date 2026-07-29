# Project Timeline
---

#### 14th July 2026
- Identify dataset and download raw fastq files from ENA.
- Learn how to adapt fastqc to 10X chromium dataset (just make sure to ignore the read 1 results as read 1 is not the actual sequence of interest)
- Dataset identified: PRJNA722257 (10X Chromium)
- Reason for not using this dataset: coupled fastq files stored in ENA (i.e. there are not R_1.fastq and R_2.fastq for each sample but rather a single R.fastq file for the sample)

#### 19th July 2026
- Identify new dataset: PRJNA722257
- Perform quality control using FastQC and MultiQC
- Build Kallisto index for mouse genome: realised the sample was from single nuclei RNA sequencing and not single cell RNA sequencing meaning the system RAM (32GB DDR5) was not sufficient enough to build index from scratch. Confirmed the system RAM was full along with the swap space by monitoring using htop

#### 29th July 2026
- Identify and download new dataset: PRJNA815362 (ensured it is a single cell RNA sequencing dataset)