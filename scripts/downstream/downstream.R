library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(DropletUtils)
library(AcidPlots)
library(scuttle)
library(scDblFinder)

### LOAD THE DATA ###

# load the gene expression matrix
undiff_mat <- readMM("data/processed-data/undifferentiated/counts_unfiltered/cells_x_genes.mtx")
# load barcode and gene information
undiff_barcodes <- readLines("data/processed-data/undifferentiated/counts_unfiltered/cells_x_genes.barcodes.txt")
undiff_genes <- readLines("data/processed-data/undifferentiated/counts_unfiltered/cells_x_genes.genes.txt")
undiff_gene_names <- readLines("data/processed-data/undifferentiated/counts_unfiltered/cells_x_genes.genes.names.txt")

### Sanity Checks ###

# Check the dimension and class of the gene expression matrix
dim(undiff_mat)
class(undiff_mat)
# check the lengths of the vectors containing information on barcodes, gene ID and names
length(undiff_barcodes)
length(undiff_genes)
length(undiff_gene_names)

### SingleCellExperiment (SCE) Object Creation ###

# Transpose the gene expression matrix so that it becomes a (Gene X Cell matrix) from (Cell X Gene matrix) that it currently is
undiff_mat <- t(undiff_mat)

rownames(undiff_mat) <- undiff_genes # Set row names of the matrix as the ensembl gene IDs
colnames(undiff_mat) <- undiff_barcodes # Set column names of the matrix as the barcodes for each droplet

# Create the SCE object
undiff_sce <- SingleCellExperiment(assays = list(counts = undiff_mat))
undiff_sce$condition <- "Undifferentiated"

# Sanity Check the SCE Object
dim(undiff_mat)
length(undiff_barcodes)
length(undiff_genes)
dim(undiff_sce)

### Cell Calling ###

# Calculate barcode ranks
undiff_ranks <- barcodeRanks(counts(undiff_sce))

plot(undiff_ranks$rank, undiff_ranks$total, log = "xy", xlab = "Barcode rank", ylab = "Total UMI count", main = "Undiff barcode rank")
o <- order(undiff_ranks$rank)
lines(undiff_ranks$rank[o], undiff_ranks$fitted[o])
abline(h = metadata(undiff_ranks)$knee, lty = 2)
abline(h = metadata(undiff_ranks)$inflection, lty = 2)

# store the values as metadata
metadata(undiff_ranks)$knee
metadata(undiff_ranks)$inflection

# Dropping empty droplets
undiff_emptydrops <- emptyDrops(counts(undiff_sce))


sum(undiff_emptydrops$FDR <= 0.01, na.rm = TRUE) # How many barcodes have been called for cell containing droplets with a FDR less than 1% (8194 for undifferentiated library)

# Check Unique Molecular Identifier (UMI) distribution of candidate cells
undiff_called <- undiff_emptydrops$FDR <= 0.01
summary(undiff_called) # Median UMI for candidate cells is 8194

undiff_called <- !is.na(undiff_emptydrops$FDR) &
                undiff_emptydrops$FDR <= 0.01

sum(undiff_called)

# In total we have 8194 droplets out of 903531 droplets that have contains cells statistically

### Quality Control ###

# Subset the original SCE object to only retain droplets that contain cell/s
undiff_sce_called <- undiff_sce[,undiff_called] # using the undiff_called logical vector we are only retaining the columns (barcodes) of the SCE object that have been statistically found to contain cells

# Calculating  cell level QC metrics (total number of genes per cell, total number of UMI counts per cell, mitochondrial gene percentage)
undiff_qc <- perCellQCMetrics(undiff_sce_called)

#undiff_qc$sum
#undiff_qc$detected

# Add the metadata to the SCE object
colData(undiff_sce_called)$total_counts <- undiff_qc$sum
colData(undiff_sce_called)$detected_genes <- undiff_qc$detected

hist(undiff_sce_called$total_counts,breaks = 100,main = "Undifferentiated: total UMI counts per called cell",xlab = "Total UMI counts")
hist(undiff_sce_called$detected_genes, break=100, main = "Undifferentiated: total number of genes per called cell", xlab = "Number of detected genes")

summary(undiff_sce_called$total_counts)
summary(undiff_sce_called$detected_genes)

# Find Mitochondrial gene percentage in the called cells

# Add the gene names from the undiff_gene_names vector we loaded to R earlier to the original SCE object
rowData(undiff_sce)$gene_name <- undiff_gene_names

# Check how many gene names start with MT (13 in this case)
sum(grepl("^MT-",rowData(undiff_sce)$gene_name))

# Create a boolean vector which will be used later to ommit the mitochondrial genes with explanation
is_mito <- grepl("^MT-",rowData(undiff_sce)$gene_name)

undiff_qc <- perCellQCMetrics(undiff_sce_called,subsets=list(mito=is_mito))
colData(undiff_sce_called)$mito_percent <- undiff_qc$subsets_mito_percent

summary(undiff_sce_called$mito_percent)

hist(undiff_sce_called$mito_percent, breaks=100, main="Undifferentiated: Mitochondrial Percentage", xlab="Mitochondrial counts (%)")
plot(undiff_sce_called$total_counts,undiff_sce_called$mito_percent, log = "x", pch = 16, cex = 0.5, xlab = "Total UMI counts", ylab = "Mitochondrial counts (%)", main = "Undifferentiated: UMI Counts vs Mitochondrial Percentage")

# Identify outliers in the QC metrics
undiff_low_umi <- isOutlier(undiff_sce_called$total_counts, log = TRUE, type = "lower", nmads = 3)
undiff_low_genes <- isOutlier(undiff_sce_called$detected_genes, log=TRUE, type="lower", nmads=3)
undiff_high_mito <- isOutlier(undiff_sce_called$mito_percent, log=TRUE, type="lower", nmads=3)

# Store these  information in the called SCE object
colData(undiff_sce_called)$low_umi <- undiff_low_umi
colData(undiff_sce_called)$low_genes <- undiff_low_genes
colData(undiff_sce_called)$high_mito <- undiff_high_mito

# Now create a new SCE object by filtering the original SCE object to only retain cells that meet the QC criteria (i.e. remove cells that have values for the QC metrics more than 3 standard deviations away from the mean of the respective qc metric)

# Create a boolean vector for retaining all cells meeting the QC criteria
undiff_keep <- !undiff_low_umi & !undiff_low_genes & !undiff_high_mito
undiff_sce_filtered <- undiff_sce_called[, undiff_keep]

dim(undiff_sce_filtered) # Filetered sce object contains 8033 called cells

# Detecting droplets containing more than one cells
undiff_sce_filtered <- scDblFinder(undiff_sce_filtered)

table(undiff_sce_filtered$scDblFinder.class)

# Only retain droplets containing single cell
undiff_sce_filtered <- undiff_sce_filtered[, undiff_sce_filtered$scDblFinder.class == "singlet"]

dim(undiff_sce_filtered) # 7677 droplets remain after removing doublets

save.image("data/10X_pipeline_undiff_session.RData")
