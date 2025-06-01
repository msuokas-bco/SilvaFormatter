## SilvaFormatter

### Overview

SilvaFormatter is a R package designed to prepare SILVA reference files for microbiome analysis. Primarily it is intended to be used for non-redundant data set.

A makeSilvaFastaNR function is based on taxonomy functions in dada2 R package (see citation). It will select Bacterial and Archaeal taxonomies. Output taxonomy can contain ranks up to the genus or the species level. If species names are included, genera and species names are matched and ambiguous species names are filtered. You can also choose number of random eukaryotic SSU sequences to be included in the output, but eukaryotic sequences are only identified at Kingdom level.

Major difference to the dada2 package is that function writes separate fasta sequence and taxonomy files. Thus, output files can be used e.g. as an alignment reference and taxonomy file helps to parse taxonomy based on alignment hits. Consistent taxonomy is beneficial when using LCA to resolve taxonomy.

If you use tool, please cite the original article.

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). **DADA2: High-resolution sample inference from Illumina amplicon data**. *Nature Methods*, **13**, 581–583. <https://doi.org/10.1038/nmeth.3869>

### 

Installation

You can install SilvaFormatter directly from GitHub with devtools:

``` r
# Install devtools if you don't have it already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install SilvaFormatter
devtools::install_github("msuokas-bco/SilvaFormatter")
```

### Dependencies

SilvaFormatter requires the following R packages: - Biostrings - dada2 (required for species level processing)

### Usage

The main function provided by this package is `makeFastaSilvaNR()`:

``` r
library(SilvaFormatter)

makeFastaSilvaNR(
  fin = "path/to/silva_sequences.fasta",
  ftax = "path/to/silva_taxonomy.tsv",
  fout_fasta = "path/to/output_sequences.fasta.gz",
  fout_taxonomy = "path/to/output_taxonomy.tsv",
  include.species = TRUE,
  compress = TRUE,
  n_euk = 100
)
```

### Function Parameters

| Parameter | Description |
|---------------------------------|---------------------------------------|
| `fin` | Path to the input SILVA FASTA file |
| `ftax` | Path to the SILVA taxonomy file with valid taxon strings and levels |
| `fout_fasta` | Path to the output FASTA file |
| `fout_taxonomy` | Path to the output taxonomy TSV file |
| `include.species` | Logical; whether to include valid species names in taxonomy (default: TRUE) |
| `compress` | Logical; whether to compress the FASTA output using gzip (default: FALSE) |
| `n_euk` | Integer; number of Eukaryota sequences to include in the final output (default: 100) |

### What the Function Does

The `makeFastaSilvaNR()` function performs the following operations:

1.  Reads a SILVA FASTA file containing sequences
2.  Extracts and parses the taxonomic annotations from the FASTA headers
3.  Validates the taxonomic lineages against a reference taxonomy file
4.  Filters out invalid or problematic annotations (e.g., "Uncultured")
5.  Properly formats species names when `include.species = TRUE`
6.  Randomly selects a specified number of Eukaryota sequences to include
7.  Outputs a cleaned FASTA file (optionally compressed) and a taxonomy TSV file

### Output Files

1.  **FASTA file**: Contains the selected sequences with cleaned headers
2.  **Taxonomy TSV file**: A tab-separated file with two columns:
    -   ID: The sequence identifier
    -   Taxonomy: The cleaned taxonomic lineage
