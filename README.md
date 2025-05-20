## SilvaFormatter

#### Overview

SilvaFormatter is a R package designed to prepare SILVA reference files for microbiome analysis. Primarily it is intended to be used for non-redundant data set.

A makeSilvaFastaNR function is based on dada2 R package. It will select Bacterial and Archaea taxonomies. Result taxonomy can contain ranks either to genus or species level. If species level is chosen, it will make sure that genus and species names are in agreement. Additionally, it will remove uncultured and unidentified taxa names. Finally, you can choose number of random eukaryotic sequences to be included in the files. Their identity is resolved only at Kingdom level.

Major difference to the dada2 package version is that function writes separate fasta sequence and taxonomy files. Sequences can be used e.g. as an alignment reference and taxonomy file helps to parse taxonomy based on mapped reads. Consistent taxonomy is beneficial when identity is predicted by the lowest common ancestor approach. That is case, when the sequence is mapped to several sequences in the reference file.

If you use tool, please cite the original article.

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). **DADA2: High-resolution sample inference from Illumina amplicon data**. *Nature Methods*, **13**, 581–583. <https://doi.org/10.1038/nmeth.3869>

#### 
Installation

You can install SilvaFormatter directly from GitHub:

``` r
# Install devtools if you don't have it already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install SilvaFormatter
devtools::install_github("msuokas-bco/SilvaFormatter")
```

#### Dependencies

SilvaFormatter requires the following R packages: - Biostrings - ShortRead - dada2

#### Usage

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

#### Function Parameters

| Parameter | Description |
|---------------------------------|---------------------------------------|
| `fin` | Path to the input SILVA FASTA file |
| `ftax` | Path to the SILVA taxonomy file with valid taxon strings and levels |
| `fout_fasta` | Path to the output FASTA file |
| `fout_taxonomy` | Path to the output taxonomy TSV file |
| `include.species` | Logical; whether to include valid species names in taxonomy (default: TRUE) |
| `compress` | Logical; whether to compress the FASTA output using gzip (default: TRUE) |
| `n_euk` | Integer; number of Eukaryota sequences to include in the final output (default: 100) |

#### What the Function Does

The `makeFastaSilvaNR()` function performs the following operations:

1.  Reads a SILVA FASTA file containing sequences
2.  Extracts and parses the taxonomic annotations from the FASTA headers
3.  Validates the taxonomic lineages against a reference taxonomy file
4.  Filters out invalid or problematic annotations (e.g., "Uncultured")
5.  Properly formats species names when `include.species = TRUE`
6.  Randomly selects a specified number of Eukaryota sequences to include
7.  Outputs a cleaned FASTA file (optionally compressed) and a taxonomy TSV file

#### Output Files

1.  **FASTA file**: Contains the selected sequences with cleaned headers
2.  **Taxonomy TSV file**: A tab-separated file with two columns:
    -   ID: The sequence identifier
    -   Taxonomy: The cleaned taxonomic lineage

#### Example Workflow

``` r
library(SilvaFormatter)

# Process the SILVA 138.1 SSU NR99 reference
makeFastaSilvaNR(
  fin = "SILVA_138.2_SSURef_NR99_tax_silva.fasta",
  ftax = "SILVA_138.2_taxonomy.tsv",
  fout_fasta = "silva_138.2_processed.fasta.gz",
  fout_taxonomy = "silva_138.2_taxonomy_processed.tsv",
  include.species = TRUE,
  compress = TRUE,
  n_euk = 100
)
```
