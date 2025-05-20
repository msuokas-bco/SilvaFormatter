#' Create SILVA-style FASTA and taxonomy files
#'
#' Reads a SILVA FASTA file, extracts taxonomic annotations, cleans and filters them,
#' and outputs a new FASTA and taxonomy TSV file. Optionally includes a set number of Eukaryota sequences
#' and filters for valid species-level annotations.
#'
#' @param fin Path to the input SILVA FASTA file (RNA sequences).
#' @param ftax Path to the SILVA taxonomy file with valid taxon strings and levels.
#' @param fout_fasta Path to the output FASTA file.
#' @param fout_taxonomy Path to the output taxonomy TSV file.
#' @param include.species Logical; whether to include valid species names in taxonomy.
#' @param compress Logical; whether to compress the FASTA output using gzip.
#' @param n_euk Integer; number of Eukaryota sequences to include in the final output.
#'
#' @return Writes a FASTA and a taxonomy file to disk. Returns nothing.
#' @export
#'
#' @examples
#' \dontrun{
#' makeFastaSilvaNR("SILVA.fasta", "SILVA_taxonomy.tsv",
#'                  "output.fasta.gz", "output_taxonomy.tsv",
#'                  include.species = TRUE, compress = TRUE, n_euk = 10000)
#'}
makeFastaSilvaNR <- function(fin, ftax, fout_fasta, fout_taxonomy,
                             include.species = TRUE,
                             compress = TRUE,
                             n_euk = 100) {
  # Read RNA sequences from FASTA file and convert to DNA
  xset <- DNAStringSet(readRNAStringSet(fin, format = "fasta"))

  # Extract the description lines (these are stored in 'names(xset)')
  descriptions <- names(xset)

  # Extract sequence IDs (first word before the space)
  ref_ids <- sapply(strsplit(descriptions, "\\s"), `[`, 1)

  # Check for duplicated IDs
  if (any(duplicated(ref_ids))) stop("Duplicated sequence IDs detected.")

  # Store the full descriptions with the sequence IDs as names
  taxl <- descriptions
  names(taxl) <- ref_ids

  # Rename sequences with their IDs for easier indexing
  names(xset) <- ref_ids

  # Extract taxonomic information (everything after the ID in the description)
  taxl <- gsub("^[A-Za-z0-9.]+\\s", "", taxl)
  taxl <- gsub(";YM;", ";", taxl)
  taxa <- strsplit(taxl, ";")

  # Read taxonomy reference data
  silva.taxa <- read.table(ftax, sep = "\t", col.names = c("Taxon", "V2", "Level", "V4", "V5"),
                           stringsAsFactors = FALSE)[, c("Taxon", "Level")]

  # Filter by kingdom
  kingdom <- sapply(taxa, `[`, 1)
  taxl.ba <- taxl[kingdom %in% c("Bacteria", "Archaea")]
  taxa.ba <- taxa[names(taxl.ba)]

  # Create taxonomy matrix for Bacteria and Archaea
  taxa.ba.mat <- matrix(sapply(taxa.ba, function(flds) {
    c(flds[1], flds[2], flds[3], flds[4], flds[5], flds[6])
  }), ncol = 6, byrow = TRUE)
  rownames(taxa.ba.mat) <- names(taxl.ba)

  # Build taxonomic strings for validation
  taxa.ba.mat.string <- matrix("UNDEF", nrow = nrow(taxa.ba.mat), ncol = ncol(taxa.ba.mat))
  rownames(taxa.ba.mat.string) <- names(taxl.ba)
  taxa.ba.mat.string[, 1] <- paste0(taxa.ba.mat[, 1], ";")
  for (col in seq(2, 6)) {
    taxa.ba.mat.string[, col] <- paste0(taxa.ba.mat.string[, col - 1], taxa.ba.mat[, col], ";")
  }
  if (any(taxa.ba.mat.string == "UNDEF")) stop("Taxon string matrix was not fully initialized.")

  # Validate taxonomy against reference
  taxa.ba.mat.is_valid <- matrix(taxa.ba.mat.string %in% silva.taxa$Taxon, ncol = 6)
  taxa.ba.mat[!taxa.ba.mat.is_valid] <- NA
  taxa.ba.mat[taxa.ba.mat %in% c("Uncultured", "uncultured")] <- NA

  # Handle species information if requested
  if (include.species) {
    taxa.ba.mat <- cbind(taxa.ba.mat, matrix(sapply(taxa.ba, `[`, 7), ncol = 1, byrow = TRUE))
    genus <- taxa.ba.mat[, 6]
    genus <- gsub("Candidatus ", "", genus)
    genus <- gsub("\\[|\\]", "", genus)

    binom <- taxa.ba.mat[, 7]
    binom <- gsub("Candidatus ", "", binom)
    binom <- gsub("\\[|\\]", "", binom)
    binom <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                   sapply(strsplit(binom, "\\s"), `[`, 2))

    gen.match <- mapply(dada2:::matchGenera, genus, binom[, 1], split.glyph = "-")
    is.NA <- apply(binom, 1, function(x) any(is.na(x)))
    is.sp <- grepl("sp\\.", binom[, 2])
    is.endo <- binom[, 1] %in% "endosymbiont" | binom[, 2] %in% "endosymbiont"
    is.uncult <- grepl("[Uu]ncultured", binom[, 1]) | grepl("[Uu]ncultured", binom[, 2])
    is.unident <- grepl("[Uu]nidentified", binom[, 1]) | grepl("[Uu]nidentified", binom[, 2])

    valid.spec <- gen.match & !is.NA & !is.sp & !is.endo & !is.uncult & !is.unident
    taxa.ba.mat[, 7] <- ifelse(valid.spec,
                               paste(genus, binom[, 2], sep = "_"),
                               NA)
  }

  # Sample Eukaryota sequences
  set.seed(100)
  euk_ids <- ref_ids[kingdom %in% "Eukaryota"]
  if (length(euk_ids) < n_euk) stop("Requested number of eukaryotic sequences exceeds available.")
  euk.keep <- sample(euk_ids, n_euk)

  # Create taxonomy matrix for Eukaryota
  taxa.euk.mat <- matrix("", nrow = n_euk, ncol = ncol(taxa.ba.mat))
  rownames(taxa.euk.mat) <- euk.keep
  taxa.euk.mat[, 1] <- "Eukaryota"
  taxa.euk.mat[, 2:ncol(taxa.euk.mat)] <- NA

  # Combine taxonomy matrices
  taxa.mat.final <- rbind(taxa.ba.mat, taxa.euk.mat)

  # Generate final taxonomy strings
  taxa.string.final <- apply(taxa.mat.final, 1, function(x) {
    tst <- paste(x, collapse = ";")
    tst <- paste0(tst, ";")
    gsub("NA;", "", tst)
  })

  # Validation checks
  if (any(is.na(names(taxa.string.final)))) stop("NA names in final taxon strings.")
  if (!all(names(taxa.string.final) %in% ref_ids)) stop("Some names don't match sequence IDs.")

  # Subset sequences for output
  xset.out <- xset[names(taxa.string.final)]

  # Prepare FASTA output with proper headers
  # First, create properly formatted FASTA headers with sequence IDs
  proper_headers <- names(xset.out)  # These are the sequence IDs

  # Create ShortRead object with proper sequence IDs
  short_read_out <- ShortRead(sread = xset.out, id = BStringSet(proper_headers))

  # Write output files
  writeFasta(short_read_out, fout_fasta, width = 20000L, compress = compress)

  taxonomy_df <- data.frame(ID = names(taxa.string.final), Taxonomy = taxa.string.final)
  write.table(taxonomy_df, fout_taxonomy, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  # Report summary statistics
  cat(length(xset.out), "reference sequences were output.\n")
  print(table(taxa.mat.final[, 1], useNA = "ifany"))
  if (include.species)
    cat(sum(!is.na(taxa.mat.final[, 7])), "entries include species names.\n")
}
