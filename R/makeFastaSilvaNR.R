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
#'                  include.species = TRUE, compress = TRUE, n_euk = 100)
#'}

# Check and install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Define the Bioconductor packages you need
bioc_packages <- c("Biostrings", "dada2")

# Check and install each Bioconductor package
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg)
  }
  # Load the package after ensuring it's installed
  library(pkg, character.only = TRUE)
}

library(Biostrings)
library(dada2)

makeFastaSilvaNR <- function(fin, ftax, fout_fasta, fout_taxonomy,
                             include.species = TRUE,
                             compress = FALSE,
                             n_euk = 100) {

  # 1. Read RNA sequences (Silva) from FASTA and convert to DNA
  xset <- DNAStringSet(readRNAStringSet(fin, format = "fasta"))

  # 2. Extract sequence descriptions and IDs
  descriptions <- names(xset)
  ref_ids <- sapply(strsplit(descriptions, "\\s"), `[`, 1)

  if (any(duplicated(ref_ids))) stop("FATAL: Duplicated sequence IDs detected in input FASTA.")

  names(xset) <- ref_ids # FASTA sequences will be named by ref_id internally

  # 3. Extract taxonomic annotation strings from descriptions
  # taxl_full_descriptions was names(xset) from input, used for parsing
  taxl_parsed <- gsub("^[A-Za-z0-9.]+\\s", "", descriptions) # Remove ID from original description string
  taxl_parsed <- gsub(";YM;", ";", taxl_parsed)
  taxa_list_from_fasta <- strsplit(taxl_parsed, ";")
  names(taxa_list_from_fasta) <- ref_ids

  # 4. Read taxonomy reference data (SILVA taxmap file)
  silva.taxa <- read.table(ftax, sep = "\t",
                           col.names = c("Taxon", "V2", "Level", "V4", "V5"),
                           stringsAsFactors = FALSE)[, c("Taxon", "Level")]

  # 5. Kingdom-level filter
  kingdom_from_fasta <- sapply(taxa_list_from_fasta, function(t) if(length(t) > 0) t[1] else NA_character_)

  is_ba_from_fasta <- kingdom_from_fasta %in% c("Bacteria", "Archaea")
  is_ba_from_fasta[is.na(is_ba_from_fasta)] <- FALSE

  ref_ids.ba <- names(taxa_list_from_fasta[is_ba_from_fasta])

  # 6. Construct taxonomy matrix for Bacteria/Archaea (up to Genus)
  num_ranks_ba <- 6 # K, P, C, O, F, G
  ba_colnames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  taxa.ba.mat <- matrix(NA_character_, nrow = length(ref_ids.ba), ncol = num_ranks_ba)
  colnames(taxa.ba.mat) <- ba_colnames

  if (length(ref_ids.ba) > 0) {
    rownames(taxa.ba.mat) <- ref_ids.ba
    current_taxa_ba_list <- taxa_list_from_fasta[ref_ids.ba]
    for (i in seq_along(current_taxa_ba_list)) {
      len <- length(current_taxa_ba_list[[i]])
      # Fill ranks, NAs will remain if original taxa string is shorter
      taxa.ba.mat[i, 1:min(num_ranks_ba, len)] <- current_taxa_ba_list[[i]][1:min(num_ranks_ba, len)]
    }
  }

  # 7. Validate Bacteria/Archaea taxonomy against SILVA reference (K to G)
  if (nrow(taxa.ba.mat) > 0) {
    taxa.ba.mat.string <- matrix(NA_character_, nrow = nrow(taxa.ba.mat), ncol = ncol(taxa.ba.mat))
    for (r in 1:nrow(taxa.ba.mat)) {
      current_path <- ""
      for (c in 1:ncol(taxa.ba.mat)) {
        if (!is.na(taxa.ba.mat[r,c]) && taxa.ba.mat[r,c] != "") {
          current_path <- paste0(current_path, taxa.ba.mat[r,c], ";")
          taxa.ba.mat.string[r,c] <- current_path
        } else {
          taxa.ba.mat[r,c] <- NA_character_
          break
        }
      }
    }

    taxa.ba.mat.is_valid <- matrix(TRUE, nrow = nrow(taxa.ba.mat.string), ncol = ncol(taxa.ba.mat.string))
    for (r_idx in 1:nrow(taxa.ba.mat.string)) {
      for (c_idx in 1:ncol(taxa.ba.mat.string)) {
        current_path_val <- taxa.ba.mat.string[r_idx, c_idx]
        if (!is.na(current_path_val)) {
          if (!(current_path_val %in% silva.taxa$Taxon)) {
            taxa.ba.mat.is_valid[r_idx, c_idx] <- FALSE
          }
        }
      }
    }

    for (r_idx in 1:nrow(taxa.ba.mat)) {
      for (c_idx in 1:ncol(taxa.ba.mat)) {
        if (!taxa.ba.mat.is_valid[r_idx, c_idx]) {
          taxa.ba.mat[r_idx, c_idx:ncol(taxa.ba.mat)] <- NA_character_
          break
        }
      }
    }

    taxa.ba.mat[grepl("uncultured|unknown|unidentified", taxa.ba.mat, ignore.case = TRUE)] <- NA_character_

    # Handle "Incertae Sedis" - remove from terminal positions and propagate backwards
    for (r_idx in 1:nrow(taxa.ba.mat)) {
      # Find the rightmost (terminal) non-NA position
      rightmost_non_na <- max(which(!is.na(taxa.ba.mat[r_idx, ])), 0)

      # Work backwards from the terminal position, removing "Incertae Sedis"
      if (rightmost_non_na > 0) {
        for (c_idx in rightmost_non_na:1) {
          if (!is.na(taxa.ba.mat[r_idx, c_idx]) &&
              grepl("Incertae Sedis", taxa.ba.mat[r_idx, c_idx], ignore.case = FALSE)) {
            taxa.ba.mat[r_idx, c_idx] <- NA_character_
          } else {
            # Stop when we hit a non-"Incertae Sedis" term
            break
          }
        }
      }
    }

    # Propagate NAs hierarchically after all modifications
    for (r_idx in 1:nrow(taxa.ba.mat)) {
      for (c_idx in 1:(ncol(taxa.ba.mat)-1)) { # Iterate up to second to last column
        if(is.na(taxa.ba.mat[r_idx, c_idx])) {
          # If a rank is NA, all subsequent ranks must also be NA
          taxa.ba.mat[r_idx, (c_idx+1):ncol(taxa.ba.mat)] <- NA_character_
        }
      }
    }
  }

  # 8. Handle species information if requested
  species_col_name <- "Species"
  final_colnames <- ba_colnames
  if (include.species) {
    final_colnames <- c(ba_colnames, species_col_name)

    # Initialize species column with NAs for all BA seqs
    species_values_col <- matrix(NA_character_, nrow = nrow(taxa.ba.mat), ncol = 1)
    colnames(species_values_col) <- species_col_name
    if(nrow(taxa.ba.mat) > 0) rownames(species_values_col) <- rownames(taxa.ba.mat)
    taxa.ba.mat <- cbind(taxa.ba.mat, species_values_col)

    if (nrow(taxa.ba.mat) > 0) {
      # Get the 7th field from the original Silva taxonomic annotation for B/A sequences
      raw_species_field <- sapply(taxa_list_from_fasta[rownames(taxa.ba.mat)],
                                  function(x) if(length(x) >= 7) x[7] else NA_character_)

      # Get validated genus from the matrix (6th column)
      genus_from_col6 <- taxa.ba.mat[, "Genus"] # This is already cleaned and validated to some extent
      genus_clean_for_match <- gsub("Candidatus ", "", genus_from_col6)
      genus_clean_for_match <- gsub("\\[|\\]", "", genus_clean_for_match)

      # Clean the raw 7th field (potential binomial string)
      binom_raw_clean <- gsub("Candidatus ", "", raw_species_field)
      binom_raw_clean <- gsub("\\[|\\]", "", binom_raw_clean)

      # Split the cleaned 7th field into two parts: potential genus and potential epithet
      # binom_parsed_parts will be a 2-column matrix
      binom_parsed_parts_list <- lapply(strsplit(binom_raw_clean, "\\s"), function(parts) {
        c( # First element is the first word (potential genus from species string)
          if(length(parts) >= 1) parts[1] else NA_character_,
          # Second element is the second word (potential specific epithet)
          if(length(parts) >= 2) parts[2] else NA_character_
        )
      })
      binom_parsed_parts <- do.call(rbind, binom_parsed_parts_list)
      colnames(binom_parsed_parts) <- c("parsed_genus", "parsed_epithet")

      # --- Validation Logic from Original Script ---
      gen.match <- rep(FALSE, nrow(taxa.ba.mat))
      # Only attempt matchGenera where both curated genus and parsed genus are available
      can_match_idx <- !is.na(genus_clean_for_match) & genus_clean_for_match != "" &
        !is.na(binom_parsed_parts[, "parsed_genus"]) & binom_parsed_parts[, "parsed_genus"] != ""

      if(any(can_match_idx)){
        gen.match[can_match_idx] <- mapply(dada2:::matchGenera,
                                           genus_clean_for_match[can_match_idx],
                                           binom_parsed_parts[can_match_idx, "parsed_genus"],
                                           MoreArgs = list(split.glyph = "-")) # Assuming default split.glyph
      }

      is.NA_epithet <- is.na(binom_parsed_parts[, "parsed_epithet"]) | binom_parsed_parts[, "parsed_epithet"] == ""
      # Use original script's checks on binom_parsed_parts (parsed_genus, parsed_epithet)
      is.sp <- grepl("sp\\.", binom_parsed_parts[, "parsed_epithet"])
      is.sp[is.na(is.sp)] <- FALSE # NA epithets are not "sp."

      # Check both parsed_genus and parsed_epithet for endo/uncult/unident
      is.endo <- (grepl("endosymbiont", binom_parsed_parts[, "parsed_genus"], ignore.case = TRUE) |
                    grepl("endosymbiont", binom_parsed_parts[, "parsed_epithet"], ignore.case = TRUE))
      is.endo[is.na(is.endo)] <- FALSE

      is.uncult <- (grepl("[Uu]ncultured", binom_parsed_parts[, "parsed_genus"]) |
                      grepl("[Uu]ncultured", binom_parsed_parts[, "parsed_epithet"]))
      is.uncult[is.na(is.uncult)] <- FALSE

      is.unident <- (grepl("[Uu]nidentified", binom_parsed_parts[, "parsed_genus"]) |
                       grepl("[Uu]nidentified", binom_parsed_parts[, "parsed_epithet"]))
      is.unident[is.na(is.unident)] <- FALSE

      valid.spec <- gen.match & !is.NA_epithet & !is.sp & !is.endo & !is.uncult & !is.unident
      valid.spec[is.na(valid.spec)] <- FALSE # Ensure no NAs

      # Store the full species name (Genus_epithet) if valid
      full_species_to_store <- character(length(valid.spec))
      full_species_to_store[] <- NA_character_ # Initialize with NAs

      # Conditions for forming the full species name:
      # 1. Species is valid according to the original script's logic (`valid.spec`)
      # 2. Both genus and epithet are not NA or empty
      can_form_full_species <- valid.spec &
        !is.na(genus_from_col6) & genus_from_col6 != "" &
        !is.na(binom_parsed_parts[, "parsed_epithet"]) & binom_parsed_parts[, "parsed_epithet"] != ""

      # Construct the full species name for valid entries
      if(any(can_form_full_species)) {
        full_species_to_store[can_form_full_species] <- paste0(
          genus_from_col6[can_form_full_species], "_",
          binom_parsed_parts[can_form_full_species, "parsed_epithet"]
        )
      }

      taxa.ba.mat[, species_col_name] <- full_species_to_store
    }
  }

  # 9. Sample Eukaryota sequences
  euk_indices <- which(kingdom_from_fasta == "Eukaryota" & !is.na(kingdom_from_fasta))
  euk_ids_available <- ref_ids[euk_indices]
  n_euk_actual <- 0
  euk.keep_ids <- character(0)

  if (n_euk > 0 && length(euk_ids_available) > 0) {
    if (length(euk_ids_available) < n_euk) {
      warning(paste("Requested", n_euk, "eukaryotic sequences, but only", length(euk_ids_available), "are available. Using all available."))
      n_euk_actual = length(euk_ids_available)
    } else {
      n_euk_actual = n_euk
    }
    set.seed(100)
    euk.keep_ids <- sample(euk_ids_available, n_euk_actual)
  }

  taxa.euk.mat <- matrix(NA_character_, nrow = n_euk_actual, ncol = length(final_colnames))
  colnames(taxa.euk.mat) <- final_colnames
  if (n_euk_actual > 0) {
    rownames(taxa.euk.mat) <- euk.keep_ids
    taxa.euk.mat[, "Kingdom"] <- "Eukaryota"
    # Other Euk ranks will be NA, correctly handled by prefixing function
  }

  # 10. Combine taxonomy matrices
  taxa.mat.final <- rbind(taxa.ba.mat, taxa.euk.mat)

  # 11. Combine sequences to keep & Perform final validation checks
  seqs.keep.ids <- rownames(taxa.mat.final) # These are the ref_ids for selected sequences
  if(length(seqs.keep.ids) > 0) { # Only perform checks if there are sequences to output
    if (any(is.na(seqs.keep.ids))) stop("FATAL: NA found in final sequence IDs before output.")
    if (!all(seqs.keep.ids %in% ref_ids)) stop("FATAL: Some final sequence IDs do not originate from input ref_ids.")
  }

  seqs.keep <- xset[seqs.keep.ids]
  # Ensure names of seqs.keep are just the IDs (already true by xset[seqs.keep.ids] as names(xset) are ref_ids)
  if(length(seqs.keep) != length(seqs.keep.ids)) stop("FATAL: Mismatch in number of sequences to keep and their IDs.")


  # 12. Write FASTA output
  if (length(seqs.keep) > 0) {
    if (compress) {
      writeXStringSet(seqs.keep, fout_fasta, format = "fasta", width = 20000L, compress = TRUE)
    } else {
      writeXStringSet(seqs.keep, fout_fasta, format = "fasta", width = 20000L)
    }
  } else {
    cat(">", file=fout_fasta) # Create empty but valid-ish fasta file
    warning("No sequences selected for FASTA output. Empty FASTA file created.")
  }


  # 13. Prepare and write Taxonomy TSV file
  if (nrow(taxa.mat.final) > 0) {
    prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

    formatted_taxonomy_strings <- apply(taxa.mat.final, 1, function(tax_row) {
      num_levels <- length(tax_row) # Will be 6 or 7
      current_prefixes <- prefixes[1:num_levels]
      parts <- character(num_levels)

      first_na_or_empty_hit = FALSE
      for (i in 1:num_levels) {
        if (first_na_or_empty_hit) {
          parts[i] <- current_prefixes[i]
        } else if (!is.na(tax_row[i]) && tax_row[i] != "") {
          parts[i] <- paste0(current_prefixes[i], tax_row[i])
        } else {
          parts[i] <- current_prefixes[i]
          first_na_or_empty_hit = TRUE
        }
      }
      paste(parts, collapse = ";")
    })

    taxonomy_output_df <- data.frame(
      ReferenceID = rownames(taxa.mat.final),
      Taxonomy = formatted_taxonomy_strings,
      stringsAsFactors = FALSE
    )

    write.table(taxonomy_output_df, file = fout_taxonomy,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    cat("ReferenceID\tTaxonomy\n", file = fout_taxonomy)
    warning("No sequences selected for taxonomy output. Taxonomy file with only headers created.")
  }

  # 14. Report summary statistics
  cat("\nProcessing Summary:\n")
  cat("--------------------\n")
  if(nrow(taxa.mat.final) > 0){
    cat(nrow(taxa.mat.final), "reference sequences were selected for output.\n")
    cat("Breakdown by Kingdom (from final processed taxonomy):\n")
    # Ensure Kingdom column exists before trying to table it
    if("Kingdom" %in% colnames(taxa.mat.final)){
      print(table(taxa.mat.final[, "Kingdom"], useNA = "ifany"))
    } else {
      cat("Kingdom column not found in final taxonomy matrix for breakdown.\n")
    }
    if (include.species && "Species" %in% colnames(taxa.mat.final)) {
      # Species column now contains only epithet or NA
      species_count <- sum(!is.na(taxa.mat.final[, "Species"]) & taxa.mat.final[, "Species"] != "")
      cat(species_count, "entries include species identity\n")
    }
  } else {
    cat("No sequences were selected for output.\n")
  }
  cat("--------------------\n\n")

  return(invisible(TRUE))
}

