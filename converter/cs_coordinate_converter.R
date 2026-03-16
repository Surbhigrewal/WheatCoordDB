# =============================================================
# cs_coordinate_converter.R
# =============================================================
# Auto-generated coordinate converter for CS v2.1 -> target assemblies
#
# Usage:
#   source("cs_coordinate_converter.R")
#   conv <- load_cs_converter("/path/to/conversion_tables")
#
#   # Single position
#   result <- cs_convert(conv, "chr1A", 500e6, "Mace")
#
#   # Region
#   result <- cs_convert_region(conv, "chr1A", 490e6, 510e6, "Mace")
#
#   # BED data frame
#   result_df <- cs_convert_bed(conv, my_bed_df, "SY_Mattis")
#   
#   # List targets
#   cs_list_targets(conv)
# =============================================================

library(data.table)

#' Load coordinate converter from tables directory
#' @param tables_dir Path to conversion_tables/ directory
#' @return Converter object (list with tables and functions)
load_cs_converter <- function(tables_dir) {
  if (!dir.exists(tables_dir)) {
    stop(paste("Conversion tables directory not found:", tables_dir))
  }
  
  conv <- list(
    tables_dir = tables_dir,
    cache = list()
  )
  class(conv) <- "CSConverter"
  return(conv)
}

#' List available target assemblies
cs_list_targets <- function(conv) {
  dirs <- list.dirs(conv$tables_dir, full.names = FALSE, recursive = FALSE)
  return(sort(dirs))
}

#' List available chromosomes for a target
cs_list_chromosomes <- function(conv, target_name) {
  tgt_dir <- file.path(conv$tables_dir, target_name)
  files <- list.files(tgt_dir, pattern = "_conversion\\.tsv\\.gz$")
  gsub("_conversion\\.tsv\\.gz", "", files)
}

# Internal: load and cache a conversion table
.load_table <- function(conv, target_name, cs_chr) {
  key <- paste(target_name, cs_chr, sep = "_")
  
  if (!is.null(conv$cache[[key]])) {
    return(conv$cache[[key]])
  }
  
  table_file <- file.path(conv$tables_dir, target_name,
                          paste0(cs_chr, "_conversion.tsv.gz"))
  
  if (!file.exists(table_file)) {
    # Try without chr prefix
    alt_chr <- gsub("^chr", "", cs_chr)
    table_file <- file.path(conv$tables_dir, target_name,
                            paste0(alt_chr, "_conversion.tsv.gz"))
    if (!file.exists(table_file)) {
      stop(paste("No conversion table for", target_name, cs_chr))
    }
  }
  
  dt <- fread(table_file)
  conv$cache[[key]] <<- dt
  return(dt)
}

#' Convert a single CS position to target coordinates
#' @param conv Converter object from load_cs_converter()
#' @param cs_chr CS chromosome name (e.g. "chr1A")  
#' @param cs_pos CS position (integer)
#' @param target_name Target assembly name
#' @return Named list: tgt_chr, tgt_pos, anchor_density
cs_convert <- function(conv, cs_chr, cs_pos, target_name) {
  dt <- .load_table(conv, target_name, cs_chr)
  
  tgt_chr <- dt$tgt_chr[1]
  tgt_pos <- as.integer(approx(dt$cs_pos, dt$tgt_pos, xout = cs_pos, 
                                rule = 2)$y)
  density <- as.integer(approx(dt$cs_pos, dt$anchor_density_5Mb, xout = cs_pos,
                                rule = 2)$y)
  
  return(list(tgt_chr = tgt_chr, tgt_pos = tgt_pos, anchor_density = density))
}

#' Convert a CS region to target coordinates
#' @return Named list: tgt_chr, tgt_start, tgt_end, mean_anchor_density
cs_convert_region <- function(conv, cs_chr, cs_start, cs_end, target_name) {
  dt <- .load_table(conv, target_name, cs_chr)
  
  tgt_chr <- dt$tgt_chr[1]
  tgt_start <- as.integer(approx(dt$cs_pos, dt$tgt_pos, xout = cs_start, rule = 2)$y)
  tgt_end   <- as.integer(approx(dt$cs_pos, dt$tgt_pos, xout = cs_end,   rule = 2)$y)
  
  sample_pts <- seq(cs_start, cs_end, length.out = 20)
  mean_density <- mean(approx(dt$cs_pos, dt$anchor_density_5Mb, xout = sample_pts, rule = 2)$y)
  
  return(list(
    tgt_chr = tgt_chr,
    tgt_start = min(tgt_start, tgt_end),
    tgt_end   = max(tgt_start, tgt_end),
    mean_anchor_density = round(mean_density, 1)
  ))
}

#' Convert a data frame of BED regions
#' @param conv Converter
#' @param bed_df Data frame with columns: chr, start, end (and optional extras)
#' @param target_name Target assembly name
#' @return Data frame with added tgt_chr, tgt_start, tgt_end columns
cs_convert_bed <- function(conv, bed_df, target_name) {
  colnames(bed_df)[1:3] <- c("cs_chr", "cs_start", "cs_end")
  
  results <- lapply(seq_len(nrow(bed_df)), function(i) {
    tryCatch({
      r <- cs_convert_region(conv,
                             bed_df$cs_chr[i],
                             bed_df$cs_start[i], 
                             bed_df$cs_end[i],
                             target_name)
      data.frame(
        cs_chr   = bed_df$cs_chr[i],
        cs_start = bed_df$cs_start[i],
        cs_end   = bed_df$cs_end[i],
        tgt_chr  = r$tgt_chr,
        tgt_start = r$tgt_start,
        tgt_end   = r$tgt_end,
        anchor_density = r$mean_anchor_density
      )
    }, error = function(e) {
      warning(paste("Failed to convert row", i, ":", conditionMessage(e)))
      NULL
    })
  })
  
  result_df <- do.call(rbind, Filter(Negate(is.null), results))
  
  # Append extra BED columns if present
  if (ncol(bed_df) > 3) {
    extra <- bed_df[, 4:ncol(bed_df), drop = FALSE]
    result_df <- cbind(result_df, extra[!sapply(results, is.null), , drop = FALSE])
  }
  
  return(result_df)
}

#' Print summary of available targets and chromosomes
print.CSConverter <- function(conv, ...) {
  targets <- cs_list_targets(conv)
  cat("CS Coordinate Converter\n")
  cat("Tables dir:", conv$tables_dir, "\n")
  cat("Available targets:", paste(targets, collapse = ", "), "\n")
}
