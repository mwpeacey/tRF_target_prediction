################################################################################
# shiny_target_browser.R
#
# Query-backed Shiny browser for annotated tRF target sites. Designed for large
# datasets by pushing filtering, sorting, and pagination into DuckDB, then
# rendering only the selected alignment card in Shiny.
#
# Run from an interactive R session:
#   source("R_scripts/shiny_target_browser.R")
#   shiny::runApp(app)
#
# Or from the shell:
#   Rscript R_scripts/run_shiny_target_browser.R
################################################################################

required_packages <- c("shiny", "DBI", "duckdb", "DT")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them with install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

data_path <- Sys.getenv(
  "TRF_TARGET_BROWSER_CSV",
  unset = "import/miranda/miranda_output_unique_annotated.csv"
)
db_path <- Sys.getenv(
  "TRF_TARGET_BROWSER_DB",
  unset = "import/miranda/target_browser.duckdb"
)
default_page_size <- as.integer(Sys.getenv("TRF_TARGET_BROWSER_PAGE_SIZE", unset = "100"))
csv_available <- nzchar(data_path) && file.exists(data_path)
db_available <- nzchar(db_path) && file.exists(db_path)
read_only_default <- if (csv_available) "false" else "true"
read_only_db <- identical(
  tolower(Sys.getenv("TRF_TARGET_BROWSER_READ_ONLY", unset = read_only_default)),
  "true"
)

if (!csv_available && !db_available) {
  stop(
    "Could not find either the target-site CSV or the DuckDB database.\n",
    "Checked CSV: ", normalizePath(data_path, mustWork = FALSE), "\n",
    "Checked DuckDB: ", normalizePath(db_path, mustWork = FALSE)
  )
}

if (!read_only_db) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
}

normalize_match_display <- function(query_aln, match_str, ref_aln) {
  strip_label <- function(x, label) {
    sub(paste0("^\\s*", label, ":\\s*"), "", x)
  }

  extract_core <- function(x) {
    x <- sub("^\\s*[35]'\\s*", "", x)
    sub("\\s*[35]'\\s*$", "", x)
  }

  edge_lowercase_width <- function(x, side = c("left", "right")) {
    side <- match.arg(side)
    pattern <- if (side == "left") "^[a-z]+" else "[a-z]+$"
    match <- regmatches(x, regexpr(pattern, x))

    if (length(match) == 0L || identical(match, character(0)) || identical(match, "")) {
      return(0L)
    }

    nchar(match, type = "chars")
  }

  leading_space_width <- function(x) {
    match <- regmatches(x, regexpr("^\\s*", x))

    if (length(match) == 0L || identical(match, character(0)) || is.na(match)) {
      return(0L)
    }

    nchar(match, type = "chars")
  }

  trailing_space_width <- function(x) {
    match <- regmatches(x, regexpr("\\s*$", x))

    if (length(match) == 0L || identical(match, character(0)) || is.na(match)) {
      return(0L)
    }

    nchar(match, type = "chars")
  }

  normalize_match_string <- function(
    match_str,
    display_width,
    label_prefix_width,
    desired_left_pad,
    desired_right_pad
  ) {
    if (is.na(match_str) || is.null(match_str)) {
      return("")
    }

    lead_spaces <- leading_space_width(match_str)
    if (label_prefix_width > 0L && lead_spaces >= label_prefix_width) {
      match_str <- substring(match_str, label_prefix_width + 1L)
    }

    match_width <- nchar(match_str, type = "chars")
    if (match_width == display_width) {
      return(match_str)
    }

    lead_spaces <- leading_space_width(match_str)
    trail_spaces <- trailing_space_width(match_str)

    if (match_width > display_width) {
      extra_width <- match_width - display_width
      trim_left <- min(extra_width, max(0L, lead_spaces - desired_left_pad))
      extra_width <- extra_width - trim_left

      trim_right <- min(extra_width, max(0L, trail_spaces - desired_right_pad))
      extra_width <- extra_width - trim_right

      if (trim_left > 0L) {
        match_str <- substring(match_str, trim_left + 1L)
      }

      if (trim_right > 0L) {
        match_str <- substring(
          match_str,
          1L,
          nchar(match_str, type = "chars") - trim_right
        )
      }

      return(match_str)
    }

    missing_width <- display_width - match_width
    left_missing <- max(0L, desired_left_pad - lead_spaces)
    right_missing <- max(0L, desired_right_pad - trail_spaces)
    left_pad <- min(left_missing, missing_width)
    right_pad <- min(right_missing, missing_width - left_pad)
    remaining_pad <- missing_width - left_pad - right_pad

    paste0(
      strrep(" ", left_pad + remaining_pad),
      match_str,
      strrep(" ", right_pad)
    )
  }

  base_bounds <- function(x) {
    base_pos <- gregexpr("[A-Za-z-]", x)[[1]]
    if (length(base_pos) == 1L && base_pos[1] == -1L) {
      return(list(left = 0L, right = 0L))
    }

    first_base <- base_pos[1]
    last_base <- base_pos[length(base_pos)]

    list(
      left = first_base - 1L,
      right = nchar(x, type = "chars") - last_base
    )
  }

  q_seq <- strip_label(query_aln, "Query")
  r_seq <- strip_label(ref_aln, "Ref")
  q_core <- extract_core(q_seq)
  r_core <- extract_core(r_seq)

  q_bounds <- base_bounds(q_seq)
  r_bounds <- base_bounds(r_seq)

  left_pad <- max(q_bounds$left, r_bounds$left)
  right_pad <- max(q_bounds$right, r_bounds$right)
  overhang_left <- max(
    edge_lowercase_width(q_core, "left"),
    edge_lowercase_width(r_core, "left")
  )
  overhang_right <- max(
    edge_lowercase_width(q_core, "right"),
    edge_lowercase_width(r_core, "right")
  )
  display_width <- max(
    nchar(q_seq, type = "chars"),
    nchar(r_seq, type = "chars")
  )
  label_prefix_width <- max(
    nchar(query_aln, type = "chars") - nchar(q_seq, type = "chars"),
    nchar(ref_aln, type = "chars") - nchar(r_seq, type = "chars")
  )
  match_str <- normalize_match_string(
    match_str = match_str,
    display_width = display_width,
    label_prefix_width = label_prefix_width,
    desired_left_pad = left_pad + overhang_left,
    desired_right_pad = right_pad + overhang_right
  )

  list(
    query_label = "tRF:",
    query_seq = q_seq,
    match_seq = match_str,
    ref_label = "target:",
    ref_seq = r_seq
  )
}

escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

color_match_html <- function(match_str) {
  if (is.na(match_str) || is.null(match_str)) {
    return("N/A")
  }

  chars <- strsplit(match_str, "")[[1]]
  colored <- vapply(chars, function(ch) {
    if (ch == "|") return('<span class="wc">|</span>')
    if (ch == ":") return('<span class="gu">:</span>')
    if (ch == " ") return(" ")
    paste0('<span class="mm">', escape_html(ch), '</span>')
  }, character(1))

  paste(colored, collapse = "")
}

as_display_string <- function(x, fallback = "\u2014") {
  if (length(x) == 0L || is.null(x) || is.na(x) || identical(x, "")) {
    fallback
  } else {
    as.character(x)
  }
}

build_alignment_card <- function(row) {
  if (is.null(row) || nrow(row) == 0L) {
    return(shiny::div(class = "empty-state", "Select a row to view the alignment."))
  }

  row <- row[1, , drop = FALSE]
  alignment <- normalize_match_display(
    row$query_alignment[[1]],
    row$match_string[[1]],
    row$ref_alignment[[1]]
  )

  ltr_info <- if (!is.na(row$LTR_family[[1]]) && !identical(row$LTR_family[[1]], "")) {
    paste0(as_display_string(row$LTR_family[[1]]), " / ", as_display_string(row$LTR_gene_id[[1]]))
  } else {
    "None"
  }

  imprint_tag <- if (!is.na(row$imprint_status[[1]]) && !identical(row$imprint_status[[1]], "")) {
    shiny::HTML(sprintf(
      '<span class="imprinted-tag %s">Imprint: %s</span>',
      if (identical(row$imprint_status[[1]], "Imprinted")) "imprinted-yes" else "imprinted-no",
      escape_html(row$imprint_status[[1]])
    ))
  } else {
    NULL
  }

  loc_text <- sprintf(
    "%s:%s-%s (%s)",
    as_display_string(row$seqnames[[1]]),
    format(as.numeric(row$start[[1]]), big.mark = ",", trim = TRUE, scientific = FALSE),
    format(as.numeric(row$end[[1]]), big.mark = ",", trim = TRUE, scientific = FALSE),
    as_display_string(row$strand[[1]])
  )

  shiny::div(
    class = "alignment-card",
    shiny::div(
      class = "meta",
      shiny::tags$span(class = "trf", as_display_string(row$tRF[[1]])),
      shiny::tags$span(class = "anticodon", as_display_string(row$tRNA_anticodon[[1]], "Anticodon N/A")),
      shiny::tags$span(class = "gene", as_display_string(row$gencode_gene_name[[1]], "\u2014")),
      shiny::tags$span(class = "loc", loc_text),
      shiny::tags$span(
        class = "score",
        sprintf("Score: %s | Energy: %s",
                as_display_string(row$alignment_score[[1]]),
                as_display_string(row$energy[[1]]))
      )
    ),
    shiny::div(
      class = "meta meta-annotations",
      shiny::tags$span(
        class = paste("ltr-tag", if (identical(row$LTR[[1]], TRUE) || identical(row$LTR[[1]], "TRUE")) "ltr-yes" else "ltr-no"),
        paste0("LTR: ", ltr_info)
      ),
      shiny::tags$span(
        class = paste("gag-tag", if (identical(row$gag_gene[[1]], TRUE) || identical(row$gag_gene[[1]], "TRUE")) "gag-yes" else "gag-no"),
        paste0("Gag: ", if (identical(row$gag_gene[[1]], TRUE) || identical(row$gag_gene[[1]], "TRUE")) "Yes" else "No")
      ),
      shiny::tags$span(
        class = paste("gencode-tag", if (tolower(as_display_string(row$gencode_location[[1]])) == "intergenic") "gencode-intergenic" else "gencode-genic"),
        paste0("GENCODE: ", as_display_string(row$gencode_location[[1]]))
      ),
      shiny::tags$span(
        class = paste("stringtie-tag", if (tolower(as_display_string(row$stringtie_location[[1]])) == "intergenic") "stringtie-intergenic" else "stringtie-genic"),
        paste0("StringTie: ", as_display_string(row$stringtie_location[[1]]))
      ),
      imprint_tag
    ),
    shiny::div(
      class = "alignment-block",
      shiny::div(
        class = "alignment-grid",
        shiny::tags$span(class = "alignment-label alignment-label-query", alignment$query_label),
        shiny::tags$span(class = "alignment-seq", alignment$query_seq),
        shiny::tags$span(class = "alignment-label", ""),
        shiny::tags$span(class = "alignment-seq alignment-match", shiny::HTML(color_match_html(alignment$match_seq))),
        shiny::tags$span(class = "alignment-label alignment-label-target", alignment$ref_label),
        shiny::tags$span(class = "alignment-seq", alignment$ref_seq)
      )
    )
  )
}

sql_escape <- function(x) {
  gsub("'", "''", x, fixed = TRUE)
}

parse_coord_range <- function(coord_text) {
  if (is.null(coord_text)) {
    return(NULL)
  }

  text <- trimws(coord_text)
  if (!nzchar(text)) {
    return(NULL)
  }

  # Normalize: tabs -> spaces, ".." -> "-", strip commas, collapse whitespace.
  text <- gsub("\\t", " ", text)
  text <- gsub("\\.\\.", "-", text)
  text <- gsub(",", "", text)
  text <- gsub("\\s+", " ", text)

  # Accept several common formats:
  #   chr1:1000-2000
  #   chr1:1000 2000
  #   chr1 1000-2000
  #   chr1 1000 2000   (BED-like)
  patterns <- c(
    "^(\\S+):\\s*([0-9]+)\\s*-\\s*([0-9]+)$",
    "^(\\S+):\\s*([0-9]+)\\s+([0-9]+)$",
    "^(\\S+)\\s+([0-9]+)\\s*-\\s*([0-9]+)$",
    "^(\\S+)\\s+([0-9]+)\\s+([0-9]+)$"
  )

  m <- NULL
  for (pat in patterns) {
    candidate <- regmatches(text, regexec(pat, text))[[1]]
    if (length(candidate) == 4L) {
      m <- candidate
      break
    }
  }

  if (is.null(m)) {
    return(list(error = "Try chr1:1,000,000-2,000,000 (also accepts 'chr1 1000 2000')."))
  }

  chrom <- m[2]
  start_val <- suppressWarnings(as.numeric(m[3]))
  end_val <- suppressWarnings(as.numeric(m[4]))

  if (is.na(start_val) || is.na(end_val)) {
    return(list(error = "Could not parse start/end positions."))
  }

  if (start_val > end_val) {
    tmp <- start_val
    start_val <- end_val
    end_val <- tmp
  }

  # Accept either "chr1" or "1" — match both forms against the data.
  alt_chrom <- if (grepl("^chr", chrom, ignore.case = TRUE)) {
    sub("^chr", "", chrom, ignore.case = TRUE)
  } else {
    paste0("chr", chrom)
  }
  seqnames <- unique(c(chrom, alt_chrom))

  list(seqname = chrom, seqnames = seqnames, start = start_val, end = end_val)
}

format_coord_status <- function(parsed) {
  if (is.null(parsed)) {
    return("")
  }

  if (!is.null(parsed$error)) {
    return(parsed$error)
  }

  width <- parsed$end - parsed$start + 1
  sprintf(
    "Showing sites overlapping %s:%s-%s (%s bp window).",
    parsed$seqname,
    format(parsed$start, big.mark = ",", trim = TRUE, scientific = FALSE),
    format(parsed$end, big.mark = ",", trim = TRUE, scientific = FALSE),
    format(width, big.mark = ",", trim = TRUE, scientific = FALSE)
  )
}

target_table_sql <- function(csv_path) {
  paste0(
    "CREATE OR REPLACE TABLE targets AS ",
    "SELECT ",
    "  row_number() OVER () AS report_row_id, ",
    "  seqnames, ",
    "  CAST(start AS BIGINT) AS start, ",
    "  CAST(\"end\" AS BIGINT) AS \"end\", ",
    "  strand, ",
    "  tRF, ",
    "  NULLIF(tRNA_anticodon, 'NA') AS tRNA_anticodon, ",
    "  CAST(alignment_score AS DOUBLE) AS alignment_score, ",
    "  CAST(energy AS DOUBLE) AS energy, ",
    "  CASE ",
    "    WHEN UPPER(COALESCE(CAST(LTR AS VARCHAR), 'FALSE')) = 'TRUE' THEN TRUE ",
    "    ELSE FALSE ",
    "  END AS LTR, ",
    "  NULLIF(LTR_family, 'NA') AS LTR_family, ",
    "  NULLIF(LTR_gene_id, 'NA') AS LTR_gene_id, ",
    "  NULLIF(gencode_gene_name, 'NA') AS gencode_gene_name, ",
    "  NULLIF(gencode_location, 'NA') AS gencode_location, ",
    "  NULLIF(stringtie_location, 'NA') AS stringtie_location, ",
    "  NULLIF(imprint_status, 'NA') AS imprint_status, ",
    "  CASE ",
    "    WHEN UPPER(COALESCE(CAST(gag_gene AS VARCHAR), 'FALSE')) = 'TRUE' THEN TRUE ",
    "    ELSE FALSE ",
    "  END AS gag_gene, ",
    "  query_alignment, match_string, ref_alignment, ",
    "  LOWER(CONCAT_WS(' ', ",
    "    COALESCE(tRF, ''), ",
    "    COALESCE(NULLIF(tRNA_anticodon, 'NA'), ''), ",
    "    COALESCE(NULLIF(gencode_gene_name, 'NA'), ''), ",
    "    COALESCE(seqnames, ''), ",
    "    COALESCE(CAST(start AS VARCHAR), ''), ",
    "    COALESCE(CAST(\"end\" AS VARCHAR), ''), ",
    "    COALESCE(NULLIF(LTR_family, 'NA'), ''), ",
    "    COALESCE(NULLIF(LTR_gene_id, 'NA'), ''), ",
    "    COALESCE(NULLIF(gencode_location, 'NA'), ''), ",
    "    COALESCE(NULLIF(stringtie_location, 'NA'), ''), ",
    "    COALESCE(NULLIF(imprint_status, 'NA'), '') ",
    "  )) AS search_text ",
    "FROM read_csv_auto('", sql_escape(normalizePath(csv_path)), "', HEADER = TRUE, SAMPLE_SIZE = -1)"
  )
}

ensure_target_table <- function(con, csv_path, db_file, csv_available) {
  has_table <- DBI::dbExistsTable(con, "targets")

  if (!csv_available) {
    if (!has_table) {
      stop(
        "The DuckDB database does not contain a 'targets' table, and no CSV was available to build it.\n",
        "DuckDB checked: ", normalizePath(db_file, mustWork = FALSE)
      )
    }

    return(invisible(NULL))
  }

  needs_refresh <- !file.exists(db_file) || file.info(db_file)$mtime < file.info(csv_path)$mtime

  if (!has_table || needs_refresh) {
    DBI::dbExecute(con, "DROP TABLE IF EXISTS targets")
    DBI::dbExecute(con, target_table_sql(csv_path))
  }
}

build_where_clause <- function(
  con,
  search,
  ltr_filter,
  gag_filter,
  imprinted_filter,
  gencode_locations,
  stringtie_locations,
  score_min,
  score_max,
  coord_range = NULL
) {
  clauses <- c()

  if (!is.null(coord_range) && is.null(coord_range$error)) {
    seqname_candidates <- coord_range$seqnames
    if (is.null(seqname_candidates) || length(seqname_candidates) == 0L) {
      seqname_candidates <- coord_range$seqname
    }
    quoted_seqnames <- vapply(
      seqname_candidates,
      function(x) as.character(DBI::dbQuoteString(con, x)),
      character(1)
    )
    clauses <- c(
      clauses,
      sprintf(
        "seqnames IN (%s) AND start <= %s AND \"end\" >= %s",
        paste(quoted_seqnames, collapse = ", "),
        format(as.numeric(coord_range$end), scientific = FALSE),
        format(as.numeric(coord_range$start), scientific = FALSE)
      )
    )
  }

  if (!is.null(search) && nzchar(trimws(search))) {
    pattern <- paste0("%", trimws(search), "%")
    pattern_sql <- as.character(DBI::dbQuoteString(con, pattern))
    clauses <- c(clauses, paste0("search_text ILIKE ", pattern_sql))
  }

  if (identical(ltr_filter, "ltr")) {
    clauses <- c(clauses, "LTR = TRUE")
  } else if (identical(ltr_filter, "no-ltr")) {
    clauses <- c(clauses, "LTR = FALSE")
  }

  if (identical(gag_filter, "gag")) {
    clauses <- c(clauses, "gag_gene = TRUE")
  } else if (identical(gag_filter, "no-gag")) {
    clauses <- c(clauses, "gag_gene = FALSE")
  }

  if (identical(imprinted_filter, "imprinted")) {
    clauses <- c(clauses, "imprint_status = 'Imprinted'")
  } else if (identical(imprinted_filter, "non-imprinted")) {
    clauses <- c(clauses, "imprint_status = 'Not imprinted'")
  }

  if (!is.null(gencode_locations) && length(gencode_locations) > 0L) {
    selected_locations <- unique(trimws(as.character(gencode_locations)))
    selected_locations <- selected_locations[nzchar(selected_locations)]

    if (length(selected_locations) > 0L) {
      quoted_locations <- vapply(
        selected_locations,
        function(x) as.character(DBI::dbQuoteString(con, x)),
        character(1)
      )
      clauses <- c(
        clauses,
        paste0("COALESCE(gencode_location, '') IN (", paste(quoted_locations, collapse = ", "), ")")
      )
    }
  }

  if (!is.null(stringtie_locations) && length(stringtie_locations) > 0L) {
    selected_locations <- unique(trimws(as.character(stringtie_locations)))
    selected_locations <- selected_locations[nzchar(selected_locations)]

    if (length(selected_locations) > 0L) {
      quoted_locations <- vapply(
        selected_locations,
        function(x) as.character(DBI::dbQuoteString(con, x)),
        character(1)
      )
      clauses <- c(
        clauses,
        paste0("COALESCE(stringtie_location, '') IN (", paste(quoted_locations, collapse = ", "), ")")
      )
    }
  }

  if (!is.null(score_min) && is.finite(score_min)) {
    clauses <- c(clauses, sprintf("alignment_score >= %s", as.numeric(score_min)))
  }
  if (!is.null(score_max) && is.finite(score_max)) {
    clauses <- c(clauses, sprintf("alignment_score <= %s", as.numeric(score_max)))
  }

  if (length(clauses) == 0L) {
    ""
  } else {
    paste("WHERE", paste(clauses, collapse = " AND "))
  }
}

normalize_table_sort <- function(table_sort) {
  if (is.null(table_sort) || identical(table_sort$mode, "original")) {
    return(NULL)
  }

  if (is.null(table_sort$col) || is.null(table_sort$dir)) {
    return(NULL)
  }

  col <- as.integer(table_sort$col)
  dir <- if (identical(table_sort$dir, "desc")) "desc" else "asc"

  valid_cols <- c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L)
  if (!(col %in% valid_cols)) {
    return(NULL)
  }

  list(col = col, dir = dir)
}

sort_clause <- function(table_sort) {
  sort_spec <- normalize_table_sort(table_sort)
  if (is.null(sort_spec)) {
    return("ORDER BY report_row_id ASC")
  }

  sort_columns <- switch(
    as.character(sort_spec$col),
    "0" = c("report_row_id"),
    "1" = c("LOWER(COALESCE(tRF, ''))"),
    "2" = c("LOWER(COALESCE(tRNA_anticodon, ''))"),
    "3" = c("LOWER(COALESCE(seqnames, ''))", "start", "\"end\"", "LOWER(COALESCE(strand, ''))"),
    "4" = c("alignment_score"),
    "5" = c("energy"),
    "6" = c("LOWER(COALESCE(gencode_gene_name, ''))"),
    "7" = c("LOWER(COALESCE(gencode_location, ''))"),
    "8" = c("LOWER(COALESCE(stringtie_location, ''))"),
    "9" = c("LOWER(CASE WHEN LTR THEN COALESCE(LTR_family, 'LTR') ELSE 'None' END)"),
    "10" = c("LOWER(COALESCE(imprint_status, ''))"),
    c("report_row_id")
  )

  dir_sql <- if (identical(sort_spec$dir, "desc")) "DESC" else "ASC"
  order_terms <- sprintf("%s %s", sort_columns, dir_sql)

  if (!("report_row_id" %in% sort_columns)) {
    order_terms <- c(order_terms, "report_row_id ASC")
  }

  paste("ORDER BY", paste(order_terms, collapse = ", "))
}

sort_label <- function(table_sort) {
  sort_spec <- normalize_table_sort(table_sort)
  if (is.null(sort_spec)) {
    return("Current sort: original order")
  }

  column_labels <- c(
    "Row ID",
    "tRF",
    "Anticodon",
    "Location",
    "Alignment score",
    "Energy",
    "GENCODE gene",
    "GENCODE region",
    "StringTie region",
    "LTR",
    "Imprint status"
  )

  dir_label <- if (identical(sort_spec$dir, "desc")) "descending" else "ascending"
  sprintf("Current sort: %s (%s)", column_labels[sort_spec$col + 1L], dir_label)
}

con <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = db_path,
  read_only = read_only_db
)
ensure_target_table(con, data_path, db_path, csv_available = csv_available)

shiny::onStop(function() {
  try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
})

score_limits <- DBI::dbGetQuery(
  con,
  "SELECT FLOOR(MIN(alignment_score)) AS min_score, CEIL(MAX(alignment_score)) AS max_score FROM targets"
)
score_min_raw <- as.numeric(score_limits$min_score[[1]])
score_max_raw <- as.numeric(score_limits$max_score[[1]])
score_min <- floor(score_min_raw / 10) * 10
score_max <- ceiling(score_max_raw / 10) * 10
gencode_location_counts <- DBI::dbGetQuery(
  con,
  paste(
    "SELECT gencode_location, COUNT(*) AS n",
    "FROM targets",
    "WHERE gencode_location IS NOT NULL",
    "GROUP BY gencode_location",
    "ORDER BY n DESC, gencode_location ASC"
  )
)
gencode_location_choices <- setNames(
  as.list(gencode_location_counts$gencode_location),
  sprintf(
    "%s (%s)",
    gencode_location_counts$gencode_location,
    format(gencode_location_counts$n, big.mark = ",")
  )
)
stringtie_location_counts <- DBI::dbGetQuery(
  con,
  paste(
    "SELECT stringtie_location, COUNT(*) AS n",
    "FROM targets",
    "WHERE stringtie_location IS NOT NULL",
    "GROUP BY stringtie_location",
    "ORDER BY n DESC, stringtie_location ASC"
  )
)
stringtie_location_choices <- setNames(
  as.list(stringtie_location_counts$stringtie_location),
  sprintf(
    "%s (%s)",
    stringtie_location_counts$stringtie_location,
    format(stringtie_location_counts$n, big.mark = ",")
  )
)

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$title("3'-tRF Target Explorer"),
    shiny::tags$style(shiny::HTML("
      :root {
        --bg: #ffffff;
        --panel: #ffffff;
        --border: #dddddd;
        --ink: #222222;
        --muted: #666666;
        --color-trf: #b22425;
        --color-trf-light: #faf0f0;
        --color-target: #e1b547;
        --color-target-light: #fdf8ec;
        --color-gencode: #36648b;
        --color-gencode-light: #edf3f8;
        --accent: #b22425;
      }
      body {
        background: var(--bg);
        color: var(--ink);
        font-family: Arial, Helvetica, 'Helvetica Neue', sans-serif;
        font-size: 14px;
        line-height: 1.5;
      }

      /* ── Header ── */
      .hero {
        background: var(--bg);
        border: none;
        border-radius: 0;
        padding: 0;
        margin: 0 0 8px 0;
        box-shadow: none;
        border-bottom: 2px solid var(--border);
      }
      .hero-inner {
        padding: 24px 0 18px 0;
      }
      .hero-title {
        margin: 0 0 10px 0;
        color: var(--ink);
        font-family: Arial, Helvetica, sans-serif;
        font-size: 26px;
        line-height: 1.25;
        font-weight: 700;
      }
      .hero-subtitle {
        color: var(--muted);
        font-size: 13.5px;
        line-height: 1.5;
        margin: 0;
        max-width: 820px;
      }
      .hero-subtitle em {
        font-style: italic;
      }

      /* ── Panels ── */
      .control-panel, .results-panel {
        background: var(--panel);
        border: 1px solid var(--border);
        border-radius: 8px;
        padding: 18px;
        margin-bottom: 16px;
        box-shadow: none;
      }
      .control-panel {
        position: sticky;
        top: 16px;
      }
      .panel-title {
        font-size: 11.5px;
        font-weight: 700;
        letter-spacing: 0.06em;
        text-transform: uppercase;
        color: var(--ink);
        margin-bottom: 14px;
        padding-bottom: 7px;
        border-bottom: 1px solid #eeeeee;
      }
      .results-summary {
        color: var(--muted);
        margin-bottom: 10px;
        font-size: 13px;
      }

      /* ── Alignment card ── */
      .alignment-card {
        background: #ffffff;
        border: 1px solid var(--border);
        border-radius: 8px;
        padding: 16px;
        margin-bottom: 16px;
        box-shadow: none;
        position: sticky;
        top: 0;
        z-index: 10;
      }
      .meta {
        margin-bottom: 4px;
        line-height: 1.85;
      }
      .meta-annotations {
        margin-bottom: 10px;
      }
      .meta span {
        display: inline-block;
        margin-right: 16px;
        font-size: 13px;
        color: var(--muted);
      }
      .trf {
        font-weight: 700;
        color: var(--ink) !important;
      }
      .anticodon {
        font-weight: 600;
        color: var(--ink) !important;
      }
      .gene {
        font-weight: 700;
        color: var(--ink) !important;
        font-style: italic;
      }
      .loc {
        font-variant-numeric: tabular-nums;
      }
      .score {
        font-variant-numeric: tabular-nums;
      }
      .ltr-tag {
        padding: 1px 7px;
        border-radius: 2px;
        font-size: 12px !important;
      }
      .ltr-yes {
        color: var(--color-trf) !important;
        font-weight: 700;
        background: var(--color-trf-light);
      }
      .ltr-no {
        color: #999999 !important;
      }
      .gag-tag {
        padding: 1px 7px;
        border-radius: 2px;
        font-size: 12px !important;
      }
      .gag-yes {
        color: #6b3a6b !important;
        font-weight: 700;
        background: #f3ecf3;
      }
      .gag-no {
        color: #999999 !important;
      }
      .imprinted-tag {
        padding: 1px 7px;
        border-radius: 2px;
        font-size: 12px !important;
      }
      .imprinted-yes {
        color: var(--ink) !important;
        font-weight: 700;
      }
      .imprinted-no {
        color: #999999 !important;
      }
      .gencode-tag, .stringtie-tag {
        padding: 1px 7px;
        border-radius: 2px;
        font-size: 12px !important;
      }
      .gencode-genic {
        color: #4a6f8a !important;
        font-weight: 700;
        background: #edf3f8;
      }
      .gencode-intergenic {
        color: #999999 !important;
      }
      .stringtie-genic {
        color: #7a6520 !important;
        font-weight: 700;
        background: #fdf6e3;
      }
      .stringtie-intergenic {
        color: #999999 !important;
      }

      /* ── Alignment block ── */
      .alignment-block {
        background: #f9f9f9;
        border: 1px solid #e5e5e5;
        padding: 12px 14px;
        border-radius: 6px;
        font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
        font-size: 13px;
        line-height: 1.65;
        overflow-x: auto;
        margin: 0;
      }
      .alignment-grid {
        display: grid;
        grid-template-columns: max-content 1fr;
        column-gap: 14px;
        row-gap: 0;
        align-items: start;
        min-width: max-content;
      }
      .alignment-label,
      .alignment-seq {
        white-space: pre;
      }
      .alignment-label {
        color: var(--muted);
        font-weight: 600;
      }
      .alignment-label-query {
        color: var(--ink);
        font-weight: 700;
      }
      .alignment-label-target {
        color: var(--ink);
        font-weight: 700;
      }
      .wc { color: var(--color-gencode); font-weight: bold; }
      .gu { color: var(--color-target); font-weight: bold; }
      .mm { color: var(--color-trf); }

      .empty-state {
        background: #fafafa;
        border: 1px dashed var(--border);
        border-radius: 8px;
        padding: 20px;
        color: var(--muted);
      }

      /* ── Page controls ── */
      .page-controls {
        display: flex;
        gap: 8px;
        align-items: center;
        flex-wrap: wrap;
        margin-bottom: 12px;
      }
      .page-controls .btn {
        background: #ffffff;
        border: 1px solid var(--border);
        color: var(--ink);
        border-radius: 6px;
        font-size: 12.5px;
        font-weight: 600;
        padding: 4px 12px;
        transition: border-color 0.15s, background 0.15s;
      }
      .page-controls .btn:hover {
        border-color: var(--accent);
        color: var(--accent);
        background: var(--color-trf-light);
      }
      .page-label {
        color: var(--muted);
        font-size: 12.5px;
        font-variant-numeric: tabular-nums;
      }
      .results-hint,
      .sort-status {
        color: #999999;
        font-size: 12px;
        margin-bottom: 6px;
      }
      .coord-status {
        color: var(--muted);
        font-size: 11.5px;
        margin-top: -6px;
        margin-bottom: 6px;
        min-height: 14px;
        font-style: italic;
      }
      .coord-status:empty {
        display: none;
      }

      /* ── Table ── */
      table.dataTable {
        font-size: 12.5px;
        border-collapse: collapse;
      }
      table.dataTable thead th {
        cursor: pointer;
        font-weight: 700;
        font-size: 11px;
        letter-spacing: 0.03em;
        text-transform: uppercase;
        color: var(--muted);
        border-bottom: 1px solid var(--border) !important;
        padding: 7px 8px !important;
      }
      table.dataTable tbody td {
        padding: 5px 8px !important;
        border-bottom: 1px solid #f0f0f0 !important;
      }
      table.dataTable tbody tr:hover {
        background: #f7f7f7 !important;
      }
      table.dataTable tbody tr.selected,
      table.dataTable tbody tr.selected td,
      table.dataTable tbody > tr.selected,
      table.dataTable tbody > tr.selected td,
      table.dataTable tbody > tr > td.selected,
      .dataTable tbody tr.active,
      .dataTable tbody tr.active td {
        background: #9cb7d4 !important;
        background-color: #9cb7d4 !important;
        color: #ffffff !important;
      }
      table.dataTable thead th.sorted-asc::after {
        content: ' \\25B2';
        color: var(--accent);
        font-size: 9px;
      }
      table.dataTable thead th.sorted-desc::after {
        content: ' \\25BC';
        color: var(--accent);
        font-size: 9px;
      }

      /* ── Panel header with export ── */
      .panel-header {
        display: flex;
        justify-content: space-between;
        align-items: baseline;
      }
      .panel-header .panel-title {
        margin-bottom: 0;
      }
      #download_csv {
        background: #ffffff;
        border: 1px solid var(--border);
        color: var(--ink);
        border-radius: 6px;
        font-size: 12px;
        font-weight: 600;
        padding: 4px 14px;
        transition: border-color 0.15s, background 0.15s;
      }
      #download_csv:hover {
        border-color: var(--accent);
        color: var(--accent);
        background: var(--color-trf-light);
      }

      /* ── Form inputs ── */
      .form-control, .selectize-input {
        border-radius: 6px !important;
        border: 1px solid #cccccc !important;
        font-size: 13px !important;
      }
      .selectize-control.multi .selectize-input,
      .selectize-control.single .selectize-input,
      .selectize-input,
      .selectize-input.items {
        border: 1px solid #cccccc !important;
        box-shadow: none !important;
      }
      .selectize-dropdown {
        border: 1px solid #cccccc !important;
      }
      .form-control:focus, .selectize-input.focus {
        border-color: var(--color-gencode) !important;
        box-shadow: 0 0 0 2px rgba(54, 100, 139, 0.12) !important;
      }
      .irs-grid-text.small {
        display: none;
      }
      .irs--shiny .irs-grid-pol.small {
        height: 4px;
      }
      .irs--shiny .irs-bar {
        background: #9cb7d4 !important;
        border-top-color: #9cb7d4 !important;
        border-bottom-color: #9cb7d4 !important;
      }
      .irs--shiny .irs-handle {
        border-color: #9cb7d4 !important;
      }
      .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
        background-color: #9cb7d4 !important;
      }
      label {
        font-size: 12.5px;
        font-weight: 600;
        color: var(--ink);
      }
      hr {
        border-top: 1px solid #eeeeee;
      }
    "))
  ),
  shiny::div(
    class = "hero",
    shiny::div(
      class = "hero-inner",
      shiny::tags$h1(
        class = "hero-title",
        shiny::HTML("3&#8242;-tRF Target Explorer")
      ),
      shiny::p(
        class = "hero-subtitle",
        shiny::HTML(paste0(
          "Interactive browser for annotated 3&#8242;-tRF target sites. ",
          "Companion to: <em>3&#8242;-tRNA Fragments target domesticated LTR-Retrotransposons</em>."
        ))
      )
    )
  ),
  shiny::fluidRow(
    shiny::column(
      width = 3,
      shiny::div(
        class = "control-panel",
        shiny::div(class = "panel-title", "Filter Target Sites"),
        shiny::textInput("search", "Search Sites", placeholder = "tRF, gene, LTR..."),
        shiny::textInput("coord_range", "Genomic Range", placeholder = "chr1:1,000,000-2,000,000"),
        shiny::div(class = "coord-status", shiny::textOutput("coord_status")),
        shiny::sliderInput("score_range", "Alignment Score", min = score_min, max = score_max, value = c(score_min, score_max), step = 5, ticks = TRUE),
        shiny::selectInput("ltr_filter", "LTR Association", choices = c("All sites" = "all", "LTR-associated only" = "ltr", "Non-LTR only" = "no-ltr")),
        shiny::selectInput("gag_filter", "Gag-Derived", choices = c("All sites" = "all", "Gag-derived only" = "gag", "Non-Gag only" = "no-gag")),
        shiny::selectInput("imprinted_filter", "Imprinting", choices = c("All genes" = "all", "Imprinted only" = "imprinted", "Not imprinted only" = "non-imprinted")),
        shiny::selectizeInput("gencode_locations", "GENCODE Region", choices = gencode_location_choices, multiple = TRUE, options = list(placeholder = "All GENCODE regions")),
        shiny::selectizeInput("stringtie_locations", "StringTie Region", choices = stringtie_location_choices, multiple = TRUE, options = list(placeholder = "All StringTie regions")),
        shiny::hr()
      )
    ),
    shiny::column(
      width = 9,
      shiny::uiOutput("alignment_detail"),
      shiny::div(
        class = "results-panel",
        shiny::div(
          class = "panel-header",
          shiny::div(class = "panel-title", "Browse Alignments"),
          shiny::downloadButton("download_csv", "Export CSV")
        ),
        shiny::div(class = "results-summary", shiny::textOutput("results_summary")),
        shiny::div(class = "sort-status", shiny::textOutput("sort_status")),
        shiny::div(
          class = "page-controls",
          shiny::actionButton("prev_page", "Previous"),
          shiny::actionButton("next_page", "Next"),
          shiny::span(class = "page-label", shiny::textOutput("page_label", inline = TRUE))
        ),
        DT::dataTableOutput("results_table")
      )
    )
  )
)

server <- function(input, output, session) {
  current_page <- shiny::reactiveVal(1L)
  current_page_data <- shiny::reactiveVal(data.frame())

  parsed_coord_range <- shiny::debounce(
    shiny::reactive({
      parse_coord_range(input$coord_range)
    }),
    400
  )

  output$coord_status <- shiny::renderText({
    format_coord_status(parsed_coord_range())
  })

  shiny::observeEvent(
    list(input$search, input$ltr_filter, input$gag_filter, input$imprinted_filter, input$gencode_locations, input$stringtie_locations, input$score_range, parsed_coord_range(), input$results_table_sort),
    {
      current_page(1L)
    },
    ignoreInit = TRUE
  )

  filtered_count <- shiny::reactive({
    where_sql <- build_where_clause(
      con = con,
      search = input$search,
      ltr_filter = input$ltr_filter,
      gag_filter = input$gag_filter,
      imprinted_filter = input$imprinted_filter,
      gencode_locations = input$gencode_locations,
      stringtie_locations = input$stringtie_locations,
      score_min = input$score_range[1],
      score_max = input$score_range[2],
      coord_range = parsed_coord_range()
    )

    sql <- paste(
      "SELECT COUNT(*) AS n",
      "FROM targets",
      where_sql
    )

    DBI::dbGetQuery(con, sql)$n[[1]]
  })

  total_pages <- shiny::reactive({
    max(1L, ceiling(filtered_count() / as.integer(default_page_size)))
  })

  shiny::observe({
    if (current_page() > total_pages()) {
      current_page(total_pages())
    }
  })

  shiny::observeEvent(input$prev_page, {
    current_page(max(1L, current_page() - 1L))
  })

  shiny::observeEvent(input$next_page, {
    current_page(min(total_pages(), current_page() + 1L))
  })

  page_query <- shiny::reactive({
    page_size <- as.integer(default_page_size)
    offset <- (current_page() - 1L) * page_size

    where_sql <- build_where_clause(
      con = con,
      search = input$search,
      ltr_filter = input$ltr_filter,
      gag_filter = input$gag_filter,
      imprinted_filter = input$imprinted_filter,
      gencode_locations = input$gencode_locations,
      stringtie_locations = input$stringtie_locations,
      score_min = input$score_range[1],
      score_max = input$score_range[2],
      coord_range = parsed_coord_range()
    )

    sql <- paste(
      "SELECT",
      "  report_row_id, tRF, tRNA_anticodon, seqnames, start, \"end\", strand,",
      "  alignment_score, energy,",
      "  COALESCE(gencode_gene_name, '—') AS gencode_gene_name,",
      "  COALESCE(gencode_location, '—') AS gencode_location,",
      "  COALESCE(stringtie_location, '—') AS stringtie_location,",
      "  CASE WHEN LTR THEN COALESCE(LTR_family, 'LTR') ELSE 'None' END AS ltr_label,",
      "  CASE WHEN gag_gene THEN 'Yes' ELSE 'No' END AS gag_label,",
      "  COALESCE(imprint_status, '—') AS imprint_status",
      "FROM targets",
      where_sql,
      sort_clause(input$results_table_sort),
      sprintf("LIMIT %d OFFSET %d", page_size, offset)
    )

    DBI::dbGetQuery(con, sql)
  })

  shiny::observe({
    current_page_data(page_query())
  })

  output$results_summary <- shiny::renderText({
    total <- filtered_count()
    if (total == 0) {
      "No target sites match the current filters."
    } else {
      sprintf("%s target sites match the current filters.", format(total, big.mark = ","))
    }
  })

  output$page_label <- shiny::renderText({
    total <- filtered_count()
    if (total == 0) {
      "Page 1 of 1"
    } else {
      page_size <- as.integer(default_page_size)
      start_idx <- (current_page() - 1L) * page_size + 1L
      end_idx <- min(current_page() * page_size, total)
      sprintf(
        "Showing %s-%s of %s (page %s of %s)",
        format(start_idx, big.mark = ","),
        format(end_idx, big.mark = ","),
        format(total, big.mark = ","),
        current_page(),
        total_pages()
      )
    }
  })

  output$sort_status <- shiny::renderText({
    sort_label(input$results_table_sort)
  })

  output$results_table <- DT::renderDataTable({
    page_df <- page_query()
    if (nrow(page_df) == 0L) {
      return(DT::datatable(data.frame()))
    }

    display_df <- page_df
    display_df$location <- paste0(
      display_df$seqnames, ":",
      format(display_df$start, big.mark = ","),
      "-",
      format(display_df$end, big.mark = ","),
      " (", display_df$strand, ")"
    )
    display_df <- display_df[, c(
      "report_row_id",
      "tRF",
      "tRNA_anticodon",
      "location",
      "alignment_score",
      "energy",
      "gencode_gene_name",
      "gencode_location",
      "stringtie_location",
      "ltr_label",
      "gag_label",
      "imprint_status"
    )]
    names(display_df) <- c(
      "Row", "tRF", "tRNA Anticodon", "Location",
      "Alignment Score", "Energy", "Gene Name",
      "GENCODE Region", "StringTie Region",
      "LTR", "Gag", "Imprint Status"
    )

    current_sort <- normalize_table_sort(input$results_table_sort)
    current_col_js <- if (is.null(current_sort)) {
      "null"
    } else {
      as.character(current_sort$col)
    }
    current_dir_js <- if (is.null(current_sort)) {
      "null"
    } else {
      sprintf("'%s'", current_sort$dir)
    }

    table_widget <- DT::datatable(
      display_df,
      rownames = FALSE,
      selection = "single",
      options = list(
        dom = "t",
        scrollX = TRUE,
        scrollY = "60vh",
        paging = FALSE,
        info = FALSE,
        lengthChange = FALSE,
        searching = FALSE,
        ordering = FALSE
      ),
      callback = DT::JS(sprintf(
        "
        var currentCol = %s;
        var currentDir = %s;
        var headers = $(table.table().header()).find('th');

        headers.removeClass('sorted-asc sorted-desc');
        if (currentCol !== null && currentDir !== null) {
          headers.eq(currentCol).addClass(currentDir === 'desc' ? 'sorted-desc' : 'sorted-asc');
        }

        headers.off('click.codexSort').on('click.codexSort', function() {
          var col = table.column(this).index();
          var payload;

          if ($(this).hasClass('sorted-asc')) {
            payload = {col: col, dir: 'desc', mode: 'sort', nonce: Date.now()};
          } else if ($(this).hasClass('sorted-desc')) {
            payload = {mode: 'original', nonce: Date.now()};
          } else {
            payload = {col: col, dir: 'asc', mode: 'sort', nonce: Date.now()};
          }

          Shiny.setInputValue('results_table_sort', payload, {priority: 'event'});
        });
        ",
        current_col_js,
        current_dir_js
      ))
    )
    table_widget <- DT::formatRound(table_widget, "Alignment Score", 0)
    DT::formatRound(table_widget, "Energy", 2)
  }, server = FALSE)

  selected_row_id <- shiny::reactive({
    page_df <- current_page_data()
    if (nrow(page_df) == 0L) {
      return(NULL)
    }

    selected <- input$results_table_rows_selected
    if (length(selected) == 0L) {
      return(page_df$report_row_id[[1]])
    }

    page_df$report_row_id[[selected[1]]]
  })

  selected_row <- shiny::reactive({
    row_id <- selected_row_id()
    if (is.null(row_id)) {
      return(NULL)
    }

    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT report_row_id, seqnames, start, \"end\", strand, tRF, tRNA_anticodon,",
          "alignment_score, energy, LTR, LTR_family, LTR_gene_id, gag_gene,",
          "gencode_gene_name, gencode_location, stringtie_location, imprint_status,",
          "query_alignment, match_string, ref_alignment",
          "FROM targets WHERE report_row_id = %d"
        ),
        as.integer(row_id)
      )
    )
  })

  output$alignment_detail <- shiny::renderUI({
    build_alignment_card(selected_row())
  })

  output$download_csv <- shiny::downloadHandler(
    filename = function() {
      parts <- c("trf_targets")

      if (!is.null(input$search) && nzchar(trimws(input$search))) {
        slug <- tolower(gsub("[^A-Za-z0-9]+", "-", trimws(input$search)))
        slug <- sub("^-+|-+$", "", slug)
        parts <- c(parts, slug)
      }

      if (identical(input$ltr_filter, "ltr")) {
        parts <- c(parts, "LTR")
      } else if (identical(input$ltr_filter, "no-ltr")) {
        parts <- c(parts, "noLTR")
      }

      if (identical(input$gag_filter, "gag")) {
        parts <- c(parts, "gag")
      } else if (identical(input$gag_filter, "no-gag")) {
        parts <- c(parts, "noGag")
      }

      if (identical(input$imprinted_filter, "imprinted")) {
        parts <- c(parts, "imprinted")
      } else if (identical(input$imprinted_filter, "non-imprinted")) {
        parts <- c(parts, "notImprinted")
      }

      if (!is.null(input$gencode_locations) && length(input$gencode_locations) > 0L) {
        gc_slug <- tolower(gsub("[^A-Za-z0-9]+", "-", paste(input$gencode_locations, collapse = "-")))
        parts <- c(parts, paste0("gencode-", gc_slug))
      }

      if (!is.null(input$stringtie_locations) && length(input$stringtie_locations) > 0L) {
        st_slug <- tolower(gsub("[^A-Za-z0-9]+", "-", paste(input$stringtie_locations, collapse = "-")))
        parts <- c(parts, paste0("stringtie-", st_slug))
      }

      score_range <- input$score_range
      if (!is.null(score_range)) {
        score_lo <- score_range[1]
        score_hi <- score_range[2]
        if (score_lo > score_min || score_hi < score_max) {
          parts <- c(parts, sprintf("score%s-%s", score_lo, score_hi))
        }
      }

      coord_parsed <- parsed_coord_range()
      if (!is.null(coord_parsed) && is.null(coord_parsed$error)) {
        coord_slug <- tolower(gsub(
          "[^A-Za-z0-9]+",
          "-",
          sprintf("%s-%s-%s", coord_parsed$seqname, coord_parsed$start, coord_parsed$end)
        ))
        coord_slug <- sub("^-+|-+$", "", coord_slug)
        parts <- c(parts, paste0("range-", coord_slug))
      }

      if (length(parts) == 1L) {
        parts <- c(parts, "all")
      }

      paste0(paste(parts, collapse = "_"), ".csv")
    },
    content = function(file) {
      where_sql <- build_where_clause(
        con = con,
        search = input$search,
        ltr_filter = input$ltr_filter,
        gag_filter = input$gag_filter,
        imprinted_filter = input$imprinted_filter,
        gencode_locations = input$gencode_locations,
        stringtie_locations = input$stringtie_locations,
        score_min = input$score_range[1],
        score_max = input$score_range[2],
        coord_range = parsed_coord_range()
      )

      sql <- paste(
        "SELECT",
        "  tRF, tRNA_anticodon, seqnames, start, \"end\", strand,",
        "  alignment_score, energy,",
        "  LTR, LTR_family, LTR_gene_id, gag_gene,",
        "  gencode_gene_name, gencode_location,",
        "  stringtie_location, imprint_status",
        "FROM targets",
        where_sql,
        sort_clause(input$results_table_sort)
      )

      result <- DBI::dbGetQuery(con, sql)
      utils::write.csv(result, file, row.names = FALSE, na = "")
    }
  )
}

app <- shiny::shinyApp(ui = ui, server = server)
