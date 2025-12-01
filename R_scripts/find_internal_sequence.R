
find_internal_sequence <- function(LTR_df, internal_name, max_gap = 200) {
  
  # 1) Split: non-internal LTR vs internal hits
  LTR_nonint <- LTR_df[!grepl("int", LTR_df$repName, ignore.case = TRUE), ]
  LTR_int    <- LTR_df[grepl(internal_name, LTR_df$repName, ignore.case = FALSE), ]
  
  # 2) GRanges
  query <- makeGRangesFromDataFrame(
    LTR_nonint,
    keep.extra.columns = TRUE,
    seqnames.field = "genoName",
    start.field    = "genoStart",
    end.field      = "genoEnd",
    strand.field   = "strand"
  )
  
  subject <- makeGRangesFromDataFrame(
    LTR_int,
    keep.extra.columns = TRUE,
    seqnames.field = "genoName",
    start.field    = "genoStart",
    end.field      = "genoEnd",
    strand.field   = "strand"
  )
  
  # 3) Upstream windows (strand-aware)
  sub_plus  <- subject[strand(subject) == "+"]
  sub_minus <- subject[strand(subject) == "-"]
  
  win_plus <- GRanges(
    seqnames = seqnames(sub_plus),
    ranges   = IRanges(start = pmax(start(sub_plus) - max_gap, 1),
                       end   = pmax(start(sub_plus) - 1, 0)),
    strand   = strand(sub_plus)
  )
  
  win_minus <- GRanges(
    seqnames = seqnames(sub_minus),
    ranges   = IRanges(start = end(sub_minus) + 1,
                       end   = end(sub_minus) + max_gap),
    strand   = strand(sub_minus)
  )
  
  # Keep ETn indices
  win_plus$etn_idx  <- seq_along(sub_plus)
  win_minus$etn_idx <- seq_along(sub_minus)
  
  # 4) Overlaps
  hits_plus  <- findOverlaps(query, win_plus,  ignore.strand = TRUE)
  hits_minus <- findOverlaps(query, win_minus, ignore.strand = TRUE)
  
  # 5) Distances
  df_plus <- as.data.frame(hits_plus) %>%
    mutate(etn_idx = win_plus$etn_idx[subjectHits],
           dist_bp = start(sub_plus)[etn_idx] - end(query[queryHits]))
  
  df_minus <- as.data.frame(hits_minus) %>%
    mutate(etn_idx = win_minus$etn_idx[subjectHits],
           dist_bp = start(query[queryHits]) - end(sub_minus)[etn_idx])
  
  # 6) Nearest upstream LTR per ETn
  pick_nearest <- function(df) {
    if (nrow(df) == 0) return(df)
    df %>% group_by(etn_idx) %>%
      slice_min(order_by = dist_bp, with_ties = FALSE) %>%
      ungroup()
  }
  df_near <- bind_rows(pick_nearest(df_plus), pick_nearest(df_minus))
  
  # 7) Counts
  if (nrow(df_near) > 0) {
    nearest_LTR_rep <- mcols(query)$repName[df_near$queryHits]
    counts <- as.data.frame(table(nearest_LTR_rep), stringsAsFactors = FALSE) %>%
      dplyr::arrange(dplyr::desc(Freq)) %>%
      dplyr::rename(LTR_repName = nearest_LTR_rep, n_upstream_hits = Freq)
  } else {
    counts <- tibble(LTR_repName = character(), n_upstream_hits = integer())
  }
  
  return(counts)
}

home_directory = '~/CSHL Dropbox Team Dropbox/Matthew Peacey/SchornLab (1)/manuscripts/Domesticated_targets'

LTR = read.csv(glue('{home_directory}/Figure_1/UCSC_mm10_LTR.csv'),
                header = TRUE, stringsAsFactors = FALSE)

LTR$repFamily = LTR$repFamily %>%
  str_remove("\\?") %>% 
  str_trim()

unique(LTR[grepl("int", LTR$repName, ignore.case = TRUE), ]$repName)

find_internal_sequence(LTR, "MMVL30-int")