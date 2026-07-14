#' Create virus feature annotations
#'
#' @param x A data frame with feature intervals, or \code{NULL} when columns are supplied separately.
#' @param feature Feature names.
#' @param start Feature start positions on the virus genome.
#' @param end Feature end positions on the virus genome.
#' @param type Optional feature type labels.
#' @return A \code{vi_virus_features} object.
#' @export
virus_features <- function(x = NULL, feature = NULL, start = NULL, end = NULL, type = NULL) {
  if (is.null(x)) {
    x <- data.frame(
      feature = feature,
      start = start,
      end = end,
      stringsAsFactors = FALSE
    )
    if (!is.null(type)) {
      x$type <- type
    }
  }

  if (!is.data.frame(x)) {
    stop("x must be a data frame or NULL when feature/start/end are supplied.", call. = FALSE)
  }

  df <- as.data.frame(x, stringsAsFactors = FALSE)
  rename_one <- function(target, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0) {
      return(NULL)
    }
    df[[target]] <<- df[[hit[1]]]
    invisible(NULL)
  }

  rename_one("feature", c("feature", "gene", "name", "label"))
  rename_one("start", c("start", "start_pos", "from"))
  rename_one("end", c("end", "stop", "end_pos", "to"))
  rename_one("type", c("type", "feature_type"))

  required_cols <- c("feature", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Virus features must contain columns: feature, start, end", call. = FALSE)
  }

  if (!"type" %in% colnames(df)) {
    df$type <- "feature"
  }

  out <- df[, c("feature", "start", "end", "type"), drop = FALSE]
  out$feature <- as.character(out$feature)
  out$start <- suppressWarnings(as.numeric(out$start))
  out$end <- suppressWarnings(as.numeric(out$end))
  out$type <- as.character(out$type)

  keep <- !is.na(out$feature) & nzchar(out$feature) &
    !is.na(out$start) & !is.na(out$end) & out$end > out$start
  if (!all(keep)) {
    warning(sum(!keep), " virus features were dropped because they were incomplete or invalid.")
  }
  out <- out[keep, , drop = FALSE]

  if (nrow(out) == 0) {
    stop("No valid virus features were found.", call. = FALSE)
  }

  out <- out[order(out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("vi_virus_features", "data.frame")
  out
}

#' Built-in HBV feature annotations
#'
#' @param version HBV annotation version. Currently only \code{"ayw"} is supported.
#' @return A \code{vi_virus_features} object.
#' @export
hbv_features <- function(version = "ayw") {
  version <- match.arg(version, "ayw")
  virus_features(data.frame(
    feature = c("preS1", "preS2", "S", "X", "preC", "C", "P", "P"),
    start = c(2850, 1, 155, 1374, 1814, 1901, 2307, 1),
    end = c(3204, 154, 835, 1838, 1900, 2452, 3215, 1623),
    type = c("surface", "surface", "surface", "regulatory", "core", "core", "polymerase", "polymerase"),
    stringsAsFactors = FALSE
  ))
}

normalize_virus_features <- function(features) {
  if (inherits(features, "vi_virus_features")) {
    return(virus_features(features))
  }
  virus_features(features)
}

layout_virus_features <- function(features, virus_length = NULL) {
  df <- as.data.frame(normalize_virus_features(features), stringsAsFactors = FALSE)
  if (!is.null(virus_length)) {
    wrap <- df$end <= df$start
    df$end[wrap] <- df$end[wrap] + virus_length
  }

  df$level <- 1L
  level_ends <- numeric()
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    for (level in seq_along(level_ends)) {
      if (df$start[i] >= level_ends[level]) {
        df$level[i] <- level
        level_ends[level] <- df$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      level_ends <- c(level_ends, df$end[i])
      df$level[i] <- length(level_ends)
    }
  }

  df
}
