#' Plot host-virus breakpoint coordinates
#'
#' Display integration records on a two-dimensional map where the x-axis is the
#' concatenated host genome and the y-axis is the virus genome. This view is
#' useful for comparing host breakpoint distribution with recurrent viral
#' breakpoint positions.
#'
#' @param integrations Integration records accepted by [as_integrations()].
#' @param host Host genome specification. Can be a built-in host name, a host
#'   genome table, or a [host_genome()] object.
#' @param virus Virus genome specification. Can be a [virus_genome()] object or
#'   a named numeric vector such as `c(HBV = 3215)`.
#' @param features Optional virus feature intervals from [virus_features()].
#' @param chrom_file Optional path to a chromosome-size table. Supply either
#'   `host` or `chrom_file`, not both.
#' @param group_by Optional column in `integrations` mapped to point color.
#' @param size_by Optional numeric column in `integrations` mapped to point size.
#' @param colors Optional vector of colors. A named vector can be used to set
#'   colors for specific groups.
#' @param point_size Point size used when `size_by` is `NULL`.
#' @param size_range Point-size range used when `size_by` is supplied.
#' @param alpha Point transparency.
#' @return A `ggplot` object.
#' @export
plot_breakpoint_map <- function(integrations, host = "hg38", virus,
                                features = NULL, chrom_file = NULL,
                                group_by = NULL, size_by = NULL,
                                colors = NULL, point_size = 2.5,
                                size_range = c(1.5, 5),
                                alpha = 0.75) {
  integrations <- as.data.frame(as_integrations(integrations), stringsAsFactors = FALSE)

  validate_single_column_arg(group_by, integrations, "group_by", allow_null = TRUE)
  validate_single_column_arg(size_by, integrations, "size_by", allow_null = TRUE)
  if (!is.null(size_by) && !is.numeric(integrations[[size_by]])) {
    stop("size_by must refer to a numeric column.", call. = FALSE)
  }
  if (!is.numeric(point_size) || length(point_size) != 1L || is.na(point_size) ||
      point_size <= 0) {
    stop("point_size must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(size_range) || length(size_range) != 2L ||
      any(is.na(size_range)) || any(size_range <= 0) || size_range[1] > size_range[2]) {
    stop("size_range must contain two positive numbers in increasing order.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a single number between 0 and 1.", call. = FALSE)
  }

  virus_obj <- if (inherits(virus, "vi_virus_genome")) {
    virus
  } else {
    virus_genome(virus)
  }
  host_df <- prepare_breakpoint_host_axis(host = host, chrom_file = chrom_file)

  chr_index <- match(integrations$host_chr, host_df$chr)
  in_host <- !is.na(chr_index)
  in_host_range <- in_host
  in_host_range[in_host] <- integrations$host_pos[in_host] >= host_df$start[chr_index[in_host]] &
    integrations$host_pos[in_host] <= host_df$end[chr_index[in_host]]
  in_virus_range <- integrations$virus_pos >= 0 & integrations$virus_pos <= virus_obj$length
  keep <- in_host & in_host_range & in_virus_range
  if (!all(keep)) {
    warning(sum(!keep), " integration records were outside the host or virus genome and were omitted.",
            call. = FALSE)
  }

  integrations <- integrations[keep, , drop = FALSE]
  chr_index <- chr_index[keep]
  if (nrow(integrations) == 0L) {
    stop("No integration records matched the host and virus genomes.", call. = FALSE)
  }

  integrations$genome_pos <- host_df$offset[chr_index] +
    integrations$host_pos - host_df$start[chr_index]
  integrations$chromosome <- factor(integrations$host_chr, levels = host_df$chr, ordered = TRUE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = host_df,
      ggplot2::aes(xmin = .data$genome_start, xmax = .data$genome_end,
                   ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = host_df$band,
      color = NA
    )

  if (!is.null(features)) {
    feature_df <- prepare_breakpoint_features(features, virus_obj$length)
    if (nrow(feature_df) > 0L) {
      p <- p +
        ggplot2::geom_rect(
          data = feature_df,
          ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = .data$start, ymax = .data$end),
          inherit.aes = FALSE,
          fill = feature_df$fill,
          alpha = 0.18,
          color = NA
        )
    }
  }

  point_mapping <- ggplot2::aes(x = .data$genome_pos, y = .data$virus_pos)
  if (!is.null(group_by) && !is.null(size_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$genome_pos,
      y = .data$virus_pos,
      color = .data[[group_by]],
      size = .data[[size_by]]
    )
  } else if (!is.null(group_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$genome_pos,
      y = .data$virus_pos,
      color = .data[[group_by]]
    )
  } else if (!is.null(size_by)) {
    point_mapping <- ggplot2::aes(
      x = .data$genome_pos,
      y = .data$virus_pos,
      size = .data[[size_by]]
    )
  }

  if (is.null(group_by) && is.null(size_by)) {
    p <- p + ggplot2::geom_point(
      data = integrations,
      mapping = point_mapping,
      color = "#0072B2",
      size = point_size,
      alpha = alpha,
      na.rm = TRUE
    )
  } else if (is.null(group_by)) {
    p <- p + ggplot2::geom_point(
      data = integrations,
      mapping = point_mapping,
      color = "#0072B2",
      alpha = alpha,
      na.rm = TRUE
    )
  } else if (is.null(size_by)) {
    p <- p + ggplot2::geom_point(
      data = integrations,
      mapping = point_mapping,
      size = point_size,
      alpha = alpha,
      na.rm = TRUE
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = integrations,
      mapping = point_mapping,
      alpha = alpha,
      na.rm = TRUE
    )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = host_df$midpoint,
      labels = host_df$chr,
      expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, virus_obj$length),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(x = "Host chromosome", y = paste0(virus_obj$name, " position")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = if (is.null(group_by) && is.null(size_by)) "none" else "right"
    )

  if (!is.null(group_by)) {
    p <- p + ggplot2::labs(color = group_by)
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_color_manual(values = colors)
    }
  }
  if (!is.null(size_by)) {
    p <- p +
      ggplot2::scale_size_continuous(range = size_range) +
      ggplot2::labs(size = size_by)
  }

  p
}

validate_single_column_arg <- function(column, df, arg, allow_null = TRUE) {
  if (is.null(column) && allow_null) {
    return(invisible(TRUE))
  }
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column)) {
    stop(arg, " must be NULL or a single column name.", call. = FALSE)
  }
  if (!column %in% colnames(df)) {
    stop("Unknown ", arg, " column: ", column, call. = FALSE)
  }
  invisible(TRUE)
}

prepare_breakpoint_host_axis <- function(host = "hg38", chrom_file = NULL) {
  host_df <- resolve_host_chrom_sizes(host = host, chrom_file = chrom_file)
  host_df <- host_df[!duplicated(host_df$chr), , drop = FALSE]
  host_df$length <- host_df$end - host_df$start
  if (any(host_df$length <= 0)) {
    stop("Host chromosome end positions must be greater than start positions.", call. = FALSE)
  }
  host_df$offset <- c(0, utils::head(cumsum(host_df$length), -1L))
  host_df$genome_start <- host_df$offset
  host_df$genome_end <- host_df$offset + host_df$length
  host_df$midpoint <- (host_df$genome_start + host_df$genome_end) / 2
  host_df$band <- rep(c("#F2F2F2", "#FFFFFF"), length.out = nrow(host_df))
  host_df
}

prepare_breakpoint_features <- function(features, virus_length) {
  feature_df <- as.data.frame(normalize_virus_features(features), stringsAsFactors = FALSE)
  feature_df$start <- pmax(0, pmin(feature_df$start, virus_length))
  feature_df$end <- pmax(0, pmin(feature_df$end, virus_length))
  feature_df <- feature_df[feature_df$end > feature_df$start, , drop = FALSE]
  if (nrow(feature_df) == 0L) {
    return(feature_df)
  }
  fill_map <- normalize_named_colors(
    grDevices::hcl.colors(length(unique(feature_df$type)), "Set 3"),
    unique(feature_df$type)
  )
  feature_df$fill <- unname(fill_map[feature_df$type])
  feature_df
}
