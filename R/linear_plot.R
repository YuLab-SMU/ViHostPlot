#' Plot integrations on linear host-virus axes
#'
#' @param integrations Integration records.
#' @param host Host genome specification. Can be a built-in host name, a host
#'   genome table, or a \code{vi_host_genome} object.
#' @param virus Virus genome specification, typically a named numeric vector of
#'   length 1 such as \code{c(HBV = 3215)}.
#' @param features Optional virus feature intervals from \code{virus_features()}.
#' @param chrom_file Optional path to a host chromosome-size table.
#' @param host_scope Host chromosomes to show. \code{"hit"} shows only chromosomes
#'   with integrations; \code{"all"} shows all host chromosomes.
#' @param link_color Link color.
#' @param point_color Point color.
#' @param feature_fill Virus feature fill color or named vector by feature type.
#' @return A \code{ggplot} object.
#' @export
plot_linear_integrations <- function(integrations, host = "hg38", virus,
                                     features = NULL, chrom_file = NULL,
                                     host_scope = c("hit", "all"),
                                     link_color = "grey70",
                                     point_color = "#2C7FB8",
                                     feature_fill = NULL) {
  host_scope <- match.arg(host_scope)
  integrations <- as_integrations(integrations)
  virus_obj <- virus_genome(virus)
  host_obj <- if (inherits(host, "vi_host_genome")) {
    host
  } else {
    host_genome(host = host, chrom_file = chrom_file)
  }

  host_df <- host_obj$data
  df <- as.data.frame(integrations, stringsAsFactors = FALSE)
  df <- df[df$host_chr %in% host_df$chr, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No integration records matched the host genome.", call. = FALSE)
  }

  if (host_scope == "hit") {
    host_df <- host_df[host_df$chr %in% unique(df$host_chr), , drop = FALSE]
  }

  host_df <- make_linear_host_axis(host_df)
  df <- merge(
    df,
    host_df[, c("chr", "offset"), drop = FALSE],
    by.x = "host_chr",
    by.y = "chr",
    all.x = TRUE,
    sort = FALSE
  )
  df$host_x <- df$offset + df$host_pos

  host_span <- max(host_df$xend, na.rm = TRUE)
  df$virus_x <- df$virus_pos / virus_obj$length * host_span

  host_rect <- data.frame(
    xmin = host_df$xstart,
    xmax = host_df$xend,
    ymin = 0.88,
    ymax = 1.12,
    label = host_df$chr,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = host_rect,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = "grey92",
      color = "white"
    ) +
    ggplot2::geom_segment(
      data = df,
      ggplot2::aes(x = .data$host_x, xend = .data$virus_x, y = 1, yend = 0),
      color = link_color,
      linewidth = 0.25,
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$host_x, y = 1),
      color = point_color,
      size = 1.5
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$virus_x, y = 0),
      color = point_color,
      size = 1.5
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 1),
      labels = c(virus_obj$name, "Host"),
      limits = c(-0.3, 1.3),
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  if (!is.null(features)) {
    feature_df <- layout_virus_features(features, virus_length = virus_obj$length)
    feature_df$xmin <- pmax(0, pmin(feature_df$start, virus_obj$length)) / virus_obj$length * host_span
    feature_df$xmax <- pmax(0, pmin(feature_df$end, virus_obj$length)) / virus_obj$length * host_span
    feature_df$ymin <- -0.2 - (feature_df$level - 1) * 0.08
    feature_df$ymax <- feature_df$ymin + 0.06

    fill_map <- if (is.null(feature_fill)) {
      normalize_named_colors(grDevices::hcl.colors(length(unique(feature_df$type)), "Set 3"), unique(feature_df$type))
    } else {
      normalize_named_colors(feature_fill, unique(feature_df$type))
    }
    feature_df$fill <- unname(fill_map[feature_df$type])

    p <- p +
      ggplot2::geom_rect(
        data = feature_df,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
        fill = feature_df$fill,
        color = "white"
      ) +
      ggplot2::geom_text(
        data = feature_df,
        ggplot2::aes(x = (.data$xmin + .data$xmax) / 2, y = (.data$ymin + .data$ymax) / 2, label = .data$feature),
        size = 2.5
      )
  }

  p
}

make_linear_host_axis <- function(host_df) {
  host_df <- host_df[order(match(host_df$chr, mixed_chr_order(host_df$chr))), , drop = FALSE]
  lengths <- host_df$end - host_df$start
  host_df$offset <- c(0, cumsum(lengths)[-length(lengths)])
  host_df$xstart <- host_df$offset
  host_df$xend <- host_df$offset + lengths
  host_df
}

mixed_chr_order <- function(chr) {
  chr <- as.character(chr)
  numeric_chr <- suppressWarnings(as.integer(chr))
  chr[order(is.na(numeric_chr), numeric_chr, chr)]
}
