test_that("plot_integrations can build first and draw once", {
  integrations <- data.frame(
    host_chr = c("1", "1", "2", "2"),
    host_pos = c(100, 600, 300, 800),
    virus_pos = c(50, 150, 450, 850),
    sample = c("A", "B", "A", "B"),
    support = c(10, 20, 15, 30),
    virus_strand = c("+", "-", "+", "-"),
    method = c("PB", "Pos", "PB", "Pos")
  )
  host <- host_genome(data.frame(
    chr = c("1", "2"),
    start = c(0, 0),
    end = c(1000, 1000)
  ))
  virus <- virus_genome(c(HBV = 1000))
  features <- virus_features(
    feature = c("geneA", "short"),
    start = c(50, 500),
    end = c(450, 560),
    type = c("gene", "gene")
  )

  p <- plot_integrations(
    integrations,
    host = host,
    virus = virus,
    tracks = list(
      track_ideogram(),
      track_sites(color = "method", label = "Sites"),
      track_virus_genes(features, label_min_width = 100),
      track_virus_density(label = "Virus density"),
      track_links(color = "method")
    ),
    draw = FALSE
  )

  expect_s3_class(p, "vi_integration_plot")
  expect_false(p$drawn)
  expect_named(
    integration_method_legend_spec(p$tracks, p$plot_df),
    c("PB", "Pos")
  )

  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 900, height = 700)
  on.exit({
    circlize::circos.clear()
    grDevices::dev.off()
  }, add = TRUE)

  expect_no_error(p <- draw_integration_plot(p))
  expect_true(p$drawn)
  expect_no_error(print(p))
})
