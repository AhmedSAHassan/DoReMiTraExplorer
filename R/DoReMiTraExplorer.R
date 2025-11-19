
# Helper: detect RNA-seq platform
is_rnaseq_platform <- function(se) {
  if (!is.null(colData(se)$Platform)) {
    platform_vals <- unique(as.character(colData(se)$Platform))
    platform_vals <- platform_vals[!is.na(platform_vals)]
    return(any(grepl("RNAseq|RNA-seq|rna", platform_vals, ignore.case = TRUE)))
  }
  FALSE
}

# Determine short data type label from Platform
determine_data_type <- function(se) {
  if (!is.null(colData(se)$Platform)) {
    pvals <- unique(as.character(colData(se)$Platform))
    pvals <- pvals[!is.na(pvals)]
    p <- tolower(pvals)
    if (length(p) > 0) {
      if (any(grepl("rna", p))) return("RNAseq")
      if (any(grepl("array|microarr|affymetrix|agilent", p))) return("Microarray")
    }
  }
  if (is_rnaseq_platform(se)) return("RNAseq")
  return("Unknown")
}

# return expression matrix for general use (VST for RNAseq, exprs for microarray)
get_expr_matrix <- function(se) {
  if (is_rnaseq_platform(se)) {
    exprs_raw <- tryCatch(as.matrix(assay(se)), error = function(e) NULL)
    if (is.null(exprs_raw)) stop("Could not extract count matrix from assay(se) for RNA-seq.")
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData = exprs_raw,
                                                   colData = colData(se),
                                                   design = ~1))
    vsd <- vst(dds, blind = TRUE)
    return(list(mat = assay(vsd), vsd = vsd, dds = dds))
  } else {
    if ("exprs" %in% names(assays(se))) {
      emat <- assays(se)$exprs
    } else {
      emat <- assay(se)
    }
    return(list(mat = as.matrix(emat), vsd = NULL, dds = NULL))
  }
}


#'Visualize radiation transcriptomic datasets in the form of SE objects from the DoReMiTra collection
#'
#' @param se SummarizedExperiment object from the DoReMiTra collection
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' se <- DoReMiTra::get_DoReMiTra_data("SE_Salah_2025_ExVivo", gene_symbol = TRUE)
#' DoReMiTra_explorer(se)
#
DoReMiTra_explorer <- function(se) {
  se_name <- deparse(substitute(se))

  expr_list <- get_expr_matrix(se)
  mat <- expr_list$mat
  vsd_obj <- expr_list$vsd
  dds_obj <- expr_list$dds

  gene_choices <- rownames(mat)

  platform_vals <- if (!is.null(colData(se)$Platform)) unique(as.character(colData(se)$Platform)) else NULL
  platform_vals <- platform_vals[!is.na(platform_vals)]
  platform_text <- if (length(platform_vals) == 0) "Unknown" else paste(platform_vals, collapse = ", ")
  is_rnaseq <- is_rnaseq_platform(se)
  data_type_text <- determine_data_type(se)  # short label based on Platform (RNAseq / Microarray / Unknown)

  # ggplot theme
  ggpt <- theme_bw() +
    theme(
      axis.title.x = element_text(size = rel(1.6)),
      axis.title.y = element_text(size = rel(1.6)),
      axis.text.x = element_text(size = rel(1.6)),
      axis.text.y = element_text(size = rel(1.6)),
      strip.text.x = element_text(size = rel(1.6)),
      strip.text.y = element_text(size = rel(1.6)),
      legend.text = element_text(size = rel(1.4)),
      legend.title = element_text(size = rel(1.4)),
      plot.title = element_text(size = rel(1.6))
    )

  ui <- dashboardPage(
    header = dashboardHeader(title = "DoReMiTra Explorer"),
    sidebar = dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("Summary", tabName = "summary", icon = icon("eye")),
        menuItem("PCA", tabName = "pca", icon = icon("chart-line")),
        menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
        menuItem("Boxplot", tabName = "boxplot", icon = icon("box-open")),
        menuItem("Gene Plot", tabName = "geneplot", icon = icon("dna"))
      ),
      conditionalPanel(
        condition = "input.tabs == 'heatmap'",
        hr(),
        h4("Heatmap Controls"),
        numericInput("topn", "Top variable genes:", value = 50, min = 5),
        numericInput("clusters", "Number of clusters:", value = 2, min = 1),
        checkboxInput("cluster_rows", "Cluster rows", value = TRUE),
        checkboxInput("cluster_cols", "Cluster columns", value = TRUE)
        # NOTE: heatmap row size
      ),
      conditionalPanel(
        condition = "input.tabs == 'pca'",
        hr(),
        pickerInput(
          "pca_color_by", "Color by:",
          choices = colnames(colData(se)),
          selected = if("Dose" %in% colnames(colData(se))) "Dose" else colnames(colData(se))[1]
        ),
        numericInput("pca_topn", "Number of top variable genes:", value = 500, min = 100)
      ),
      conditionalPanel(
        condition = "input.tabs == 'boxplot'",
        hr(),
        pickerInput("box_group_var", "Group by:",
                    choices = colnames(colData(se)),
                    selected = colnames(colData(se))[1])
      ),
      conditionalPanel(
        condition = "input.tabs == 'geneplot'",
        hr(),
        selectizeInput("gene", "Select gene(s) to plot:",
                       choices = NULL, multiple = TRUE, options = list(placeholder = 'Type gene name...')),
        pickerInput("geneplot_group", "Group by:",
                    choices = if("Dose" %in% colnames(colData(se))) "Dose" else colnames(colData(se)),
                    selected = if("Dose" %in% colnames(colData(se))) "Dose" else colnames(colData(se))[1]),
        checkboxInput("geneplot_labels", "Show sample labels", value = FALSE),
        uiOutput("genecard_links"),
        downloadButton("download_geneplot", "Download Gene Plot")
      )
    ),
    body = dashboardBody(
      tags$head(tags$style(HTML("
        .dorem-base-container { max-width: 980px; margin: 0 auto; }
        .dorem-base-container .panel { padding-top: 10px; padding-bottom: 20px; }
      "))),
      tabItems(
        tabItem(tabName = "summary", div(class = "dorem-base-container panel", uiOutput("summary_ui"))),
        tabItem(tabName = "pca", div(class = "dorem-base-container panel", plotOutput("pca_plot", height = "420px"))),
        tabItem(tabName = "heatmap", div(class = "dorem-base-container panel", plotOutput("heatmap_plot"))),
        tabItem(tabName = "boxplot", div(class = "dorem-base-container panel", plotOutput("boxplot", height = "420px"))),
        tabItem(tabName = "geneplot", div(class = "dorem-base-container panel", uiOutput("gene_plot_controls"), plotOutput("geneplots_combined")))
      )
    )
  )

  server <- function(input, output, session) {


    observe({
      updateSelectizeInput(session, "gene",
                           choices = gene_choices,
                           selected = NULL,
                           server = TRUE)
    })

    output$summary_ui <- renderUI({
      n_samples <- ncol(se)
      n_genes <- nrow(se)
      organism <- if ("Organism" %in% colnames(colData(se))) unique(as.character(colData(se)$Organism)) else "NA"
      radiation_type <- if ("Radiation_type" %in% colnames(colData(se))) unique(as.character(colData(se)$Radiation_type)) else "NA"
      exp_setting <- if ("Exp_setting" %in% colnames(colData(se))) unique(as.character(colData(se)$Exp_setting)) else "NA"
      tagList(
        tags$b("Dataset: "), se_name, tags$br(),
        tags$b("Platform: "), platform_text, tags$br(),
        tags$b("Organism(s): "), paste(organism, collapse = ", "), tags$br(),
        tags$b("Radiation Type: "), paste(radiation_type, collapse = ", "), tags$br(),
        tags$b("Experiment Setting: "), paste(exp_setting, collapse = ", "), tags$br(),
        tags$b("Number of Genes: "), n_genes, tags$br(),
        tags$b("Number of Samples: "), n_samples, tags$br()
      )
    })

    # PCA
    output$pca_plot <- renderPlot({
      topn <- min(as.numeric(input$pca_topn), nrow(mat))
      gene_vars <- matrixStats::rowVars(mat)
      top_genes <- order(gene_vars, decreasing = TRUE)[1:topn]
      mat_top <- mat[top_genes, , drop = FALSE]
      mat_top <- mat_top[apply(mat_top, 1, function(x) all(is.finite(x))), , drop = FALSE]
      pca_res <- prcomp(t(mat_top), scale. = TRUE)
      percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

      df <- data.frame(pca_res$x[, 1:2], colData(se))

      # determine valid color column name (or NULL)
      color_col <- if (!is.null(input$pca_color_by) && input$pca_color_by %in% colnames(colData(se))) {
        input$pca_color_by
      } else {
        NULL
      }

      if (!is.null(color_col)) {
        p <- ggplot(df, aes_string(x = "PC1", y = "PC2", color = color_col)) +
          geom_point(size = 3)
      } else {
        p <- ggplot(df, aes(x = PC1, y = PC2)) +
          geom_point(size = 3)
      }

      p + xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggpt
    }, res = 96)

    # Heatmap
    output$heatmap_plot <- renderPlot({

      req(input$topn)
      topn <- min(as.numeric(input$topn), nrow(mat))

      # Fixed pixels per row
      row_px <- 16L

      # compute dynamic height based on rows
      extra_h <- 140L
      raw_h <- topn * row_px + extra_h
      max_height <- 2400L
      plot_h <- min(max(420L, raw_h), max_height)

      # Prepare annotation
      dose <- if ("Dose" %in% colnames(colData(se))) factor(colData(se)$Dose) else factor(rep(NA, ncol(mat)))
      dose_levels <- levels(dose)
      if (length(dose_levels) == 0) {
        dose_colors <- structure(character(0), names = character(0))
      } else {
        pal_len <- max(3, length(dose_levels))
        dose_colors <- RColorBrewer::brewer.pal(pal_len, "Set1")[1:length(dose_levels)]
        names(dose_colors) <- dose_levels
      }
      ha <- ComplexHeatmap::HeatmapAnnotation(Dose = dose, col = list(Dose = dose_colors))

      var_genes_idx <- head(order(matrixStats::rowVars(mat), decreasing = TRUE), topn)
      var_mat <- mat[var_genes_idx, , drop = FALSE]

      # Enforced values
      fontsize_row <- 5
      show_row_names_flag <- TRUE

      # Automatically determine max rowname width
      max_name_length <- max(nchar(rownames(var_mat)))
      rowname_maxwidth_cm <- max(4, min(12, max_name_length * 0.12))

      # Detect draw() formals for compatibility
      draw_formals <- names(formals(ComplexHeatmap::draw))
      draw_supports_row_names_max_width <- "row_names_max_width" %in% draw_formals
      draw_supports_row_names_side <- "row_names_side" %in% draw_formals

      if (!draw_supports_row_names_max_width && show_row_names_flag) {
        approx_chars_per_cm <- 3.5
        max_chars <- max(3, floor(rowname_maxwidth_cm * approx_chars_per_cm))
        rn <- rownames(var_mat)
        rn_trunc <- ifelse(nchar(rn) > max_chars, paste0(substr(rn, 1, max_chars - 1), "…"), rn)
        rownames(var_mat) <- rn_trunc
      }

      # Prepare row splitting (k-means)
      k <- as.integer(if (!is.null(input$clusters)) input$clusters else 1)
      if (!is.na(k) && k > 1) {
        set.seed(123)
        km <- tryCatch(kmeans(var_mat, centers = k), error = function(e) NULL)
        if (!is.null(km) && length(km$cluster) == nrow(var_mat)) {
          row_split <- factor(km$cluster)
          row_title <- rep("", k)
        } else {
          row_split <- NULL
          row_title <- NULL
        }
      } else {
        row_split <- NULL
        row_title <- NULL
      }

      # Build heatmap
      ht <- ComplexHeatmap::Heatmap(
        var_mat - rowMeans(var_mat),
        name = "Scaled Expression",
        heatmap_width = unit(12, "cm"),
        row_names_max_width = unit(rowname_maxwidth_cm, "cm"),
        cluster_rows = input$cluster_rows,
        cluster_columns = input$cluster_cols,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_split = row_split,
        row_title = row_title,
        top_annotation = ha,
        show_row_names = show_row_names_flag,
        row_names_gp = grid::gpar(fontsize = fontsize_row),
        column_names_gp = grid::gpar(fontsize = 7),
        row_names_side = "left"
      )

      # Draw heatmap
      if (draw_supports_row_names_max_width) {
        ComplexHeatmap::draw(
          ht,
          heatmap_legend_side = "right",
          row_names_side = if (draw_supports_row_names_side) "left" else NULL,
          row_names_max_width = grid::unit(rowname_maxwidth_cm, "cm"),
          padding = unit(c(5, 10, 5, 10), "mm")
        )
      } else {
        ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
      }

    },
    height = function() {
      topn <- min(as.numeric(input$topn), nrow(mat))
      row_px <- 16L
      extra_h <- 140L
      raw_h <- topn * row_px + extra_h
      max_height <- 2400L
      plot_h <- min(max(420L, raw_h), max_height)
      plot_h
    },
    width = 1400,
    res = 150)


    # Boxplot
    output$boxplot <- renderPlot({
      df <- data.frame(
        Expression = as.vector(mat),
        Sample = rep(colnames(mat), each = nrow(mat)),
        Group = rep(colData(se)[[input$box_group_var]], each = nrow(mat))
      )
      ggplot(df, aes(x = Group, y = Expression)) +
        geom_boxplot() +
        labs(title = "Overall Gene Expression", x = input$box_group_var, y = "Expression") +
        ggpt
    }, res = 96)

    # Gene plots

    gene_plot_reactive <- reactive({
      req(input$gene)
      genes <- input$gene
      group_var <- input$geneplot_group
      color_by <- if (!is.null(input$geneplot_color_by)) input$geneplot_color_by else group_var
      show_labels <- input$geneplot_labels

      transform_axis_default <- is_rnaseq

      plots <- lapply(genes, function(gene) {
        transform_axis <- transform_axis_default

        if (is_rnaseq) {
          if (is.null(dds_obj)) stop("dds object not available for RNA-seq dataset.")
          dds_tmp <- estimateSizeFactors(dds_obj)
          norm_mat <- counts(dds_tmp, normalized = TRUE)
          if (!(gene %in% rownames(norm_mat))) stop(sprintf("Gene '%s' not found in normalized count matrix.", gene))
          df <- data.frame(exp_value = as.numeric(norm_mat[gene, ]), sample_id = colnames(norm_mat))
          if (transform_axis) {
            df$exp_value_plot <- df$exp_value + 1
            y_label <- "Normalized counts"
          } else {
            df$exp_value_plot <- df$exp_value
            y_label <- "Normalized counts"
          }
        } else {
          if ("exprs" %in% names(assays(se))) {
            exprs_mat <- tryCatch(as.matrix(assays(se)$exprs), error = function(e) as.matrix(assay(se)))
          } else {
            exprs_mat <- as.matrix(assay(se))
          }
          if (!(gene %in% rownames(exprs_mat))) stop(sprintf("Gene '%s' not found in expression matrix.", gene))
          df <- data.frame(exp_value = as.numeric(exprs_mat[gene, ]), sample_id = colnames(exprs_mat))
          df$exp_value_plot <- df$exp_value
          y_label <- "Expression"
          transform_axis <- FALSE
        }

        if (!is.null(group_var) && group_var %in% colnames(colData(se))) {
          df$Group <- factor(colData(se)[[group_var]])
        } else {
          df$Group <- factor(rep(NA, nrow(df)))
        }
        if (!is.null(color_by) && color_by %in% colnames(colData(se))) {
          df$color_var <- colData(se)[[color_by]]
        } else {
          df$color_var <- df$Group
        }

        p <- ggplot(df, aes(x = factor(Group), y = exp_value_plot, color = color_var)) +
          geom_jitter(position = position_jitter(width = 0.2), size = 3) +
          xlab("Experimental group") +
          labs(title = gene)

        if (show_labels) {
          p <- p + ggrepel::geom_text_repel(aes(label = sample_id), min.segment.length = 0)
        }

        color_vec <- df$color_var
        if (is.numeric(color_vec)) {
          p <- p + scale_color_gradient(name = color_by, low = "blue", high = "red")
        } else {
          p <- p + scale_color_brewer(name = color_by, palette = "Set1")
        }

        if (transform_axis && is_rnaseq) {
          p <- p + scale_y_log10(name = paste0(y_label, " (log10 scale)"))
        } else {
          p <- p + scale_y_continuous(name = y_label)
        }

        p + ggpt
      })

      plots
    })

    # controls for gene plots layout & color_by selection
    output$gene_plot_controls <- renderUI({
      req(input$gene)
      tagList(
        fluidRow(
          column(6, helpText("Gene plots are arranged side-by-side. Choose columns per row:")),
          column(6, selectInput("gene_ncol", "Columns per row", choices = c(1,2,3), selected = min(2, max(1, length(input$gene)))))
        ),
        fluidRow(
          column(6,
                 selectInput(
                   "geneplot_color_by",
                   "Color by:",
                   choices = colnames(colData(se)),
                   selected = if("Dose" %in% colnames(colData(se))) "Dose" else colnames(colData(se))[1]
                 )
          ),
          column(6,
                 helpText(paste("X-axis grouping (Experimental group):", input$geneplot_group))
          )
        )
      )
    })

    # combine gene plots side-by-side with dynamic ncol; compute height from rows
    output$geneplots_combined <- renderPlot({
      plots <- gene_plot_reactive()
      if (length(plots) == 0) return(invisible())

      ncol_choice <- as.integer(if (!is.null(input$gene_ncol)) input$gene_ncol else min(3, length(plots)))
      ncol_choice <- max(1, min(3, ncol_choice))
      combined <- wrap_plots(plots, ncol = min(ncol_choice, length(plots)))
      print(combined)
    }, height = function() {
      n_genes <- max(1, length(input$gene))
      ncol_choice <- as.integer(if (!is.null(input$gene_ncol)) input$gene_ncol else min(3, n_genes))
      ncol_choice <- max(1, min(3, ncol_choice))
      nrows <- ceiling(n_genes / ncol_choice)
      per_panel_height <- if (is_rnaseq) 420L else 380L
      total <- per_panel_height * nrows
      total <- min(total, 3000L)
      max(400L, total)
    }, res = 96)

    # GeneCards links
    output$genecard_links <- renderUI({
      req(input$gene)
      genes <- input$gene
      # Use gene symbols directly for GeneCards URLs
      gene_ids <- genes
      links <- lapply(gene_ids, function(g) {
        tags$a(href = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", g),
               target = "_blank", g)
      })
      tagList(tags$h5("GeneCards links:"), tags$ul(lapply(links, tags$li)))
    })

    # Download gene plot
    output$download_geneplot <- downloadHandler(
      filename = function() paste0(paste(input$gene, collapse = "_"), "_geneplot.png"),
      content = function(file) {
        plots <- gene_plot_reactive()
        if (length(plots) == 0) return(NULL)
        ncol_choice <- as.integer(if (!is.null(input$gene_ncol)) input$gene_ncol else min(3, length(plots)))
        ncol_choice <- max(1, min(3, ncol_choice))
        nrows <- ceiling(length(plots) / ncol_choice)
        width_in <- 4 * ncol_choice
        height_in <- 4 * nrows
        ggsave(file, wrap_plots(plots, ncol = ncol_choice), width = width_in, height = height_in, dpi = 300, limitsize = FALSE)
      }
    )

  }

  shinyApp(ui, server)
}
