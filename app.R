# app.R
options(shiny.maxRequestSize = 1024*1024^2)  # adjust if needed

library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(plotly)
library(DT)
library(processx)
library(blastula)
library(conflicted)
library(pheatmap)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(stringr)
library(forcats)


conflicted::conflicts_prefer(
  dplyr::intersect,
  dplyr::select,
  dplyr::filter,
  dplyr::arrange,
  dplyr::mutate,
  IRanges::union,
  dplyr::union,
  base::union
)

suppressPackageStartupMessages({
  if ("conflicted" %in% .packages(all.available = TRUE)) {
    conflicted::conflicts_prefer(
      dplyr::select,
      dplyr::filter,
      dplyr::arrange,
      dplyr::mutate,
      dplyr::intersect,
      base::setdiff
    )
  }
})

# ---------- Helpers ----------
chr_map_plink19 <- function(x){
  x <- toupper(as.character(x))
  x <- sub("^CHR", "", x)
  x[x=="X"] <- "23"; x[x=="Y"] <- "24"; x[x %in% c("MT","M","MTDNA")] <- "26"
  suppressWarnings(as.integer(x))
}
parse_p_robust <- function(p){
  if (is.numeric(p)) return(as.numeric(p))
  p_chr <- trimws(as.character(p))
  p_chr <- gsub("\\s+", "", p_chr)
  # if no dot but comma is present -> convert comma to dot (1,30E-06 -> 1.30E-06)
  needs_swap <- !grepl("\\.", p_chr) & grepl(",", p_chr)
  p_chr[needs_swap] <- gsub(",", ".", p_chr[needs_swap], fixed = TRUE)
  suppressWarnings(as.numeric(p_chr))
}

# safe OR operator (non-vectorized)
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) b else a
}

# always use base::intersect to avoid IRanges/dplyr conflicts
b_intersect <- function(x, y) base::intersect(x, y)
pick_first  <- function(x) if (length(x) > 0) x[[1]] else NULL
find_col_ci <- function(cols, patterns) {
  for (p in patterns) {
    hits <- grep(p, cols, ignore.case = TRUE, value = TRUE, perl = TRUE)
    if (length(hits)) return(hits[[1]])
  }
  NULL
}

ui <- navbarPage(
  title = div(
    style = "font-weight:700; font-size:22px; color:#1A4E8A;",
    HTML("🧬 NonSyn Variant Inspector")
  ),
  id = "topnav",
  
  # 👉 CSS personalitzat
  header = tags$style(HTML("
    /* TOTS els h3 del sidebarPanel */
    .sidebarPanel h3 {
      color: #1A4E8A;       /* blau */
      font-size: 20px;      /* mida */
      font-weight: 700;     /* gruix */
      margin-top: 20px;
      margin-bottom: 10px;
    }

    /* TOTS els h4 del sidebarPanel */
    .sidebarPanel h4 {
      color: #14507C;
      font-size: 17px;
      font-weight: 600;
      margin-top: 15px;
      margin-bottom: 8px;
    }
  ")),
  
  # ================= TAB 1: Analysis (Steps 1–4) =================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Analysis (Steps 1–4)</span>"),
    tags$style(HTML("
      .tab-content { padding-top: 0 !important; }
      .main-panel h3:first-child { margin-top: 0 !important; }
    ")),
    sidebarLayout(
      sidebarPanel(width = 4,
                   h3("Step 1 · Reading & selection",
                      style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"),
                   actionButton("info_00", "ℹ"), # Info about input file
                   fileInput("file", "Upload p-value table (TSV/CSV)", accept = c(".tsv",".txt",".csv")),
                   checkboxInput("header", "First row contains column names", TRUE),
                   radioButtons("sep", "Separator",
                                c("Tab \\t" = "\t", "Comma ," = ",", "Semicolon ;" = ";"),
                                selected = "\t"),
                   sliderInput("pthr", "-log10(P) threshold", min = 2, max = 20, value = 7.3, step = 0.1),
                   numericInput("flank", "Flank (+/- bp)", value = 10000, min = 0, step = 1000),
                   actionButton(
                     "build_ranges",
                     "➊ Generate selected_intervals",
                     style = "background-color: #ffdd57; color: black; font-weight: bold;"
                   ),
                   h4("Preview selected_intervals.range"),
                   verbatimTextOutput("ranges_preview"),
                   tags$hr(),
                   ####################
                   h3("Step 2 · Fill intervals with mapped SNP's from dbSNP catalog",
                      style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"),
                   helpText("Subsets the dbSNP VCF by intervals and creates a site-level VCF for dbNSFP."),
                   textInput("bcftools_path", "Path to bcftools (paste path)", value = ""),
                   verbatimTextOutput("bcftools_detected"),
                   textInput("sites_vcf", "Path to dbSNP file (.vcf.gz) (paste path)",
                             value = ""),
                   checkboxInput("sites_need_chr", "Force 'chr' prefix in intervals (if the VCF uses 'chr1')", FALSE),
                   checkboxInput("force_vcf40", "Force VCFv4.0 header (legacy compatibility)", TRUE),
                   actionButton(
                     "fill_intervals_sites",
                     "➋ Fill intervals (dbSNP + bcftools)",
                     style = "background-color: #ffdd57; color: black; font-weight: bold;"
                   ),

                   tags$hr(),
                   
                   h3("Step 3 · Annotate with dbNSFP",
                      style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"),
                   textInput("dbnsfp_dir",  "Path to dbNSFP folder (contains *.gz and search tool)(paste path)",
                             value = ""),
                   textInput("dbnsfp_tool", "Name of search tool (Java class), default: 'search_dbNSFP50a'", value = "search_dbNSFP50a"),
                   numericInput("java_gb",   "Java memory (GB)", value = 6, min = 2, step = 1),
                   checkboxInput("show_dbnsfp_warning", "Show duration warning", TRUE),
                   actionButton(
                     "run_dbnsfp",
                     "➌ Run dbNSFP",
                     style = "background-color: #ffdd57; color: black; font-weight: bold;"
                   ),
                 #  downloadButton("dl_ann",  "Download dbNSFP.out"),
                   tags$hr(),
                   
                   h3("Step 4 · Format dbNSFP output",
                      style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"),
                   actionButton(
                     "format_dbnsfp",
                     "➍ Format dbNSFP output",
                     style = "background-color: #ffdd57; color: black; font-weight: bold;"
                   ),
                   helpText("Converts numeric columns, deduplicates ';', and renames headers if matched."),
                   conditionalPanel(
                     condition = "output.norm_ready_flag",
                     tagList(
                       downloadButton("dl_dbnsfp_csv", "Download normalized CSV", style = "background-color: #ffdd57; color: black; font-weight: bold;"),
                       downloadButton("dl_dbnsfp_rds", "Download normalized RDS", style = "background-color: #ffdd57; color: black; font-weight: bold;")
                     )
                   ),
                   tags$hr(),
                   h4("Temporary folder"),
                   verbatimTextOutput("workdir_txt")
      ),
      mainPanel(
        h3("Manhattan plot"),
        plotlyOutput("manhattan", height = 520),
        tags$hr(),
        h4("SNPs ≥ threshold (select rows)"),
        DT::DTOutput("hits_tbl"),
        tags$hr(),
        h4("Extract VCF – bcftools logs"),
        verbatimTextOutput("sites_log"),
        tags$hr(),
        h4("Annotate with dbNSFP50a – logs"),
        verbatimTextOutput("dbnsfp_log")
      )
    )
  ),
  
  # ================= TAB 2: Visualization (Step 5) =================
  tabPanel(
    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>Visualization (Step 5)</span>"),
    sidebarLayout(
      sidebarPanel(width = 4,
                   h3("Step 5 · Showing results",
                      style = "color:#1A4E8A; font-size:22px; font-weight:700; margin-top:10px;"),
                   actionButton("info_01", "ℹ"), # Info about input file
                   radioButtons("viz_source", NULL,
                                choices = c("Upload external file" = "upload",
                                            "Use generated file in the latest session" = "generated"),
                                selected = "upload"),
                   
                   fileInput("viz_file", "Uploaded file",
                             accept = c(".csv"), buttonLabel = "Choose CSV",
                             placeholder = "No file selected"),
                   tags$hr(),
                   uiOutput("viz_controls")
      ),
      mainPanel(
        uiOutput("viz_tabs")  # tabs: Heatmap / GO / KEGG
      )
    )
  )
)

# ===========================================================================
# ---------- SERVER ----------
# ===========================================================================

server <- function(input, output, session){
  
  source("R/format_dbnsfp.R", local = TRUE)
  source("R/clinical_metrics.R", local = TRUE)
  
  suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(AnnotationDbi)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(plotly)
    library(dplyr)
  })
  
 # info files
  observeEvent(input$info_00, {
    showModal(modalDialog(
      title = "Info: Hits to plot",
      HTML('
<p style="text-align: justify;">
  Input file format [(*) needed columns]
</p>

<table style="border-collapse: collapse; width: 100%;" border="1">
  <thead>
    <tr>
      <th>CHR(*)</th>
      <th>SNP(*)</th>
      <th>BP(*)</th>
      <th>A1</th>
      <th>TEST</th>
      <th>NMISS</th>
      <th>OR</th>
      <th>STAT</th>
      <th>P(*)</th>
    </tr>
  </thead>

  <tbody>
    <tr><td>14</td><td>14:33967212:C:A</td><td>33967212</td><td>C</td><td>ADD</td><td>720</td><td>2.03</td><td>4.84</td><td>1.30E-06</td></tr>
    <tr><td>18</td><td>rs12955421</td><td>76002716</td><td>A</td><td>ADD</td><td>720</td><td>2.393</td><td>4.76</td><td>1.93E-06</td></tr>
    <tr><td>3</td><td>rs1145036</td><td>2328334</td><td>G</td><td>ADD</td><td>720</td><td>1.816</td><td>4.62</td><td>3.83E-06</td></tr>

    <tr><td>2</td><td>rs59874571</td><td>29998814</td><td>A</td><td>ADD</td><td>720</td><td>0.4769</td><td>-4.558</td><td>5.17E-06</td></tr>

    <tr><td>3</td><td>3:2328258:T:C</td><td>2328258</td><td>C</td><td>ADD</td><td>720</td><td>1.804</td><td>4.534</td><td>5.78E-06</td></tr>
    <tr><td>3</td><td>rs4263314</td><td>75290033</td><td>A</td><td>ADD</td><td>720</td><td>2.255</td><td>4.52</td><td>6.20E-06</td></tr>

    <tr><td>18</td><td>18:75996882:G:T</td><td>75996882</td><td>T</td><td>ADD</td><td>720</td><td>2.208</td><td>4.509</td><td>6.50E-06</td></tr>
    <tr><td>18</td><td>rs12969521</td><td>75996460</td><td>A</td><td>ADD</td><td>720</td><td>2.208</td><td>4.509</td><td>6.50E-06</td></tr>
    <tr><td>18</td><td>rs71361173</td><td>76000450</td><td>G</td><td>ADD</td><td>720</td><td>2.208</td><td>4.509</td><td>6.50E-06</td></tr>
    <tr><td>18</td><td>rs7229435</td><td>75996609</td><td>C</td><td>ADD</td><td>720</td><td>2.208</td><td>4.509</td><td>6.50E-06</td></tr>
    <tr><td>18</td><td>rs7243650</td><td>75996566</td><td>T</td><td>ADD</td><td>720</td><td>2.208</td><td>4.509</td><td>6.50E-06</td></tr>
    <tr><td>18</td><td>rs7244629</td><td>75997161</td><td>T</td><td>ADD</td><td>720</td><td>2.208</td><td>4.509</td><td>6.50E-06</td></tr>

    <tr><td>5</td><td>5:134520666:CT:C</td><td>134520666</td><td>C</td><td>ADD</td><td>720</td><td>0.5442</td><td>-4.499</td><td>6.81E-06</td></tr>
    <tr><td>10</td><td>rs4948491</td><td>61937130</td><td>G</td><td>ADD</td><td>720</td><td>0.5648</td><td>-4.492</td><td>7.05E-06</td></tr>
    <tr><td>3</td><td>rs1143978</td><td>2317302</td><td>C</td><td>ADD</td><td>720</td><td>1.775</td><td>4.464</td><td>8.05E-06</td></tr>
  </tbody>
</table>
'),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$info_01, {
    showModal(modalDialog(
      title = "Info: Format exterval file",
      p(style = "text-align: justify;",
        "Load formated dbNSFP output file as performet at step 4 (dbnsfp_normalized.csv)",
      ),
      easyClose = TRUE,
      footer = NULL
    ))
  })

############################################################### 
  
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(dplyr::select, dplyr::filter, dplyr::arrange,
                                 dplyr::mutate, base::setdiff, base::union, base::intersect)
  }
  
  # Safe FDR helper
  get_qcut <- function(x, default = 0.05) {
    if (!is.null(x) && length(x) == 1 && is.finite(x)) as.numeric(x) else default
  }
  plot_placeholder <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 5) +
      ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::theme_void()
  }
  
  
  # Per-session Shiny working folder (safe for downloads)
  workdir <- file.path(tempdir(), paste0("gwas_inspector_", Sys.getpid()))
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  session$onSessionEnded(function(...) {
    if (dir.exists(workdir)) unlink(workdir, recursive = TRUE, force = TRUE)
  })
  output$workdir_txt <- renderText(workdir)
  
  # ---------- Step 1: Reading & Manhattan ----------
  dat_raw <- reactive({
    req(input$file)
    read_delim(
      file   = input$file$datapath,
      delim  = input$sep,
      col_names = input$header,
      locale = locale(decimal_mark = ".", grouping_mark = ","),
      show_col_types = FALSE,
      na = c("", "NA")
    )
  })
  
  gwas_df <- reactive({
    df <- dat_raw()
    names(df) <- tolower(names(df))
    req(all(c("chr","bp","p") %in% names(df)))
    if (!"snp" %in% names(df)) df$snp <- NA_character_
    
    BP <- if (is.numeric(df$bp)) df$bp else readr::parse_number(as.character(df$bp))
    Pval <- parse_p_robust(df$p)
    CHR <- chr_map_plink19(df$chr)
    
    df <- df %>%
      mutate(
        CHR = CHR,
        BP = as.numeric(BP),
        Pval = as.numeric(Pval),
        logp = -log10(Pval)
      ) %>%
      filter(!is.na(CHR), !is.na(BP), !is.na(Pval), Pval > 0)
    
    if (!"a1" %in% names(df)) df$a1 <- NA_character_
    df
    
    df
  })
  
  output$manhattan <- renderPlotly({
    df <- gwas_df(); req(nrow(df) > 0)
    # Manhattan positions
    chr_lengths <- df %>% dplyr::group_by(CHR) %>%
      dplyr::summarise(chr_len = max(BP, na.rm=TRUE), .groups="drop") %>%
      dplyr::arrange(CHR) %>% dplyr::mutate(tot = cumsum(chr_len) - chr_len)
    
    dfp <- df %>% dplyr::inner_join(chr_lengths, by = "CHR") %>%
      dplyr::arrange(CHR, BP) %>% dplyr::mutate(pos = BP + tot, colgrp = CHR %% 2)
    
    axis_df <- dfp %>% dplyr::group_by(CHR) %>% dplyr::summarise(center = median(pos), .groups="drop")
    
    thr_y <- input$pthr
    g <- ggplot(dfp, aes(x=pos, y=logp,
                         text = paste0(
                           "CHR: ", CHR,
                           "<br>BP: ", BP,
                           "<br>P: ", signif(Pval,3),
                           ifelse(!is.na(snp), paste0("<br>SNP: ", snp), "")
                         ))) +
      geom_point(aes(color = factor(colgrp)), size = 0.6) +
      geom_hline(yintercept = thr_y, linetype = "dashed") +
      scale_color_manual(values = c("#377eb8", "#ff7f00")) +
      scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
      labs(x = "Chromosome", y = "-log10(P)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text")
  })
  
  # Hit table by threshold and row selection
  hits_df <- reactive({
    df <- gwas_df(); req(nrow(df)>0)
    df %>%
      dplyr::filter(logp >= input$pthr) %>%
      dplyr::arrange(desc(logp)) %>%
      dplyr::select(CHR, BP, snp, p = Pval, logp, any_of("a1"))
  })
  
  output$hits_tbl <- renderDT({
    datatable(hits_df(), selection = "multiple", rownames = FALSE, options = list(pageLength = 10))
  })
  
  # --- helper: PLINK chromosome number -> label (1..22, X, Y, MT)
  chr_label_plink <- function(chr_num) {
    out <- as.character(chr_num)
    out[out == "23"] <- "X"
    out[out == "24"] <- "Y"
    out[out == "26"] <- "MT"
    out
  }
  
  # ---------- Step 1 -> generate selected_intervals (4 columns) ----------
  selected_ranges_path <- reactiveVal(NULL)
  
  observeEvent(input$build_ranges, {
    h <- hits_df(); req(nrow(h) > 0)
    
    sel_idx <- input$hits_tbl_rows_selected
    picks <- if (length(sel_idx) > 0) {
      h[sel_idx, , drop = FALSE]
    } else {
      showNotification("You did not select any rows: I will use ALL SNPs ≥ threshold.", type = "warning")
      h
    }
    
    flank <- as.integer(if (!is.null(input$flank) && !is.na(input$flank)) input$flank else 0)
    
    # Safe columns
    CHR <- as.integer(picks$CHR)
    BP  <- as.integer(picks$BP)
    
    ranges <- data.frame(
      CHR_lab = chr_label_plink(CHR),
      BP1     = pmax(1L, BP - flank),
      BP2     = BP + flank,
      stringsAsFactors = FALSE
    )
    
    # LABEL: rsID if present; otherwise "C{CHR}:{BP}±{flank}"
    base_label <- ifelse(!is.na(picks$snp) & nzchar(picks$snp) & grepl("^rs", picks$snp, ignore.case = TRUE),
                         picks$snp,
                         paste0("C", ranges$CHR_lab, ":", BP))
    label <- if (flank > 0) paste0(base_label, "±", flank) else base_label
    
    # Order and deduplicate labels if repeated
    ord <- order(ranges$CHR_lab, ranges$BP1, ranges$BP2, na.last = NA)
    ranges <- ranges[ord, , drop = FALSE]
    label  <- label[ord]
    label  <- make.unique(label, sep = "_")
    
    # Remove rows with NA
    keep <- !is.na(ranges$CHR_lab) & !is.na(ranges$BP1) & !is.na(ranges$BP2)
    ranges <- ranges[keep, , drop = FALSE]
    label  <- label[keep]
    
    # Build lines with SPACES (no tabs) → "CHR BP1 BP2 LABEL"
    lines <- sprintf("%s %d %d %s", ranges$CHR_lab, ranges$BP1, ranges$BP2, label)
    lines <- unique(lines)
    lines <- lines[nchar(lines) > 0]
    
    # Write file
    rng_path <- file.path(workdir, "selected_intervals.range")
    writeLines(lines, rng_path, useBytes = TRUE)
    
    # Check first lines (≥4 tokens)
    head5 <- readLines(rng_path, n = 5, warn = FALSE)
    ok_tokens <- length(head5) == 0 || all(sapply(strsplit(head5, "\\s+"), function(x) length(x) >= 4))
    validate(need(ok_tokens, "The interval file does not match the CHR BP1 BP2 LABEL format."))
    
    selected_ranges_path(rng_path)
    
    # Preview in UI + console
    # Show compact preview in the UI (only first 10 lines) + total summary
    output$ranges_preview <- renderText({
      n_total <- length(lines)
      n_show  <- min(10, n_total)
      head10  <- utils::head(lines, n_show)
      
      paste0(
        "Total intervals: ", n_total, "\n",
        if (n_total > 10) "Showing the first 10:\n" else "Showing all:\n",
        paste(head10, collapse = "\n")
      )
    })
  })
  
  ## ======================
  ##   Step 2 ====== Extract dbSNP SNPs in the intervals
  ## ======================
  # -- chromosome normalizers --
  .norm_chr <- function(x){
    x <- as.character(x)
    x <- sub("^chr","", x, ignore.case = TRUE)
    x[x == "23"] <- "X"; x[x == "24"] <- "Y"
    # homogenize MT -> M (later to prefix as chrM)
    x[x %in% c("M","m","Mt","MT","MTDNA","mtdna")] <- "M"
    x
  }
  .add_chr <- function(x){
    x <- as.character(x)
    ifelse(grepl("^chr", x, ignore.case = TRUE), x, paste0("chr", x))
  }
  
  
  vcf_out_path <- reactiveVal(NULL)
  
  observeEvent(input$fill_intervals_sites, {
    req(selected_ranges_path())
    bcftools <- input$bcftools_path %||% "bcftools"
    sites_vcf <- "/Volumes/LaCie_Drive/ANNOTATE_RS_CODES_Resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    validate(need(file.exists(sites_vcf), paste("Does not exist:", sites_vcf)))
    
    rng_in  <- selected_ranges_path()
    workdir <- dirname(rng_in)
    
    # --- Robust reading of ranges (CHR START END [LABEL]) ---
    rng <- try(read.table(rng_in, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "", comment.char = ""),
               silent = TRUE)
    if (inherits(rng, "try-error") || ncol(rng) < 3) {
      rng <- read.table(rng_in, header = FALSE, sep = "",
                        stringsAsFactors = FALSE, quote = "", comment.char = "")
    }
    validate(need(ncol(rng) >= 3, "selected_intervals.range must have at least 3 columns (CHR START END)."))
    
    # Rename first 3 columns without touching extras
    nm <- names(rng); nm[1:3] <- c("CHR","START","END"); names(rng) <- nm
    
    rng$CHR   <- .norm_chr(rng$CHR)
    rng$START <- suppressWarnings(as.integer(rng$START))
    rng$END   <- suppressWarnings(as.integer(rng$END))
    validate(need(all(is.finite(rng$START)) && all(is.finite(rng$END)), "START/END must be integers."))
    
    # --- Merge overlaps by CHR ---
    merge_ranges <- function(df){
      df <- df[order(df$CHR, df$START, df$END), ]
      out <- list()
      cur <- df[1, , drop = FALSE]
      if (nrow(df) > 1) {
        for (i in 2:nrow(df)) {
          if (identical(df$CHR[i], cur$CHR) && df$START[i] <= cur$END) {
            cur$END <- max(cur$END, df$END[i])
          } else {
            out[[length(out) + 1]] <- cur
            cur <- df[i, , drop = FALSE]
          }
        }
      }
      out[[length(out) + 1]] <- cur
      do.call(rbind, out)
    }
    
    rng$CHR   <- .norm_chr(rng$CHR)
    rng_merged <- merge_ranges(rng[, c("CHR","START","END")])
    
    # --- Detect if VCF uses 'chr' prefix by reading header ---
    hdr <- tryCatch(system2(bcftools, c("view","-h", normalizePath(sites_vcf)),
                            stdout = TRUE, stderr = TRUE), error = function(e) character(0))
    vcf_uses_chr <- any(grepl("^##contig=<ID=chr", hdr))
    
    # Apply prefix only if VCF uses it
    if (isTRUE(vcf_uses_chr)) {
      rng_merged$CHR <- .add_chr(rng_merged$CHR)
    }
    
    # Write reduced regions file
    regions <- file.path(workdir, "selected_intervals.regions.txt")
    write.table(rng_merged, regions, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # --- Ensure tabix index (.tbi) for -R ---
    idx <- paste0(normalizePath(sites_vcf), ".tbi")
    if (!file.exists(idx)) {
      system2(bcftools, c("index","-t", normalizePath(sites_vcf)), stdout = TRUE, stderr = TRUE)
    }
    
    # --- Extract variants with bcftools -R ---
    out_vcf_plain <- file.path(workdir, "selected_sites.vcf")
    
    args_plain <- c("view", "-R", regions, "-Ov", "-o", out_vcf_plain, normalizePath(sites_vcf))
    out_plain  <- tryCatch(system2(bcftools, args = args_plain, stdout = TRUE, stderr = TRUE),
                           error = function(e) e)
    
    # ====== ROBUST VARIANT COUNT ======
    count_variants_fast <- function(vcf_path){
      # try data.table::fread (only column 1) skipping to #CHROM
      n <- tryCatch(
        nrow(data.table::fread(vcf_path, sep="\t", header=TRUE, skip="#CHROM",
                               select = 1, showProgress = FALSE)),
        error = function(e) NA_integer_
      )
      if (is.na(n)) {
        # fallback: count non-header lines by streaming
        con <- file(vcf_path, open = "r")
        on.exit(close(con), add = TRUE)
        tot <- 0L
        repeat {
          ln <- readLines(con, n = 50000, warn = FALSE)
          if (!length(ln)) break
          tot <- tot + sum(substr(ln, 1, 1) != "#")
        }
        n <- tot
      }
      as.integer(n)
    }
    
    snp_n <- if (file.exists(out_vcf_plain)) count_variants_fast(out_vcf_plain) else 0L
    
    # (optional) add bcftools stats summary to log
    bcf_stats <- character(0)
    if (file.exists(out_vcf_plain) && snp_n > 0) {
      bcf_stats <- tryCatch(system2(bcftools, c("stats", out_vcf_plain), stdout = TRUE),
                            error = function(e) character(0))
    }
    
    # Force VCFv4.0 header (dbNSFP-friendly), as you were already doing
    if (file.exists(out_vcf_plain)) {
      txt <- readLines(out_vcf_plain, warn = FALSE)
      i <- grep("^##fileformat=VCFv4\\.", txt)
      if (length(i)) {
        txt[i[1]] <- "##fileformat=VCFv4.0"
        writeLines(txt, out_vcf_plain)
      }
    }
    
    # Flag + notification
    if (file.exists(out_vcf_plain) && snp_n > 0L) {
      vcf_out_path(out_vcf_plain)
      showNotification(sprintf("✅ %s (ranges: %d → %d) ready for dbNSFP — SNPs extracted: %d",
                               basename(out_vcf_plain), nrow(rng), nrow(rng_merged), snp_n),
                       type = "message")
    } else {
      showNotification("❌ No variants were extracted. Check build (hg38) and 'chr' prefix.",
                       type = "error")
    }
    
    # Log (includes snp_n and, if present, stats)
    output$sites_log <- renderText({
      paste0(
        "Ranges IN : (", nrow(rng), ") ", rng_in, "\n",
        "Ranges OUT: (merged: ", nrow(rng_merged), ") ", regions, "\n",
        "Sites VCF: ", normalizePath(sites_vcf), "\n",
        "VCF uses 'chr': ", vcf_uses_chr, "\n",
        "Index (.tbi): ", if (file.exists(paste0(normalizePath(sites_vcf), ".tbi"))) "OK" else "NO", "\n\n",
        "SNPs extracted: ", snp_n, "\n",
        "CMD: ", bcftools, " ", paste(args_plain, collapse = " "), "\n",
        if (length(bcf_stats)) paste(c(bcf_stats), collapse = "\n") else "",
        "\nOutput: ", out_vcf_plain, " [", if (snp_n > 0) "with variants" else "empty", "]\n"
      )
    })
    
  })
  
 # observe({
 #   path <- Sys.which("bcftools")
 #   output$bcftools_detected <- renderText({
 #     paste("Detected path:", path)
 #   })
 # })
  
  ## ======================
  ##   Step 3 · Annotation with dbNSFP
  ## ======================
 
  # --- state/outputs ---
  ann_path        <- reactiveVal(NULL)   # path to generated .out
  ann_ready       <- reactiveVal(FALSE)  # flag for conditional UI
  dbnsfp_log_text <- reactiveVal("")     # log buffer
  
  # log renderer (in UI: verbatimTextOutput("dbnsfp_log"))
  output$dbnsfp_log <- renderText(dbnsfp_log_text())
  
  # helper to append lines to log
  append_dbnsfp_log <- function(..., sep = "\n") {
    old <- dbnsfp_log_text()
    dbnsfp_log_text(paste(c(old, paste(..., collapse = sep)), collapse = sep))
  }
  
  # ------------- VCF utilities -------------
  ensure_plain_vcf <- function(vcf_path, workdir) {
    # If it is .vcf.gz → decompress to .vcf
    if (grepl("\\.vcf\\.gz$", vcf_path, ignore.case = TRUE)) {
      vcf_plain <- file.path(workdir, "selected_intervals.vcf")
      if (requireNamespace("R.utils", quietly = TRUE)) {
        R.utils::gunzip(vcf_path, destname = vcf_plain, overwrite = TRUE, remove = FALSE)
      } else {
        zz <- system2("bgzip", c("-d","-c", shQuote(vcf_path)), stdout = TRUE, stderr = TRUE)
        writeLines(zz, vcf_plain, useBytes = TRUE)
      }
      return(vcf_plain)
    }
    vcf_path
  }
  
  force_vcf_version <- function(vcf_in, vcf_out, target = "VCFv4.0") {
    stopifnot(file.exists(vcf_in))
    lines <- readLines(vcf_in, warn = FALSE)
    validate(need(length(lines) > 0, "The VCF is empty."))
    # Ensure fileformat header is VCFv4.0 and reference is hg38
    if (grepl("^##fileformat=VCFv", lines[1])) {
      lines[1] <- sub("^##fileformat=VCFv[0-9.]+", paste0("##fileformat=", target), lines[1])
    } else {
      lines <- c(paste0("##fileformat=", target), lines)
    }
    if (!any(grepl("^##reference=", lines)))
      lines <- append(lines, "##reference=hg38", after = 1)
    writeLines(lines, vcf_out, useBytes = TRUE)
    vcf_out
  }
  
  # ----------- actual dbNSFP run -----------
  run_dbnsfp_real <- function() {
    req(vcf_out_path())
    vcf_in <- normalizePath(vcf_out_path(), mustWork = TRUE)
    
    db_dir  <- input$dbnsfp_dir
    tool_in <- input$dbnsfp_tool   # "search_dbNSFP50a" or .class path
    validate(need(dir.exists(db_dir), paste0("Folder does not exist: ", db_dir)))
    
    # Minimal check of database files
    chr_files <- list.files(db_dir, pattern="^dbNSFP5\\.0a_variant\\.chr(\\d+|[XYM])\\.gz$", full.names = TRUE)
    validate(need(length(chr_files) >= 24 && file.exists(file.path(db_dir, "dbNSFP5.0_gene.gz")),
                  "Some database files are missing (chr*.gz / dbNSFP5.0_gene.gz)."))
    
    java_mem <- as.integer(input$java_gb %||% 6)
    
    # Main class name (remove .class suffix if present)
    main_class <- basename(tool_in)
    main_class <- sub("\\.class$", "", main_class, ignore.case = TRUE)
    
    # → Prepare VCF v4.0
    vcf_plain <- ensure_plain_vcf(vcf_in, workdir)
    vcf_v40   <- file.path(workdir, "selected_intervals.v40.vcf")
    vcf_v40   <- force_vcf_version(vcf_plain, vcf_v40, "VCFv4.0")
    
    # count variants
    n_var <- sum(!grepl("^#", readLines(vcf_v40, warn = FALSE)))
    validate(need(n_var > 0, "The VCF (v4.0) does not contain variants to annotate."))
    
    out_file <- file.path(workdir, "dbnsfp.out")
    cmd <- list(
      cmd  = "java",
      args = c(paste0("-Xmx", java_mem, "g"),
               main_class,
               "-i", normalizePath(vcf_v40, mustWork = TRUE),
               "-o", normalizePath(out_file, mustWork = FALSE))
    )
    
    # Start log with header
    dbnsfp_log_text(paste0(
      "DB dir: ", db_dir, "\n",
      "Tool:   ", main_class, "  (running with wd = db_dir)\n",
      "VCF:    ", vcf_v40, "\n\n",
      "CMD:\n", cmd$cmd, " ", paste(cmd$args, collapse = " "), "\n\n"
    ))
    
    withProgress(message = "Running dbNSFP", value = 0.5, {
      res <- tryCatch({
        old <- getwd(); on.exit(setwd(old), add = TRUE)
        setwd(db_dir)  # run inside database folder
        system2(cmd$cmd, args = cmd$args, stdout = TRUE, stderr = TRUE, wait = TRUE)
      }, error = function(e) { attr(e, "status") <- 1L; e })
      
      status <- attr(res, "status"); if (is.null(status)) status <- 0L
      out_exists <- file.exists(out_file) && file.info(out_file)$size > 0
      
      # Append result and preview (if exists) to log
      append_dbnsfp_log(
        paste0("exit=", status, "\n"),
        paste(res, collapse = "\n"),
        if (out_exists) {
          prev <- paste(utils::head(readLines(out_file, warn = FALSE), 20), collapse = "\n")
          paste0("\n\n✔ Annotation completed.\nPreview (20 lines):\n", prev)
        } else {
          "\n✖ No output was generated."
        }
      )
      
      validate(need(out_exists, "dbNSFP did not generate any output; check the log above."))
      ann_path(out_file)
      ann_ready(TRUE)
      
      # (optional) add line with generated file size
      sz <- file.info(out_file)$size
      sz_txt <- format(sz, big.mark = ".", decimal.mark = ",")
      append_dbnsfp_log("", paste0("File ready: ", out_file, " (", sz_txt, " bytes)"))
      
      setProgress(1)
      showNotification(paste0("dbNSFP OK: ", basename(out_file)), type = "message")
    })
  }
  
  # --------- pre-run warning/confirmation ---------
  observeEvent(input$run_dbnsfp, {
    ann_ready(FALSE)
    ann_path(NULL)
    if (isTRUE(input$show_dbnsfp_warning)) {
      showModal(modalDialog(
        title = "This process may take time (minutes or hours)",
        tagList(
          p("dbNSFP annotation can be slow, especially from an external drive."),
          p("You can leave the app running; you will see the progress log here.")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_run_dbnsfp", "Understood, continue", class = "btn-primary")
        ),
        easyClose = FALSE
      ))
    } else {
      run_dbnsfp_real()
    }
  })
  
  observeEvent(input$confirm_run_dbnsfp, {
    removeModal()
    run_dbnsfp_real()
  })
  
  # --------- safe download of .out ---------
  output$dl_ann <- downloadHandler(
    filename = function() paste0("dbnsfp_annotation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".out"),
    contentType = "text/plain",
    content = function(file) {
      p <- ann_path()
      validate(need(!is.null(p) && file.exists(p), "There is no annotation available to download."))
      ok <- file.copy(p, file, overwrite = TRUE)
      validate(need(isTRUE(ok), "Could not copy the annotation file."))
    }
  )
  
  # (Optional) flag to condition UI
  output$ann_ready_flag <- reactive({ isTRUE(ann_ready()) })
  outputOptions(output, "ann_ready_flag", suspendWhenHidden = FALSE)
  
  
  #  ================================================================
  
  #  ============= Step 4
  dbnsfp_norm_path_csv <- reactiveVal(NULL)
  dbnsfp_norm_path_rds <- reactiveVal(NULL)
  
  observeEvent(input$format_dbnsfp, {
    req(ann_path())
    res <- tryCatch(
      normalize_dbnsfp(in_file = ann_path(), out_dir = workdir),
      error = function(e) { showNotification(e$message, type = "error"); NULL }
    )
    req(!is.null(res))
    
    dbnsfp_norm_path_csv(res$csv)
    dbnsfp_norm_path_rds(res$rds)
    
    output$dbnsfp_norm_tbl <- DT::renderDT({
      DT::datatable(res$preview, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    showNotification(
      paste0("dbNSFP normalized: ", basename(res$csv), " (", res$n_rows, " x ", res$n_cols, ")"),
      type = "message"
    )
  })
  
  # --- Flag to show download buttons only when files are ready
  output$norm_ready_flag <- reactive({
    p <- dbnsfp_norm_path_csv()
    !is.null(p) && file.exists(p)
  })
  outputOptions(output, "norm_ready_flag", suspendWhenHidden = FALSE)
  
  # --- Download normalized CSV
  output$dl_dbnsfp_csv <- downloadHandler(
    filename = function() paste0("dbnsfp_normalized_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    contentType = "text/csv",
    content = function(file) {
      p <- dbnsfp_norm_path_csv()
      validate(need(!is.null(p) && file.exists(p), "There is no normalized CSV to download."))
      ok <- file.copy(p, file, overwrite = TRUE)
      validate(need(isTRUE(ok), "Could not copy the normalized CSV."))
    }
  )
  
  # 4.0) Normalized dataset for Visualization (independent)
  dbnsfp_norm_df <- reactive({
    src <- input$viz_source %||% "upload"
    
    # a) Upload file mode
    if (identical(src, "upload")) {
      req(input$viz_file)
      f <- input$viz_file$datapath
      validate(need(file.exists(f), "The uploaded file does not exist."))
      df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
      validate(need(!is.null(df) && nrow(df) > 0, "Could not read the uploaded CSV."))
      attr(df, "source") <- "upload"
      return(df)
    }
    
    # b) Use the file generated by the app (Step 3.5)
    if (identical(src, "generated")) {
      path_csv <- dbnsfp_norm_path_csv()
      validate(need(!is.null(path_csv) && file.exists(path_csv),
                    "There is no dbnsfp_normalized.csv generated in this session."))
      df <- tryCatch(read.csv(path_csv, check.names = FALSE), error = function(e) NULL)
      validate(need(!is.null(df) && nrow(df) > 0, "The generated CSV is empty or unreadable."))
      attr(df, "source") <- "generated"
      return(df)
    }
    
    # fallback
    validate(need(FALSE, "Select a data source for visualization."))
  })
  
  observeEvent(dbnsfp_norm_df(), {
    df <- dbnsfp_norm_df(); req(is.data.frame(df))
    nms <- names(df)
    
    # common candidates
    cand <- c("chr","CHR","Chromosome","chromosome","chrom","CHROM","#CHROM")
    cand <- cand[cand %in% nms]
    
    # if there are several, prioritize 'chr'/'CHR' and then variants
    if (length(cand)) {
      ord <- c("chr","CHR","Chromosome","chromosome","chrom","CHROM","#CHROM")
      cand <- cand[order(match(cand, ord))]
      session$userData$chr_col <- cand[1]
    } else {
      session$userData$chr_col <- NULL
    }
  }, ignoreInit = FALSE)
  
  
  # === Visualization base: apply chr filter and guarantee rsid_final ===
  viz_base_df <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    snp_col <- session$userData$snp_col %||% "SNP"
    chr_col <- session$userData$chr_col
    bp_col  <- session$userData$bp_col
    
    # chr normalizer: "chr1" -> "1"
    norm_chr <- function(x) sub("^chr", "", tolower(trimws(as.character(x))))
    
    # chromosome filter (if applicable)
    if (!is.null(chr_col) && chr_col %in% names(df) &&
        !is.null(input$hm_chr) && input$hm_chr != "all") {
      key <- norm_chr(input$hm_chr)
      df  <- df[norm_chr(df[[chr_col]]) == key, , drop = FALSE]
    }
    
    # build rsid_final column (100% present)
    rs <- if (!is.null(snp_col) && snp_col %in% names(df)) as.character(df[[snp_col]]) else NA_character_
    empty <- is.na(rs) | rs == ""
    if (any(empty)) {
      if (!is.null(chr_col) && chr_col %in% names(df) && !is.null(bp_col) && bp_col %in% names(df)) {
        rs[empty] <- paste0(df[[chr_col]][empty], ":", df[[bp_col]][empty])
      } else {
        rs[empty] <- paste0("var_", which(empty))
      }
    }
    df$rsid_final <- rs
    
    # annotate which is the effective ID column for the rest of modules
    attr(df, "snp_col_effective") <- "rsid_final"
    attr(df, "chr_col") <- chr_col
    attr(df, "bp_col")  <- bp_col
    df
  })
  
  # when normalized data becomes available for visualization
  observeEvent(dbnsfp_norm_df(), {
    df <- dbnsfp_norm_df(); req(nrow(df) > 0)
    
    # use the chromosome column detected previously
    chr_col <- session$userData$chr_col
    if (is.null(chr_col) || !(chr_col %in% names(df))) {
      updateSelectInput(session, "hm_chr", choices = c("All" = "all"), selected = "all")
      return()
    }
    
    # normalize labels like "chr1" -> "1"
    norm_chr <- function(x) sub("^chr", "", tolower(trimws(as.character(x))))
    vals <- unique(df[[chr_col]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    vals_n <- norm_chr(vals)
    
    # show labels as "chrN" but keep normalized value
    # separate numeric and special for better ordering
    is_num <- suppressWarnings(!is.na(as.integer(vals_n)))
    nums   <- sort(unique(as.integer(vals_n[is_num])))
    specs  <- setdiff(unique(toupper(vals_n[!is_num])), character(0))
    # prioritize X,Y,MT if present
    spec_order <- c("X","Y","MT","M","MTDNA")
    specs <- c(intersect(spec_order, specs), setdiff(specs, spec_order))
    
    choices <- character(0)
    if (length(nums))  choices <- c(choices, stats::setNames(as.character(nums), paste0("chr", nums)))
    if (length(specs)) choices <- c(choices, stats::setNames(specs, paste0("chr", specs)))
    
    if (!length(choices)) {
      updateSelectInput(session, "hm_chr", choices = c("All"="all"), selected="all")
    } else {
      updateSelectInput(session, "hm_chr",
                        choices = c("All" = "all", choices),
                        selected = "all")
    }
  }, ignoreInit = FALSE)
  
  # --- Download normalized RDS (optional)
  output$dl_dbnsfp_rds <- downloadHandler(
    filename = function() paste0("dbnsfp_normalized_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    contentType = "application/octet-stream",
    content = function(file) {
      p <- dbnsfp_norm_path_rds()
      validate(need(!is.null(p) && file.exists(p), "There is no normalized RDS to download."))
      ok <- file.copy(p, file, overwrite = TRUE)
      validate(need(isTRUE(ok), "Could not copy the normalized RDS."))
    }
  )
  
  # ===========================================================================
  #  Step 5 · Visualization
  # ===========================================================================
  
  # ==========================
  # ===== Heatmap data (equivalent to score_input) most rellevant scores/ranks=====
  # ==========================
  
  heatmap_df <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    # key columns (using your normalized fields)
    snp_col <- session$userData$snp_col %||% "SNP"
    chr_col <- session$userData$chr_col %||% "chr"
    bp_col  <- session$userData$bp_col  %||% "BP"
    
    # default candidates
    default_scores <- c(
      "SIFT_score","SIFT4G_score",
      "Polyphen2_HDIV_score","Polyphen2_HVAR_score",
      "MutationTaster_score","MutationAssessor_score",
      "PROVEAN_score","REVEL_score",
      "GERP_91_mammals","phyloP17way_primate","phastCons17way_primate"
    )
    default_ranks <- c(
      "SIFT_converted_rankscore","SIFT4G_converted_rankscore",
      "Polyphen2_HDIV_rankscore","Polyphen2_HVAR_rankscore",
      "MutationAssessor_rankscore","PROVEAN_converted_rankscore",
      "REVEL_rankscore","GERP_91_mammals_rankscore",
      "phyloP17way_primate_rankscore","phastCons17way_primate_rankscore"
    )
    
    # UI selection
    if (identical(input$hm_mode, "Rankscore")) {
      cols_metrics <- input$hm_ranks
      if (length(cols_metrics) < 2)
        cols_metrics <- base::intersect(default_ranks, names(df))
    } else {
      cols_metrics <- input$hm_scores
      if (length(cols_metrics) < 2)
        cols_metrics <- base::intersect(default_scores, names(df))
    }
    
    validate(need(length(cols_metrics) >= 2,
                  "There are not at least 2 metric columns available for the heatmap."))
    
    # base columns
    cols_base <- base::intersect(
      c(snp_col, chr_col, bp_col, "ref","alt","genename"),
      names(df)
    )
    
    # --------------------------------------------------------
    # 1) df_part = només per fer el heatmap
    # --------------------------------------------------------
    df_part <- dplyr::bind_cols(
      df[, cols_base, drop = FALSE],
      df[, cols_metrics, drop = FALSE]
    )
    
    # --------------------------------------------------------
    # 2) Recuperar TOTES les columnes clíniques i anotacions
    # --------------------------------------------------------
    full <- dbnsfp_norm_df()
    validate(need(snp_col %in% names(full),
                  "Missing SNP column in the annotated dataset"))
    
    # left_join per SNP = annexa totes les dades disponibles
    out <- dplyr::left_join(df_part, full, by = snp_col)
    
    # --------------------------------------------------------
    # 3) Filters (chromosome)
    # --------------------------------------------------------
    if (!is.null(input$hm_chr) &&
        input$hm_chr != "all" &&
        chr_col %in% names(out)) {
      
      key <- suppressWarnings(as.numeric(input$hm_chr))
      key <- if (!is.na(key)) key else input$hm_chr
      
      out <- out[out[[chr_col]] == key, , drop = FALSE]
    }
    
    # --------------------------------------------------------
    # 4) Order + top N
    # --------------------------------------------------------
    if (chr_col %in% names(out) && bp_col %in% names(out)) {
      out <- out %>% dplyr::arrange(.data[[chr_col]], .data[[bp_col]])
    }
    
    n_show <- input$hm_topn %||% 100
    out <- utils::head(out, n_show)
    
    # --------------------------------------------------------
    # 5) Rename metric columns for heatmap aesthetics
    # --------------------------------------------------------
    new_names <- names(out)
    met_idx <- which(names(out) %in% cols_metrics)
    new_names[met_idx] <- gsub("_", "\n", new_names[met_idx])
    names(out) <- new_names
    
    # attributes needed later
    attr(out, "snp_col") <- snp_col
    attr(out, "chr_col") <- chr_col
    attr(out, "bp_col")  <- bp_col
    
    out
  })
  
 
  ################
  # Debug of detected schema (optional)
  output$dbnsfp_schema <- renderText({ "" })
  
  # Column detection (identical logic, but reading from dbnsfp_norm_df())
  observe({
    req(dbnsfp_norm_df())
    df   <- dbnsfp_norm_df()
    cols <- colnames(df); cols_lc <- tolower(cols)
    
    # SNP / rsid
    snp_lc <- pick_first(b_intersect(cols_lc, c("snp","rsid","rs_id","variant","variant_id")))
    snp_col <- if (!is.null(snp_lc)) cols[ match(snp_lc, cols_lc) ] else
      find_col_ci(cols, c("^snp$", "^rs.?id$", "^variant(_id)?$"))
    
    # BP / POS
    bp_lc <- pick_first(b_intersect(cols_lc, c("bp","pos","position","pos_1_based","hg38_pos_1_based","hg19_pos_1_based")))
    bp_col <- if (!is.null(bp_lc)) cols[ match(bp_lc, cols_lc) ] else
      find_col_ci(cols, c("^bp$", "^pos(ition)?$", "^pos_?1?_?based$", "^hg38_?pos", "^hg19_?pos"))
    
    # CHR
    chr_lc <- pick_first(b_intersect(cols_lc, c("chr","chromosome","#chr","hg38_chr","hg19_chr")))
    chr_col <- if (!is.null(chr_lc)) cols[ match(chr_lc, cols_lc) ] else
      find_col_ci(cols, c("^#?chr$", "^chrom(osome)?$", "^hg38_?chr$", "^hg19_?chr$"))
    
    # GENENAME
    gene_lc <- pick_first(b_intersect(cols_lc, c("genename","gene","symbol")))
    gene_col <- if (!is.null(gene_lc)) cols[ match(gene_lc, cols_lc) ] else
      find_col_ci(cols, c("^gene(name)?$", "^symbol$", "^hgnc"))
    
    # Metrics
    score_cols <- cols[ grepl("(_score$|^cadd_?phred$|^gerp(\\+\\+)?_rs$|^phylop[0-9]+way|^phastcons[0-9]+way)", cols_lc, perl=TRUE) ]
    rank_cols  <- cols[ grepl("rankscore$", cols_lc) ]
    pred_cols  <- cols[ grepl("_pred$", cols_lc) ]
    
    # GO / KEGG
    go_cols   <- b_intersect(cols, c("GO_biological_process","GO_cellular_component","GO_molecular_function"))
    if (length(go_cols) == 0)
      go_cols <- cols[ grep("^go_.*(process|function|component)", cols_lc) ]
    
    kegg_cols <- b_intersect(cols, c("Pathway_KEGG_id","Pathway_KEGG_full","Pathway_KEGG"))
    if (length(kegg_cols) == 0)
      kegg_cols <- cols[ grep("kegg", cols_lc) ]
    
    # store
    session$userData$snp_col    <- snp_col
    session$userData$bp_col     <- bp_col
    session$userData$chr_col    <- chr_col
    session$userData$gene_col   <- gene_col
    session$userData$pred_cols  <- sort(unique(pred_cols))
    session$userData$go_cols    <- sort(unique(go_cols))
    session$userData$kegg_cols  <- sort(unique(kegg_cols))
   
    session$userData$score_cols <- unique(score_cols)
    session$userData$rank_cols  <- unique(rank_cols)
    
    output$dbnsfp_schema <- renderText(paste(
      sprintf("Source: %s", attr(df, "source") %||% "unknown"),
      sprintf("Rows: %s, Columns: %s", nrow(df), ncol(df)),
      paste("SNP:",  ifelse(is.null(snp_col),"—", snp_col)),
      paste("CHR:",  ifelse(is.null(chr_col),"—", chr_col)),
      paste("BP :",  ifelse(is.null(bp_col) ,"—", bp_col)),
      paste("GENE:", ifelse(is.null(gene_col),"—", gene_col)),
      sep = "\n"
    ))
  })
  
  # small internal helper
  .norm_chr_vals <- function(v) {
    v <- as.character(v)
    v <- trimws(v)
    v <- sub("^chr", "", v, ignore.case = TRUE)
    v[v %in% c("23")] <- "X"
    v[v %in% c("24")] <- "Y"
    v[v %in% c("M","m","Mt","MTDNA","mtDNA")] <- "MT"
    v
  }
  
  # order 1..22, X, Y, MT + others
  .order_chr_vals <- function(u) {
    wanted <- c(as.character(1:22), "X", "Y", "MT")
    unique(c(intersect(wanted, u), setdiff(sort(u), wanted)))
  }
  
  output$viz_controls <- renderUI({
    req(dbnsfp_norm_df())
    
    tagList(
      # ================================================================
      # 1) HEATMAP CONTROLS
      # ================================================================
      conditionalPanel(
        condition = "input.viz_tabs == 'dbNSFP Heatmap'",
        
        h4("Heatmap (most rellevants scores/scoreranks)"),
        
        radioButtons(
          "scoreorrank", "Type",
          choices  = c("Score", "Rankscore"),
          selected = "Score",
          inline   = TRUE
        ),
        uiOutput("hm_chr_selector"),
        numericInput(
          "hm_topn", "Max rows (SNPs) to show",
          value = 50, min = 10, max = 5000, step = 5
        ),
        checkboxInput("hm_scale_cols", "Scale 0–1 per column", TRUE),
        downloadButton("dl_heatmap_png", "Download heatmap (PNG)")
      ),
      # ================================================================
      # 1) TABLE CONTROLS
      # ================================================================
      
      conditionalPanel(
        condition = "input.viz_tabs == 'Table dbNSFP'",
        h4("Full table dbNSFP scores/scoreranks"),
        uiOutput("hm_chr_selector2"),
      ),
      # ================================================================
      # 2) MANHATTAN & GENES CONTROLS
      # ================================================================
      conditionalPanel(
        condition = "input.viz_tabs == 'Manhattan & genes'",
        
        h4("Manhattan & genes"),
        
        radioButtons(
          "scoreorrank", "Type",
          choices  = c("Score", "Rankscore"),
          selected = "Score",
          inline   = TRUE
        ),
        
        h4("Metrics"),
        uiOutput("ms_metric_picker")
      ),
      
      # ================================================================
      # 3) GO CONTROLS
      # ================================================================
      conditionalPanel(
        condition = "input.viz_tabs == 'GO (BP/CC/MF)'",
        
        h4("GO"),
        sliderInput(
          "go_kegg_pcut", "FDR threshold (BH)",
          min = 0, max = 0.25, value = 0.05, step = 0.005
        ),
        downloadButton("dl_go_png", "Download GO (PNG)")
        
      ),
      
      # ================================================================
      # 4) KEGG CONTROLS
      # ================================================================
      conditionalPanel(
        condition = "input.viz_tabs == 'KEGG'",
        
        h4("KEGG"),
        sliderInput(
          "go_kegg_pcut", "FDR threshold (BH)",
          min = 0, max = 0.25, value = 0.05, step = 0.005
        ),
        numericInput(
          "kegg_top", "Top KEGG to show",
          value = 15, min = 5, max = 50, step = 1
        ),
        downloadButton("dl_kegg_png", "Download KEGG (PNG)")
      )
    )
  })
  
  # ================================================================
  # =============== CHROMOSOME SELECTOR (1) ================
  # ================================================================
  output$hm_chr_selector <- renderUI({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    # 1) Detectar columna CHR si no está registrada
    chr_col <- session$userData$chr_col
    if (is.null(chr_col) || !chr_col %in% names(df)) {
      cand <- grep("^#?chr$|^#?chrom$|chromosome",
                   names(df), ignore.case = TRUE, value = TRUE)
      chr_col <- cand[1]
    }
    
    # Si no existe → mensaje informativo
    if (is.null(chr_col) || !chr_col %in% names(df)) {
      return(tags$em("No chromosome column found; filter disabled."))
    }
    
    # 2) Normalizar cromosomas
    norm_chr <- function(x) {
      x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
      x[x == "23"] <- "X"
      x[x == "24"] <- "Y"
      x[x %in% c("M","MT","m","Mt","MTDNA","mtdna")] <- "MT"
      toupper(x)
    }
    
    chr_vals <- norm_chr(df[[chr_col]])
    chr_vals <- chr_vals[nzchar(chr_vals)]  # quitar vacíos
    chr_unique <- unique(chr_vals)
    
    # 3) Orden correcto "humano"
    order_chr <- function(v) {
      standard <- c(as.character(1:22), "X", "Y", "MT")
      v[match(v, standard, nomatch = NA)] <- standard[na.omit(match(v, standard, nomatch = NA))]
      intersect(standard, v)
    }
    
    chr_ordered <- order_chr(chr_unique)
    
    if (!length(chr_ordered)) {
      return(tags$em("Chromosome column has no valid values."))
    }
    
    # 4) Construir selector con etiquetas "chr1", "chr2", ...
    choices <- setNames(chr_ordered, paste0("chr", chr_ordered))
    
    # Mantener selección previa si existe
    selected_val <- isolate(input$hm_chr)
    if (is.null(selected_val) || !(selected_val %in% c("all", chr_ordered))) {
      selected_val <- "all"
    }
    
    # 5) Output final
    tagList(
      selectInput(
        "hm_chr",
        "Chromosome (optional filter)",
        choices = c("All" = "all", choices),
        selected = selected_val
      )
    )
  })
 
  # ================================================================
  # =============== CHROMOSOME SELECTOR (2) ================
  # ================================================================
  
  output$hm_chr_selector2 <- renderUI({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    # 1) Detectar columna CHR si no está registrada
    chr_col <- session$userData$chr_col
    if (is.null(chr_col) || !chr_col %in% names(df)) {
      cand <- grep("^#?chr$|^#?chrom$|chromosome",
                   names(df), ignore.case = TRUE, value = TRUE)
      chr_col <- cand[1]
    }
    
    # Si no existe → mensaje informativo
    if (is.null(chr_col) || !chr_col %in% names(df)) {
      return(tags$em("No chromosome column found; filter disabled."))
    }
    
    # 2) Normalizar cromosomas
    norm_chr <- function(x) {
      x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
      x[x == "23"] <- "X"
      x[x == "24"] <- "Y"
      x[x %in% c("M","MT","m","Mt","MTDNA","mtdna")] <- "MT"
      toupper(x)
    }
    
    chr_vals <- norm_chr(df[[chr_col]])
    chr_vals <- chr_vals[nzchar(chr_vals)]  # quitar vacíos
    chr_unique <- unique(chr_vals)
    
    # 3) Orden correcto "humano"
    order_chr <- function(v) {
      standard <- c(as.character(1:22), "X", "Y", "MT")
      v[match(v, standard, nomatch = NA)] <- standard[na.omit(match(v, standard, nomatch = NA))]
      intersect(standard, v)
    }
    
    chr_ordered <- order_chr(chr_unique)
    
    if (!length(chr_ordered)) {
      return(tags$em("Chromosome column has no valid values."))
    }
    
    # 4) Construir selector con etiquetas "chr1", "chr2", ...
    choices <- setNames(chr_ordered, paste0("chr", chr_ordered))
    
    # Mantener selección previa si existe
    selected_val <- isolate(input$hm_chr)
    if (is.null(selected_val) || !(selected_val %in% c("all", chr_ordered))) {
      selected_val <- "all"
    }
    
    # 5) Output final
    tagList(
      selectInput(
        "hm_chr",
        "Chromosome (optional filter)",
        choices = c("All" = "all", choices),
        selected = selected_val
      )
    )
  })
  
  # ================================================================
  # ================== VISUALIZATION STEP 5 ====================
  # ================================================================
  ########################################################
  # HEATMAP (STEP 5.1)
  ########################################################
   
  # 5.1) Visualization tabs
  output$viz_tabs <- renderUI({
    req(dbnsfp_norm_df())
    tagList(
    #  h3("Visualization"),
    #  p(em(sprintf("CSV: %s", attr(dbnsfp_norm_df(), "source") %||% "—"))),
      tabsetPanel(id = "viz_tabs",
                  tabPanel(title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🔥 dbNSFP Heatmap</span>"),
                           value = "dbNSFP Heatmap",
                           uiOutput("hm_chr_selector"),
                           uiOutput("dynamic_heatmap_plot"),
                  ),
                  tabPanel(title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📄 Table dbNSFP</span>"),
                           value = "Table dbNSFP",
                           uiOutput("hm_chr_selector2"),
                           radioButtons(
                             "clinical_class",
                             "Functional class",
                             choices = c(
                               "Default",
                               "ClinVar",
                               "OMIM / Orphanet",
                               "Pathogenicity",
                               "Conservation",
                               "Functional impact",
                               "LoF / Haploinsufficiency",
                               "Gene damage"
                             ),
                             selected = "Default",
                             inline = TRUE
                           ),
                           DT::DTOutput("heatmap_table")
                  ),
                  tabPanel(
                    title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>📊 Manhattan & genes</span>"),
                    value = "Manhattan & genes",
                    plotlyOutput("manhattan_scores", height = 420),
                  #  textOutput("rg_pos_txt"),
                    h4("UCSC Genome Browser"),
                    selectInput(
                      "ucsc_mirror", "UCSC mirror",
                      choices = c(
                        "Europe (recommended)" = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks",
                        "USA"                  = "https://genome.ucsc.edu/cgi-bin/hgTracks",
                        "Asia"                 = "https://genome-asia.ucsc.edu/cgi-bin/hgTracks"
                      ),
                      selected = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks"
                    ),
                    uiOutput("ucsc_viewer")
                  ),
                  tabPanel(title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧠 GO (BP/CC/MF)</span>"),
                           value = "GO (BP/CC/MF)",
                           plotOutput("go_bar"),
                           tags$hr(),
                           tags$hr(),
                           h4("GO table"),
                           DT::DTOutput("go_table")),
                  tabPanel(title = HTML("<span style='font-size:16px; font-weight:600; color:#1A4E8A;'>🧪 KEGG pathways</span>"),
                           value = "KEGG",
                           plotOutput("kegg_bar"), 
                           tags$hr(), 
                           h4("KEGG table"), 
                           DT::DTOutput("kegg_table")),
      )
    )
  })
  
  last_heatmap <- reactiveVal(NULL)
  last_go <- reactiveVal(NULL)
  last_kegg <- reactiveVal(NULL)
  
  # ==== Heatmap source (equivalent to your score_input()) ====
  score_input <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    # base columns and identifiers (tolerant to names)
    snp_col <- session$userData$snp_col %||% "SNP"
    chr_col <- session$userData$chr_col %||% "chr"
    bp_col  <- session$userData$bp_col  %||% "BP"
    
    wanted <- c(
      snp_col, chr_col, bp_col, "ref", "alt", "genename",
      "SIFT_score", "SIFT_pred", "SIFT_converted_rankscore",
      "SIFT4G_score", "SIFT4G_pred", "SIFT4G_converted_rankscore",
      "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore", "Polyphen2_HDIV_pred",
      "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore", "Polyphen2_HVAR_pred",
      "MutationTaster_score", "MutationTaster_pred",
      "MutationAssessor_score", "MutationAssessor_rankscore", "MutationAssessor_pred",
      "PROVEAN_score", "PROVEAN_converted_rankscore", "PROVEAN_pred",
      "REVEL_score", "REVEL_rankscore", "GERP++_RS_rankscore",
      "GERP_91_mammals", "GERP_91_mammals_rankscore",
      "phyloP17way_primate", "phyloP17way_primate_rankscore",
      "phastCons17way_primate", "phastCons17way_primate_rankscore"
    )
    
    keep <- base::intersect(wanted, names(df))
    validate(need(length(keep) >= 5, "The CSV does not contain enough of the expected columns for the heatmap."))
    
    out <- df[, keep, drop = FALSE]
    
    # chromosome filter if applicable
    if (!is.null(input$hm_chr) && input$hm_chr != "all" && chr_col %in% names(out)) {
      key <- suppressWarnings(as.numeric(input$hm_chr))
      key <- if (!is.na(key)) key else input$hm_chr
      out <- out[out[[chr_col]] == key, , drop = FALSE]
    }
    
    # order and limit top-N
    if (all(c(chr_col, bp_col) %in% names(out))) {
      out <- dplyr::arrange(out, .data[[chr_col]], .data[[bp_col]])
    }
    n_show <- input$hm_topn %||% 100
    out <- utils::head(out, n_show)
    
    attr(out, "snp_col") <- snp_col
    attr(out, "chr_col") <- chr_col
    attr(out, "bp_col")  <- bp_col
    out
  })
  
  output$heatmap_plot <- renderPlot({
    req(score_input())
    df <- score_input()
    validate(need(nrow(df) > 0, "⚠️ score_input is NULL or empty"))
    
    # 1) effective columns (always prefer 'rsid_final' if it exists)
    snp_col <- if ("rsid_final" %in% names(df)) "rsid_final" else (attr(df, "snp_col") %||% "SNP")
    chr_col <- attr(df, "chr_col") %||% "chr"
    bp_col  <- attr(df, "bp_col")  %||% "BP"
    
    # 2) metrics according to mode
    nms <- names(df)
    if (identical(input$scoreorrank, "Score")) {
      metric_cols <- c(
        nms[grepl("_score$", nms, perl = TRUE)],
        base::intersect(c("GERP_91_mammals","phyloP17way_primate","phastCons17way_primate"), nms)
      )
      title_txt <- "Heatmap Scores"
    } else {
      metric_cols <- c(
        nms[grepl("_rankscore$", nms, perl = TRUE)],
        base::intersect("GERP++_RS_rankscore", nms)
      )
      title_txt <- "Heatmap Rank-scores"
    }
    
    # 3) identifier columns (explicitly include rsid_final if present)
    id_cols <- base::intersect(
      unique(c(snp_col, "rsid_final", chr_col, bp_col, "ref", "alt", "genename")),
      nms
    )
    
    # 4) avoid duplicates/collisions
    metric_cols <- setdiff(unique(metric_cols), id_cols)
    cols_to_take <- unique(c(id_cols, metric_cols))
    validate(need(length(metric_cols) >= 2, "Choose at least 2 columns for the heatmap."))
    
    # 5) build df for heatmap
    df_hm <- df[, cols_to_take, drop = FALSE]
    
    # 6) rename ONLY metric columns (never IDs)
    new_names <- names(df_hm)
    m_idx <- match(metric_cols, new_names)
    new_names[m_idx] <- gsub("_", "\n", new_names[m_idx], fixed = TRUE)
    names(df_hm) <- new_names
    metric_cols_renamed <- names(df_hm)[m_idx]
    
    # 7) convert metrics to numeric
    df_hm[metric_cols_renamed] <- lapply(
      df_hm[metric_cols_renamed],
      function(x) suppressWarnings(as.numeric(as.character(x)))
    )
    
    # 8) matrix and row labels (rsid_final whenever available)
    score_matrix <- as.matrix(df_hm[, metric_cols_renamed, drop = FALSE])
    
    rs_vec <- if ("rsid_final" %in% names(df_hm)) df_hm[["rsid_final"]] else {
      if (snp_col %in% names(df_hm)) df_hm[[snp_col]] else NA
    }
    if (is.null(rs_vec) || all(is.na(rs_vec) | rs_vec == "")) {
      if (all(c(chr_col, bp_col) %in% names(df_hm))) {
        rs_vec <- paste0(df_hm[[chr_col]], ":", df_hm[[bp_col]])
      } else {
        rs_vec <- paste0("var_", seq_len(nrow(df_hm)))
      }
    }
    rownames(score_matrix) <- as.character(rs_vec)
    
    # 9) protected normalization
    score_matrix_normalized <- apply(score_matrix, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || (rng[2] - rng[1] == 0)) return(rep(0, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    })
    
    # 10) if large, do not show numbers
    show_numbers <- nrow(score_matrix) <= 60 && ncol(score_matrix) <= 20
    
    # ⬅️ store exactly what is plotted
    hd <- list(
      mat       = score_matrix,
      mat_norm  = score_matrix_normalized,
      rs        = rs_vec,
      title     = title_txt,
      show_numbers = show_numbers
    )
    last_heatmap(hd)
    
    heatmap_plot <- ggplotify::as.ggplot(function() {
      pheatmap::pheatmap(
        hd$mat_norm,
        cluster_rows = FALSE, cluster_cols = FALSE,
        display_numbers = if (hd$show_numbers) hd$mat else FALSE,
        labels_row = hd$rs,
        show_rownames = TRUE,
        show_colnames = TRUE,
        number_format = "%.2f",
        main = hd$title,
        fontsize_number = 10, fontsize_row = 8, fontsize_col = 10
      )
    })
    print(heatmap_plot)
  })
  
  
  # ===== Heatmap dynamic height wrapper =====
  output$dynamic_heatmap_plot <- renderUI({
    df <- score_input()
    n_rows <- if (!is.null(df)) nrow(df) else 0
    plot_height <- max(400, n_rows * 20)  # 20 px per row
    plotOutput("heatmap_plot", height = paste0(plot_height, "px"))
  })
  
  # 4.4) HEATMAP
  output$heatmap_scores <- renderPlot({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    # detected key columns
    snp_col <- session$userData$snp_col
    validate(need(!is.null(snp_col) && snp_col %in% names(df),
                  "Identifier column (SNP/rsID) not found."))
    
    chr_col <- session$userData$chr_col
    bp_col  <- session$userData$bp_col %||% NA
    
    # columns chosen by the user
    cols_sel <- if (identical(input$hm_mode, "Score")) input$hm_scores else input$hm_ranks
    
    # fallback: if <2, pick first useful numeric columns
    if (length(cols_sel) < 2) {
      num_ok <- names(df)[vapply(df, is.numeric, logical(1))]
      num_ok <- base::setdiff(num_ok, c("chr","CHR","BP","POS","Position","hg19_pos_1_based","hg38_pos_1_based"))
      cols_sel <- utils::head(num_ok, 4)
    }
    validate(need(length(cols_sel) >= 2, "Choose at least 2 columns for the heatmap"))
    
    # chromosome filter if applicable
    df2 <- df
    if (!is.null(chr_col) && chr_col %in% names(df2) && !is.null(input$hm_chr) && input$hm_chr != "all") {
      key <- suppressWarnings(as.numeric(input$hm_chr))
      key <- if (!is.na(key)) key else input$hm_chr
      df2 <- df2[df2[[chr_col]] == key, , drop = FALSE]
    }
    
    # order by chr + BP if present
    if (!is.na(bp_col) && bp_col %in% names(df2) && !is.null(chr_col) && chr_col %in% names(df2)) {
      df2 <- df2 %>% dplyr::arrange(.data[[chr_col]], .data[[bp_col]])
    }
    
    # limit top N
    n_show <- min(nrow(df2), input$hm_topn %||% 100)
    df2 <- utils::head(df2, n_show)
    
    # numeric matrix
    mat <- df2[, cols_sel, drop = FALSE]
    mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(as.character(x))))
    
    # scale 0–1 per column if requested
    if (isTRUE(input$hm_scale_cols)) {
      mat <- apply(as.matrix(mat), 2, function(x){
        rng <- range(x, na.rm = TRUE)
        if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
        (x - rng[1]) / (rng[2]-rng[1])
      })
    } else {
      mat <- as.matrix(mat)
    }
    
    # row labels
    rownames(mat) <- df2[[snp_col]]
    
    # plot
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE, cluster_cols = FALSE,
      fontsize_row = 8, fontsize_col = 9,
      main = if (identical(input$hm_mode,"Score")) "Heatmap (Scores)" else "Heatmap (Rankscores)",
      display_numbers = FALSE,
      silent = TRUE
    )
  }, height = function(){
    # 20 px per row, minimum 400
    df <- try(dbnsfp_norm_df(), silent = TRUE); if (inherits(df, "try-error")) return(400)
    chr_col <- session$userData$chr_col
    df2 <- df
    if (!is.null(chr_col) && chr_col %in% names(df2) && !is.null(input$hm_chr) && input$hm_chr != "all") {
      key <- suppressWarnings(as.numeric(input$hm_chr))
      key <- if (!is.na(key)) key else input$hm_chr
      df2 <- df2[df2[[chr_col]] == key, , drop = FALSE]
    }
    n <- min(nrow(df2), input$hm_topn %||% 100)
    max(400, n * 20)
  })
  
  ########################################################
  # Table (STEP 5.2)
  ########################################################
  
  # Heatmap Table
  # internal helper: try to use score_input(); otherwise fall back to dbnsfp_norm_df()
  .get_heatmap_df <- function() {
    df <- NULL
    if (exists("score_input", inherits = TRUE)) {
      try(df <- score_input(), silent = TRUE)
    }
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      df <- dbnsfp_norm_df()
    }
    df
  }
 
  output$heatmap_table <- DT::renderDT({
    
    df <- dbnsfp_norm_df()
    req(df)
    validate(need(nrow(df) > 0, "No data available."))
    
    # -----------------------------
    # 1) Filter by chromosome
    # -----------------------------
    chr_col <- session$userData$chr_col %||% "chr"
    
    if (!is.null(input$hm_chr) &&
        input$hm_chr != "all" &&
        chr_col %in% names(df)) {
      
      norm_chr <- function(x) {
        x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
        x[x=="23"] <- "X"
        x[x=="24"] <- "Y"
        x[x %in% c("M","MT","m","mT","mtdna","MTDNA")] <- "MT"
        toupper(x)
      }
      
      df[[chr_col]] <- norm_chr(df[[chr_col]])
      df <- df[df[[chr_col]] == input$hm_chr, , drop = FALSE]
    }
    
    validate(need(nrow(df) > 0, "No SNPs match chromosome filter."))
    
    # -----------------------------
    # 2) Base + Clinical columns
    # -----------------------------
    base_cols <- c("SNP", "chr", "BP", "ref", "alt", "genename")
    
    classes <- list(
      "Default" = c(
        "SIFT_score","SIFT4G_score",
        "Polyphen2_HDIV_score","Polyphen2_HVAR_score",
        "MutationTaster_score","MutationAssessor_score",
        "PROVEAN_score","REVEL_score",
        "GERP_91_mammals","phyloP17way_primate","phastCons17way_primate",
        "SIFT_converted_rankscore","SIFT4G_converted_rankscore",
        "Polyphen2_HDIV_rankscore","Polyphen2_HVAR_rankscore",
        "MutationAssessor_rankscore","PROVEAN_converted_rankscore",
        "REVEL_rankscore","GERP_91_mammals_rankscore",
        "phyloP17way_primate_rankscore","phastCons17way_primate_rankscore"
      ),
      
      "ClinVar" = c(
        "clinvar_id","clinvar_clnsig","clinvar_trait",
        "clinvar_review","clinvar_hgvs","clinvar_var_source",
        "clinvar_MedGen_id","clinvar_OMIM_id","clinvar_Orphanet_id"
      ),
      
      "OMIM / Orphanet" = c(
        "MIM_id","OMIM_id","MIM_phenotype_id","MIM_disease",
        "Orphanet_disorder_id","Orphanet_disorder",
        "Orphanet_association_type","Disease_description",
        "Trait_association_GWAS"
      ),
      
      "Pathogenicity" = c(
        "SIFT_score","SIFT_pred","Polyphen2_HDIV_score","Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score","Polyphen2_HVAR_pred","MutationTaster_score",
        "MutationTaster_pred","MutationAssessor_score","MutationAssessor_pred",
        "PROVEAN_score","PROVEAN_pred","VEST4_score","MetaSVM_score",
        "MetaSVM_pred","MetaLR_score","MetaLR_pred","MetaRNN_score",
        "MetaRNN_pred","M_CAP_score","M_CAP_pred","REVEL_score",
        "MutPred_score","MVP_score","gMVP_score","MPC_score",
        "PrimateAI_score","PrimateAI_pred","DEOGEN2_score","DEOGEN2_pred",
        "BayesDel_addAF_score","BayesDel_addAF_pred",
        "BayesDel_noAF_score","BayesDel_noAF_pred",
        "ClinPred_score","ClinPred_pred",
        "LIST_S2_score","LIST_S2_pred",
        "VARITY_R_score","VARITY_ER_score",
        "ESM1b_score","ESM1b_pred","AlphaMissense_score",
        "AlphaMissense_pred","PHACTboost_score","MutFormer_score",
        "MutScore_score"
      ),
      
      "Conservation" = c(
        "GERP++_RS","phyloP100way_vertebrate",
        "phyloP470way_mammalian","phyloP17way_primate",
        "phastCons100way_vertebrate","phastCons470way_mammalian",
        "phastCons17way_primate"
      ),
      
      "Functional impact" = c(
        "CADD_raw","CADD_phred",
        "DANN_score","Eigen_raw_coding","Eigen_phred_coding"
      ),
      
      "LoF / Haploinsufficiency" = c(
        "P_HI","HIPred","GHIS","ClinGen_Haploinsufficiency_Score",
        "ClinGen_Haploinsufficiency_Description",
        "ClinGen_Haploinsufficiency_PMID",
        "ClinGen_Haploinsufficiency_Disease","P_rec",
        "Known_rec_info","LoF_FDR_ExAC","ExAC_pLI","ExAC_pRec",
        "ExAC_pNull","gnomAD_pLI","gnomAD_pRec","gnomAD_pNull",
        "ExAC_del_score","ExAC_dup_score","ExAC_cnv_score","ExAC_cnv_flag"
      ),
      
      "Gene damage" = c(
        "GDI","GDI_Phred",
        "Gene_damage_prediction_all_disease_causing_genes",
        "Gene_damage_prediction_all_Mendelian_disease_causing_genes",
        "Gene_damage_prediction_Mendelian_AD_disease_causing_genes",
        "Gene_damage_prediction_Mendelian_AR_disease_causing_genes",
        "Gene_damage_prediction_all_PID_disease_causing_genes",
        "Gene_damage_prediction_PID_AD_disease_causing_genes",
        "Gene_damage_prediction_PID_AR_disease_causing_genes",
        "Gene_damage_prediction_all_cancer_disease_causing_genes",
        "Gene_damage_prediction_cancer_recessive_disease_causing_genes",
        "Gene_damage_prediction_cancer_dominant_disease_causing_genes",
        "LoFtool_score","Essential_gene",
        "Essential_gene_CRISPR","Essential_gene_CRISPR2",
        "Essential_gene_gene_trap","Gene_indispensability_score",
        "Gene_indispensability_pred"
      )
    )
    
    choice <- input$clinical_class
    cols <- c(base_cols, intersect(classes[[choice]], names(df)))
    df <- df[, unique(cols), drop = FALSE]
    
    # =====================================================
    # 3) MAKE ALL INTERACTIVE LINKS
    # =====================================================
    
    ## --- GeneCards ---
    if ("genename" %in% names(df)) {
      df$genename <- ifelse(
        is.na(df$genename) | df$genename == "",
        "",
        sprintf(
          "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s'
           target='_blank' style='color:#2a4bd7;'>%s</a>",
          df$genename, df$genename
        )
      )
    }
    
    ## --- dbSNP ---
    if ("SNP" %in% names(df)) {
      df$SNP <- ifelse(
        grepl("^rs", df$SNP, ignore.case = TRUE),
        sprintf(
          "<a href='https://www.ncbi.nlm.nih.gov/snp/%s'
         target='_blank' style='color:#0066cc;'>%s</a>",
          df$SNP, df$SNP
        ),
        df$SNP
      )
    }
    
    ## --- ClinVar ---
    if ("clinvar_id" %in% names(df)) {
      df$clinvar_id <- ifelse(
        !is.na(df$clinvar_id) & df$clinvar_id != "",
        sprintf(
          "<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/%s'
           target='_blank' style='color:#d9534f;'>%s</a>",
          df$clinvar_id, df$clinvar_id
        ),
        ""
      )
    }
    
    ## --- OMIM ---
    if ("MIM_id" %in% names(df)) {
      df$MIM_id <- ifelse(
        !is.na(df$MIM_id) & df$MIM_id != "",
        sprintf(
          "<a href='https://omim.org/entry/%s'
         target='_blank' style='color:#0066cc;'>%s</a>",
          df$MIM_id, df$MIM_id
        ),
        ""
      )
    }
    
    if ("OMIM_id" %in% names(df)) {
      df$OMIM_id <- ifelse(
        !is.na(df$OMIM_id) & df$OMIM_id != "",
        sprintf(
          "<a href='https://omim.org/entry/%s'
         target='_blank' style='color:#0066cc;'>%s</a>",
          df$OMIM_id, df$OMIM_id
        ),
        ""
      )
    }
    
    ## --- UCSC: chr:BP ---
    if (all(c("chr","BP") %in% names(df))) {
      norm_chr_ucsc <- function(x) {
        gsub("^chr","",x,ignore.case=TRUE)
      }
      
      df$UCSC <- sprintf(
        "<a href='https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr%s:%s-%s'
          target='_blank' style='color:#0066cc;'>view</a>",
        norm_chr_ucsc(df$chr),
        df$BP,
        df$BP
      )
      
      df <- df[, c(colnames(df), "UCSC")]  # add as last column
    }
    
    # -----------------------------------------------------
    # 4) Render final interactive DT
    # -----------------------------------------------------
    DT::datatable(
      df,
      escape = FALSE,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })
  
  # ================================================================
  
  output$dl_heatmap_png <- downloadHandler(
    filename = function() sprintf("heatmap_%s.png", tolower(input$scoreorrank %||% "score")),
    contentType = "image/png",
    content = function(file) {
      # 1) use the last plotted heatmap
      hd <- last_heatmap()
      
      # 2) if it does not exist (e.g. tab not opened), rebuild with same logic
      if (is.null(hd) || is.null(hd$mat) || is.null(hd$mat_norm)) {
        df <- tryCatch(score_input(), error = function(e) NULL)
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
          # reconstruction with the same rules as renderPlot
          snp_col <- if ("rsid_final" %in% names(df)) "rsid_final" else (attr(df, "snp_col") %||% "SNP")
          chr_col <- attr(df, "chr_col") %||% "chr"
          bp_col  <- attr(df, "bp_col")  %||% "BP"
          nms <- names(df)
          if (identical(input$scoreorrank, "Score")) {
            metric_cols <- c(
              nms[grepl("_score$", nms, perl = TRUE)],
              base::intersect(c("GERP_91_mammals","phyloP17way_primate","phastCons17way_primate"), nms)
            )
            title_txt <- "Heatmap Scores"
          } else {
            metric_cols <- c(
              nms[grepl("_rankscore$", nms, perl = TRUE)],
              base::intersect("GERP++_RS_rankscore", nms)
            )
            title_txt <- "Heatmap Rank-scores"
          }
          id_cols <- base::intersect(unique(c(snp_col, "rsid_final", chr_col, bp_col, "ref","alt","genename")), nms)
          metric_cols <- setdiff(unique(metric_cols), id_cols)
          if (length(metric_cols) >= 2) {
            df_hm <- df[, c(id_cols, metric_cols), drop = FALSE]
            for (cc in metric_cols) df_hm[[cc]] <- suppressWarnings(as.numeric(as.character(df_hm[[cc]])))
            mat <- as.matrix(df_hm[, metric_cols, drop = FALSE])
            rs_vec <- if ("rsid_final" %in% names(df_hm)) df_hm[["rsid_final"]] else {
              if (snp_col %in% names(df_hm)) df_hm[[snp_col]] else NA
            }
            if (is.null(rs_vec) || all(is.na(rs_vec) | rs_vec == "")) {
              if (all(c(chr_col, bp_col) %in% names(df_hm))) rs_vec <- paste0(df_hm[[chr_col]], ":", df_hm[[bp_col]])
              else rs_vec <- paste0("var_", seq_len(nrow(df_hm)))
            }
            rownames(mat) <- as.character(rs_vec)
            mat_norm <- apply(mat, 2, function(x) {
              rng <- range(x, na.rm = TRUE)
              if (!is.finite(rng[1]) || !is.finite(rng[2]) || (rng[2] - rng[1] == 0)) return(rep(0, length(x)))
              (x - rng[1]) / (rng[2] - rng[1])
            })
            hd <- list(mat = mat, mat_norm = mat_norm, rs = rs_vec, title = title_txt,
                       show_numbers = (nrow(mat) <= 60 && ncol(mat) <= 20))
          }
        }
      }
      
      # 3) if still no data, emit informative PNG
      if (is.null(hd) || is.null(hd$mat) || is.null(hd$mat_norm) ||
          nrow(hd$mat) == 0 || ncol(hd$mat) == 0) {
        png(file, width = 1400L, height = 420L, res = 120)
        par(mar = c(0,0,0,0)); plot.new()
        text(0.5, 0.55, "No data available to download the heatmap.", cex = 1.5)
        text(0.5, 0.35, "Open the heatmap tab to generate it and then download again.", cex = 1.1)
        dev.off()
        return(invisible())
      }
      
      # 4) dynamic size like the plot
      height_in <- max(6, nrow(hd$mat) * 0.20)
      width_in  <- max(9, ncol(hd$mat) * 0.50)
      
      # 5) direct saving with pheatmap
      pheatmap::pheatmap(
        hd$mat_norm,
        cluster_rows = FALSE, cluster_cols = FALSE,
        display_numbers = if (isTRUE(hd$show_numbers)) hd$mat else FALSE,
        labels_row = hd$rs,
        show_rownames = TRUE,
        main = hd$title,
        fontsize_number = 10, fontsize_row = 8, fontsize_col = 10,
        filename = file,
        width = width_in,
        height = height_in
      )
    }
  )
  
  ########################################################
  # Manhattan and genes (STEP 5.3)
  ########################################################
  
  # --- Helper: chromosome lengths hg38 (GRCh38, primary) ---
  .chr_lengths_hg38 <- function() {
    lens <- c(
      `1`=248956422, `2`=242193529, `3`=198295559, `4`=190214555, `5`=181538259,
      `6`=170805979, `7`=159345973, `8`=145138636, `9`=138394717, `10`=133797422,
      `11`=135086622, `12`=133275309, `13`=114364328, `14`=107043718, `15`=101991189,
      `16`=90338345,  `17`=83257441,  `18`=80373285,  `19`=58617616,  `20`=64444167,
      `21`=46709983,  `22`=50818468,  `X`=156040895,  `Y`=57227415,  `MT`=16569
    )
    ord <- c(as.character(1:22), "X", "Y", "MT")
    df <- data.frame(
      chr = factor(names(lens), levels = ord),
      len = as.numeric(lens),
      stringsAsFactors = FALSE
    )
    df$chr_cum <- cumsum(df$len) - df$len
    df$center  <- df$chr_cum + df$len/2
    df
  }
  
  .ref_hg38 <- .chr_lengths_hg38()
  
  # --- Manhattan of multiple metrics; X axis ALWAYS spans 1–22, X, Y, MT ---
  output$manhattan_scores <- renderPlotly({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    chr_col <- session$userData$chr_col %||% "chr"
    bp_col  <- session$userData$bp_col  %||% "BP"
    snp_col <- session$userData$snp_col %||% "SNP"
    
    validate(need(all(c(chr_col, bp_col) %in% names(df)),
                  "Missing chr/BP columns."))
    
    # ===================================================
    # 1) Select: Score / Rankscore 
    # ===================================================
    # 1) Score vs Rankscore
    user_type <- input$scoreorrank %||% "Score"
    
    if (user_type == "Score") {
      metric_cols <- grep("_score$", names(df), value = TRUE)
    } else {
      metric_cols <- grep("_rankscore$", names(df), value = TRUE)
    }
    
    metric_cols <- metric_cols[colSums(!is.na(df[metric_cols])) > 0]
    
    validate(need(length(metric_cols) > 0,
                  paste("No", user_type, "columns found with valid values.")))
    
    # 2) Aplicar selecció del picker
    selected <- input$ms_metrics
    
    # 👉 Si l’usuari no selecciona res, usar només la primera
    if (is.null(selected) || !length(selected)) {
      selected <- metric_cols[1]
    }
    
    metric_cols <- intersect(metric_cols, selected)
    
    validate(need(length(metric_cols) > 0,
                  paste("Selected", user_type, "metrics contain no valid values.")))
    
    # ===================================================
    # 2) Applay selection at checkboxGroupInput
    # ===================================================
    
    selected <- input$ms_metrics
    
    # Si no hi ha cap selecció → agafar la PRIMERA mètrica disponible
    if (is.null(selected) || length(selected) == 0) {
      selected <- metric_cols[1]     # <<--- CLAR, SIMPLE i FIABLE
    }
    
    metric_cols <- intersect(metric_cols, selected)
    validate(need(length(metric_cols) > 0,
                  "No valid metrics selected for plotting."))
    # ===================================================
    # 3) Chromosome normalization
    # ===================================================
    norm_chr <- function(x) {
      x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
      x[x=="23"] <- "X"
      x[x=="24"] <- "Y"
      x[x %in% c("M","MT","m","mT","mtdna","MTDNA")] <- "MT"
      toupper(x)
    }
    
    ref <- .ref_hg38
    
    # ===================================================
    # 4) Long format
    # ===================================================
    dfl <- df |>
      dplyr::select(dplyr::all_of(c(chr_col, bp_col, snp_col, metric_cols))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(metric_cols),
        names_to = "metric", values_to = "value"
      ) |>
      dplyr::mutate(
        chrN = norm_chr(.data[[chr_col]]),
        BP   = suppressWarnings(as.numeric(.data[[bp_col]]))
      ) |>
      dplyr::filter(!is.na(BP), !is.na(value))
    
    validate(need(nrow(dfl) > 0,
                  "No valid data for Manhattan plot."))
    
    # ===================================================
    # 5) Cummulative position
    # ===================================================
    dfl <- dfl |>
      dplyr::mutate(chrN = factor(chrN, levels = levels(ref$chr))) |>
      dplyr::inner_join(ref |> dplyr::select(chr, chr_cum),
                        by = c("chrN" = "chr")) |>
      dplyr::mutate(BPcum = BP + chr_cum)
    
    validate(need(!any(is.na(dfl$BPcum)),
                  "Invalid cumulative BP positions."))
    
    # ===================================================
    # ===================================================
    axis_df <- dfl |>
      dplyr::group_by(chrN) |>
      dplyr::summarise(
        center = mean(BPcum, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(center)) |>
      dplyr::arrange(center)
    
    axis_breaks <- axis_df$center
    axis_labels <- paste0("chr", axis_df$chrN)
    
    validate(need(length(axis_breaks) == length(axis_labels),
                  "Internal X-axis mismatch."))
    
    # ===================================================
    # 7) Plot Manhattan
    # ===================================================
    levs <- sort(unique(dfl$metric))
    pal  <- viridisLite::viridis(length(levs))
    
    p <- ggplot(dfl, aes(
      x = BPcum, y = value, color = metric,
      text = paste0(
        "SNP: ", .data[[snp_col]],
        "<br>Chr: ", chrN,
        "<br>Pos: ", BP,
        "<br>", metric, ": ", signif(value, 3)
      )
    )) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = stats::setNames(pal, levs)) +
      scale_x_continuous(
        breaks = axis_breaks,
        labels = axis_labels,
        expand = c(0.01, 0.01)
      ) +
      labs(
        x = "Genome (1–22, X, Y, MT)",
        y = user_type,
        title = paste("dbNSFP — Manhattan (", user_type, ")", sep = "")
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "right") +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    p_pl <- plotly::ggplotly(p, tooltip = "text", source = "ms")
    plotly::event_register(p_pl, "plotly_click")
    p_pl
  })
  
  output$ms_metric_picker <- renderUI({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    type <- input$scoreorrank %||% "Score"
    
    # 1) Llista de mètriques segons Tipus
    all_metrics <- if (type == "Score") {
      grep("_score$", names(df), value = TRUE)
    } else {
      grep("_rankscore$", names(df), value = TRUE)
    }
    
    # Només columnes numèriques amb algun valor
    all_metrics <- all_metrics[vapply(df[all_metrics], is.numeric, logical(1))]
    all_metrics <- all_metrics[colSums(!is.na(df[all_metrics])) > 0]
    
    if (!length(all_metrics)) return(NULL)
    
    # 2) Mapes de classes per SCORE i per RANKSCORE
    class_map_score <- list(
      "Default" = c(
        "SIFT_score","SIFT4G_score",
        "Polyphen2_HDIV_score","Polyphen2_HVAR_score",
        "MutationTaster_score","MutationAssessor_score",
        "PROVEAN_score","REVEL_score"
      ),
      "Pathogenicity" = c(
        "VEST4_score","MetaSVM_score","MetaLR_score","MetaRNN_score",
        "M_CAP_score","MutPred_score","MVP_score","gMVP_score","MPC_score",
        "PrimateAI_score","DEOGEN2_score","BayesDel_addAF_score",
        "BayesDel_noAF_score","ClinPred_score","LIST_S2_score",
        "VARITY_R_score","VARITY_ER_score","ESM1b_score",
        "AlphaMissense_score","PHACTboost_score","MutFormer_score",
        "MutScore_score"
      ),
      "Functional impact" = c("CADD_raw","CADD_phred","DANN_score"),
      "LoF / Haploinsufficiency" = c(
        "HIPred_score","ExAC_del_score","ExAC_dup_score",
        "ExAC_cnv_score","LoFtool_score"
      ),
      "Gene damage" = c(
        "GDI","GDI_Phred","Gene_indispensability_score"
      )
    )
    
    class_map_rank <- list(
      "Default" = c(
        "SIFT_converted_rankscore","SIFT4G_converted_rankscore",
        "Polyphen2_HDIV_rankscore","Polyphen2_HVAR_rankscore",
        "MutationAssessor_rankscore","PROVEAN_converted_rankscore",
        "REVEL_rankscore","GERP_91_mammals_rankscore",
        "phyloP17way_primate_rankscore","phastCons17way_primate_rankscore"
      ),
      "Pathogenicity" = c(
        "VEST4_rankscore","MetaSVM_rankscore","MetaLR_rankscore",
        "MetaRNN_rankscore","M_CAP_rankscore","MutPred_rankscore",
        "MVP_rankscore","gMVP_rankscore","MPC_rankscore",
        "PrimateAI_rankscore","DEOGEN2_rankscore",
        "BayesDel_addAF_rankscore","BayesDel_noAF_rankscore",
        "ClinPred_rankscore","LIST_S2_rankscore",
        "VARITY_R_rankscore","VARITY_ER_rankscore",
        "ESM1b_rankscore","AlphaMissense_rankscore",
        "PHACTboost_rankscore","MutFormer_rankscore",
        "MutScore_rankscore"
      ),
      "Functional impact" = c("CADD_raw_rankscore","DANN_rankscore"),
      "LoF / Haploinsufficiency" = c(
        "HIPred_score",           # HIPred no té "_rankscore" en dbNSFP
        "ExAC_del_score","ExAC_dup_score","ExAC_cnv_score"
      ),
      "Gene damage" = c(
        "GDI_Phred"               # l’índex “rank” ja està implícit
      )
    )
    
    class_map <- if (type == "Score") class_map_score else class_map_rank
    
    # 3) Construir choices agrupats
    choices_list <- list()
    used <- character(0)
    
    for (grp in names(class_map)) {
      present <- intersect(class_map[[grp]], all_metrics)
      if (length(present)) {
        choices_list[[grp]] <- stats::setNames(present, present)
        used <- c(used, present)
      }
    }
    
    remaining <- setdiff(all_metrics, used)
    if (length(remaining)) {
      choices_list[["Other"]] <- stats::setNames(remaining, remaining)
    }
    
    # 4) Per defecte: primera mètrica de la llista
    default_metric <- if (length(all_metrics)) all_metrics[1] else NULL
    
    selectizeInput(
      "ms_metrics",
      "Select metrics to plot",
      choices  = choices_list,
      selected = default_metric,
      multiple = TRUE,
      options  = list(plugins = list("remove_button"))
    )
  })
  
  
  # --- click on the Manhattan ---
  ms_click <- reactiveVal(NULL)
  
  observeEvent(plotly::event_data("plotly_click", source = "ms"), {
    ed <- plotly::event_data("plotly_click", source = "ms")
    if (is.null(ed) || is.null(ed$x)) return()
    bp_cum <- as.numeric(ed$x)
    
    ref <- .ref_hg38
    idx <- max(which(ref$chr_cum <= bp_cum))
    if (length(idx) == 0 || !is.finite(idx)) return()
    
    chr <- as.character(ref$chr[idx])
    pos <- round(bp_cum - ref$chr_cum[idx])
    
    # store the “raw” region; refine it later with the window
    ms_click(list(chr = chr, pos = pos))
  }, ignoreInit = TRUE)
  
  # if you change the window and there is no click yet, set a default one
  observeEvent(input$rg_window, ignoreInit = FALSE, {
    if (is.null(ms_click())) {
      ms_click(list(chr = "1", pos = floor(.ref_hg38$len[.ref_hg38$chr == "1"]/2)))
    }
  })
  
#  output$rg_pos_txt <- renderText({
#    p <- ms_click(); wkb <- as.integer(input$rg_window %||% 250)
#    if (is.null(p)) "Click on the Manhattan plot to center the region."
#    else sprintf("Center: chr%s:%s  |  Window: ±%d kb",
#                 p$chr, format(p$pos, big.mark = ","), wkb)
#  })
  
  
  # UCSC viewer as clickable link
  output$ucsc_viewer <- renderUI({
    p <- ms_click()
    wkb <- as.integer(input$rg_window %||% 250)
    
    if (is.null(p)) {
      return(div(class = "text-muted",
                 "🔎 Click on the Manhattan plot to generate the UCSC (hg38) link."))
    }
    
    chr  <- p$chr
    pos  <- p$pos
    start <- max(1L, pos - wkb * 1000L)
    end   <- pos + wkb * 1000L
    
    # European mirror (usually most stable)
    base_url <- "https://genome-euro.ucsc.edu/cgi-bin/hgTracks"
    
    # 🔥  highlight = chr:start-end (vertical line)
    ucsc_url <- sprintf(
      "%s?db=hg38&position=chr%s:%d-%d&highlight=hg38.chr%s:%d-%d",
      base_url, chr, start, end, chr, pos, pos
    )
    
    tags$div(
      tags$p("📍 Selected region: ",
             tags$b(sprintf("chr%s:%d-%d", chr, start, end))),
      
      tags$small("Open this region in UCSC: "),
      tags$a(
        href   = ucsc_url,
        target = "_blank",
        class  = "btn btn-primary btn-sm",
        " UCSC"
      ),
      tags$br(), tags$br(),
      tags$code(ucsc_url)
    )
  })
  
  
  # ==========================
  #  GO Enrichment (Step 5.4)
  # ==========================
  
  # --- helper to shorten terms (same as before) ---
  shorten_go <- function(x, max = 45) {
    if (is.factor(x)) x <- as.character(x)
    x <- gsub("^GO:[0-9]+\\s*;?\\s*", "", x)
    x <- gsub("\\bregulation of\\b", "reg. of", x)
    x <- gsub("\\bpositive\\b", "pos.", x)
    x <- gsub("\\bnegative\\b", "neg.", x)
    x <- gsub("\\bcellular\\b",  "cell.", x)
    x <- gsub("\\bprocess\\b",   "proc.", x)
    x <- gsub("\\bcomponent\\b", "comp.", x)
    x <- gsub("\\bactivity\\b",  "act.", x)
    x <- stringr::str_squish(x)
    x <- stringr::str_to_sentence(x)
    x <- ifelse(nchar(x) > max, paste0(substr(x, 1, max - 1), "…"), x)
    x
  }
  
  # --- data ready for plotting (reactive) ---
  go_top_df <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    go_cols <- session$userData$go_cols %||% character(0)
    validate(need(length(go_cols) >= 1, "No GO columns detected."))
    
    pcut <- if (!is.null(input$go_kegg_pcut) && length(input$go_kegg_pcut) == 1)
      input$go_kegg_pcut else 0.05
    
    df_go <- df[, go_cols, drop = FALSE] |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "ontology", values_to = "term_raw") |>
      dplyr::filter(!is.na(term_raw) & term_raw != "") |>
      tidyr::separate_rows(term_raw, sep = ";", convert = FALSE) |>
      dplyr::mutate(term_full = stringr::str_squish(term_raw)) |>
      dplyr::filter(nzchar(term_full))
    
    validate(need(nrow(df_go) > 0, "No GO terms available."))
    
    go_counts <- df_go |>
      dplyr::group_by(ontology, term_full) |>
      dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(ontology) |>
      dplyr::mutate(
        pval = 1 / rank(-count, ties.method = "first"),
        padj = p.adjust(pval, method = "BH")
      ) |>
      dplyr::ungroup()
    
    go_top <- go_counts |>
      dplyr::filter(padj <= pcut) |>
      dplyr::group_by(ontology) |>
      dplyr::slice_max(order_by = count, n = 10) |>
      dplyr::ungroup()
    
    validate(need(nrow(go_top) > 0, "No significant GO terms for the chosen FDR."))
    
    go_top |>
      dplyr::mutate(
        ontology = dplyr::recode(ontology,
                                 "GO_biological_process" = "BP",
                                 "GO_molecular_function" = "MF",
                                 "GO_cellular_component" = "CC"),
        term_short = shorten_go(term_full, max = 45),
        term_short = stringr::str_wrap(term_short, width = 28)
      ) |>
      dplyr::group_by(ontology) |>
      dplyr::arrange(count, .by_group = TRUE) |>
      dplyr::mutate(term_f = forcats::fct_inorder(term_short)) |>
      dplyr::ungroup()
  })
  
  
  output$go_bar <- renderPlot({
    dat <- go_top_df()  # <- your reactive
    validate(need(!is.null(dat) && nrow(dat) > 0, "There are no GO terms to plot."))
    
    # cache for download
    last_go(dat)
    
    ggplot(dat, aes(x = term_f, y = count, fill = ontology)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~ontology, scales = "free_y") +
      labs(x = NULL, y = "Gene count",
           title = "GO terms (abbreviated) — simulated FDR ≤ threshold") +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(lineheight = 0.9)
      )
  },
  height = function() {
    n <- tryCatch(nrow(go_top_df()), error = function(e) 0)
    max(450, 10 * n)
  })
  
  output$dl_go_png <- downloadHandler(
    filename = function() "GO_bar.png",
    contentType = "image/png",
    content = function(file) {
      dat <- isolate(last_go())
      if (is.null(dat) || !nrow(dat)) {
        # rebuild in case the tab was never opened
        dat <- tryCatch(isolate(go_top_df()), error = function(e) NULL)
      }
      
      if (is.null(dat) || !nrow(dat)) {
        png(file, width = 1200, height = 400, res = 120)
        par(mar = c(0,0,0,0)); plot.new()
        text(0.5, 0.5, "There are no GO data to download.", cex = 1.4)
        dev.off()
        return(invisible())
      }
      
      # dynamic height so labels fit
      h <- max(600, 40 * nrow(dat))
      png(file, width = 1600, height = h, res = 150)
      print(
        ggplot(dat, aes(x = term_f, y = count, fill = ontology)) +
          geom_col() +
          coord_flip() +
          facet_wrap(~ontology, scales = "free_y") +
          labs(x = NULL, y = "Gene count",
               title = "GO terms (abbreviated) — simulated FDR ≤ threshold") +
          theme_minimal(base_size = 16) +
          theme(legend.position = "none",
                strip.text = element_text(face = "bold"),
                axis.text.y = element_text(lineheight = 0.9))
      )
      dev.off()
    }
  )
  
  
  output$go_table <- DT::renderDT({
    dat <- go_top_df()
    # Selection of “nice” columns
    out <- dat |>
      dplyr::select(
        Ontology = ontology,
        `Term (abbreviated)` = term_short,
        `Full term` = term_full,
        Count = count,
        FDR = padj
      )
    
    DT::datatable(
      out,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 15,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        order = list(list(3, "desc")) # order by Count
      ),
      extensions = c("Buttons")
    ) |>
      DT::formatSignif("FDR", digits = 3)
  })

  # ==========================
  #  KEGG Enrichment (Step 5.5)
  # ==========================
  
  # 5.5) KEGG enrichment with clusterProfiler
  # Requires: org.Hs.eg.db, clusterProfiler
  kegg_result <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    gene_col <- session$userData$gene_col
    validate(need(!is.null(gene_col) && gene_col %in% names(df),
                  "There is no gene symbol column (genename)."))
    
    genes <- df[[gene_col]]
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(NULL)
    
    genes <- unique(trimws(unlist(strsplit(genes, ";", fixed = TRUE))))
    genes <- genes[nzchar(genes)]
    validate(need(length(genes) > 0, "There are no valid gene symbols."))
    
    suppressWarnings({
      gmap <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db)
    })
    validate(need(!is.null(gmap) && nrow(gmap) > 0, "Could not map genes to ENTREZ."))
    
    ek <- clusterProfiler::enrichKEGG(
      gene = unique(gmap$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 1
    )
    # normalize to a data.frame-like result in case the method changes
    if (is.null(ek) || is.null(ek@result) || !nrow(ek@result)) return(NULL)
    ek
  })
  
  output$kegg_bar <- renderPlot({
    ek <- kegg_result()
    validate(need(!is.null(ek) && nrow(ek@result) > 0, "No KEGG results."))
    
    df <- ek@result
    # FDR + top-N
    df <- df %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      dplyr::filter(p.adjust <= (input$go_kegg_pcut %||% 0.05)) %>%
      dplyr::slice_head(n = input$kegg_top %||% 15)
    
    validate(need(nrow(df) > 0, "No KEGG terms below the selected FDR."))
    
    pcut <- if (!is.null(input$go_kegg_pcut) && length(input$go_kegg_pcut) == 1) input$go_kegg_pcut else 0.05
    
    # Cache for download
    last_kegg(list(data = df, pcut = pcut))
    
    ggplot2::ggplot(df, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust))) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "KEGG enrichment (top)", x = NULL, y = "-log10(FDR)") +
      ggplot2::theme_minimal(base_size = 12)
  })
  
  output$dl_kegg_png <- downloadHandler(
    filename = function() {
      sprintf("KEGG_bar_FDR_%.3f.png", input$go_kegg_pcut %||% 0.05)
    },
    contentType = "image/png",
    content = function(file) {
      dat <- isolate(last_kegg())
      # rebuild if the plot has not been drawn yet
      if (is.null(dat) || is.null(dat$data) || nrow(dat$data) == 0) {
        ek <- tryCatch(isolate(kegg_result()), error = function(e) NULL)
        if (!is.null(ek) && nrow(ek@result) > 0) {
          tmp <- ek@result %>%
            dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
            dplyr::filter(p.adjust <= (input$go_kegg_pcut %||% 0.05)) %>%
            dplyr::slice_head(n = input$kegg_top %||% 15)
          if (nrow(tmp) > 0) dat <- list(data = tmp, pcut = input$go_kegg_pcut %||% 0.05)
        }
      }
      
      if (is.null(dat) || is.null(dat$data) || nrow(dat$data) == 0) {
        png(file, width = 1200, height = 420, res = 120)
        par(mar = c(0,0,0,0)); plot.new()
        text(0.5, 0.55, "There are no KEGG data to download.", cex = 1.5)
        text(0.5, 0.35, "Adjust the FDR threshold or the top-N.", cex = 1.1)
        dev.off()
        return(invisible())
      }
      
      df   <- dat$data
      pcut <- dat$pcut
      # Dynamic height (more bars = taller)
      h <- max(600, 40 * nrow(df))
      
      png(file, width = 1600, height = h, res = 150)
      print(
        ggplot2::ggplot(df, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust))) +
          ggplot2::geom_col(fill = "steelblue") +
          ggplot2::coord_flip() +
          ggplot2::labs(
            title = sprintf("KEGG enrichment (top) — FDR ≤ %.3f", pcut),
            x = NULL, y = "-log10(FDR)"
          ) +
          ggplot2::theme_minimal(base_size = 16)
      )
      dev.off()
    }
  )
  
  # ================================================================
  
  output$kegg_table <- DT::renderDT({
    ek <- kegg_result()
    validate(need(!is.null(ek) && !is.null(ek@result) && nrow(ek@result) > 0,
                  "No KEGG results."))
    
    res <- as.data.frame(ek@result, stringsAsFactors = FALSE)
    
    # Ensure FDR column:
    if (!"p.adjust" %in% names(res) && "pvalue" %in% names(res)) {
      res$p.adjust <- p.adjust(res$pvalue, method = "BH")
    }
    # Remove rows without FDR
    res <- res[!is.na(res$p.adjust), , drop = FALSE]
    
    # FDR threshold from the slider
    pcut <- if (!is.null(input$go_kegg_pcut) && length(input$go_kegg_pcut) == 1)
      input$go_kegg_pcut else 0.05
    
    res <- res[res$p.adjust <= pcut, , drop = FALSE]
    validate(need(nrow(res) > 0, "No KEGG pathways under the selected FDR."))
    
    # Clickable KEGG link + “nice” columns
    df <- res |>
      dplyr::mutate(
        KEGG = sprintf("<a href='https://www.kegg.jp/pathway/%s' target='_blank'>%s</a>", ID, ID)
      ) |>
      dplyr::arrange(p.adjust) |>
      dplyr::transmute(
        KEGG,                                      # link
        Description   = Description,
        `Gene ratio`  = GeneRatio,
        `Background (BgRatio)` = BgRatio,
        Count         = Count,
        `p-value`     = pvalue,
        FDR           = p.adjust,
        qvalue        = qvalue,
        Genes         = geneID
      )
    
    DT::datatable(
      df,
      escape   = FALSE,   # keep links
      rownames = FALSE,
      filter   = "top",
      extensions = c("Buttons"),
      options = list(
        pageLength = 15,
        scrollX    = TRUE,
        dom        = "Bfrtip",
        buttons    = c("copy", "csv", "excel"),
        order      = list(list(6, "asc")) # order by FDR (col 7 in 1-index → 6 in 0-index)
      )
    ) |>
      DT::formatSignif(c("p-value", "FDR", "qvalue"), digits = 3)
  })

}

shinyApp(ui, server)

