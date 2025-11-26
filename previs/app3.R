# app.R
options(shiny.maxRequestSize = 1024*1024^2)  # ajústalo si lo necesitas

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
  # si no hay punto y sí hay coma -> cambia coma por punto (1,30E-06 -> 1.30E-06)
  needs_swap <- !grepl("\\.", p_chr) & grepl(",", p_chr)
  p_chr[needs_swap] <- gsub(",", ".", p_chr[needs_swap], fixed = TRUE)
  suppressWarnings(as.numeric(p_chr))
}

# operador OR seguro (no vectoriza)
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) b else a
}

# Siempre usa intersect de base para evitar conflictos con IRanges/dplyr
# Helpers (si no los tienes ya)
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
  title = "GWAS Inspector",
  id = "topnav",
  
  # ================= TAB 1: Análisis (Etapas 1–3) =================
  tabPanel(
    "Análisis (Etapas 1–3)",
    tags$style(HTML("
      .tab-content { padding-top: 0 !important; }
      .main-panel h3:first-child { margin-top: 0 !important; }
    ")),
    sidebarLayout(
      sidebarPanel(width = 4,
                   h3("Etapa 1 · Lectura & selección"),
                   fileInput("file", "Carga GWAS (TSV/CSV)", accept = c(".tsv",".txt",".csv")),
                   checkboxInput("header", "Primera fila con nombres", TRUE),
                   radioButtons("sep", "Separador",
                                c("Tab \\t" = "\t", "Coma ," = ",", "Punto y coma ;" = ";"),
                                selected = "\t"),
                   sliderInput("pthr", "Umbral -log10(P)", min = 2, max = 20, value = 7.3, step = 0.1),
                   numericInput("flank", "Flanco (+/- bp)", value = 10000, min = 0, step = 1000),
                   actionButton("build_ranges", "➊ Generar selected_intervals"),
                   h4("Preview selected_intervals.range"),
                   verbatimTextOutput("ranges_preview"),
                   tags$hr(),
                   
                   h3("Etapa 2 · VCF con PLINK 1.9"),
                   helpText("Coloca Plink y los ref.* (bed/bim/fam) en la carpeta www/\n(Opcional incluye la ruta a Plink y a la carpeta con el binario de referencia)"),
                   textInput("plink_path", "Binario PLINK \n(Opcional: Añadir ruta a la carpeta con el binario de referencia)", value = "www/plink19"),
                   textInput("ref_bfile", "Prefijo --bfile (auto si vacío)", value = ""),
                   numericInput("threads19", "Hilos", 4, min = 1),
                   checkboxInput("out_chr_mt", "Etiquetas chr con MT (1..22,X,Y,MT)", TRUE),
                   checkboxInput("use_rosetta", "Usar Rosetta (Mac ARM con PLINK x86_64)", FALSE),
                   div(
                     actionButton("run_plink_simple", "➋ Extraer VCF (PLINK)"),
                     downloadButton("dl_vcf", "Descargar selected_intervals.vcf.gz")
                   ),
                   tags$hr(),
                   
                   h3("Etapa 3 · Anotar con dbNSFP"),
                   textInput("dbnsfp_dir",  "Carpeta dbNSFP (contiene *.gz y el buscador)\n(Opcional: Añadir ruta a la carpeta dbNSFP)",
                             value = "/Volumes/LaCie_Drive/dbNSFP5.0a"),
                   textInput("dbnsfp_tool", "Buscador (clase Java)", value = "search_dbNSFP50a"),
                   numericInput("java_gb",   "Memoria Java (GB)", value = 6, min = 2, step = 1),
                   checkboxInput("show_dbnsfp_warning", "Mostrar aviso de duración", TRUE),
                   actionButton("run_dbnsfp", "➌ Ejecutar dbNSFP"),
                   downloadButton("dl_ann",  "Descargar dbNSFP.out"),
                   tags$hr(),
                   
                   h3("Etapa 3.5 · Formatear salida dbNSFP"),
                   actionButton("format_dbnsfp", "➍ Formatear salida dbNSFP"),
                   helpText("Convierte columnas numéricas y renombra cabeceras (si coinciden)."),
                   downloadButton("dl_dbnsfp_csv", "Descargar CSV normalizado"),
                   downloadButton("dl_dbnsfp_rds", "Descargar RDS normalizado"),
                   tags$hr(),
                   h4("Carpeta temporal"),
                   verbatimTextOutput("workdir_txt")
      ),
      mainPanel(
        h3("Manhattan (Etapa 1)"),
        plotlyOutput("manhattan", height = 520),
        tags$hr(),
        h4("SNPs ≥ umbral (selecciona filas)"),
        DT::DTOutput("hits_tbl"),
        tags$hr(),
        h4("Extraer vcf - Log PLINK"),
        verbatimTextOutput("plink_log"),
        tags$hr(),
        verbatimTextOutput("Anortar con dbNSFP50a - dbnsfp_log")
        
      )
    )
  ),
  
  # ================= TAB 2: Visualización (Etapa 4) =================
  tabPanel(
    "Visualización (Etapa 4)",
    sidebarLayout(
      sidebarPanel(width = 4,
                   h3("Fuente de datos (dbnsfp_normalized.csv)"),
                   radioButtons("viz_source", NULL,
                                choices = c("Subir archivo externo" = "upload",
                                            "Usar el último generado en la app" = "generated"),
                                selected = "upload"),
                   fileInput("viz_file", "Sube dbnsfp_normalized.csv",
                             accept = c(".csv"), buttonLabel = "Elegir CSV",
                             placeholder = "Sin archivo"),
                   tags$hr(),
                   checkboxInput("enable_viz", "Activar visualización", FALSE),
                   radioButtons("scoreorrank", "Tipo", choices = c("Score","Rankscore"),
                                selected = "Score", inline = TRUE),
                   numericInput("hm_topn", "Máx. SNPs (filas)", value = 100, min = 10, step = 10),
                   selectInput("hm_chr", "Cromosoma (filtro opcional)",
                               choices = c("Todos" = "all"), selected = "all"),
                   tags$hr(),
                   sliderInput("go_kegg_qcut", "FDR (q-value) máximo",
                               min = 0, max = 0.25, value = 0.05, step = 0.005),
                   tags$hr(),
                   downloadButton("dl_heatmap_png", "Descargar heatmap (PNG)")
      ),
      mainPanel(
        # Las pestañas de Etapa 4 se generan SOLO cuando enable_viz = TRUE y hay datos
        uiOutput("viz_tabs")
      )
    )
  )
)

# ===========================================================================  
# ---------- SERVER ----------
# ===========================================================================  

server <- function(input, output, session){
  
  source("R/format_dbnsfp.R", local = TRUE)
  
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(dplyr::select, dplyr::filter, dplyr::arrange,
                                 dplyr::mutate, base::setdiff, base::union, base::intersect)
  }
  
  # FDR seguro
  get_qcut <- function(x, default = 0.05) {
    if (!is.null(x) && length(x) == 1 && is.finite(x)) as.numeric(x) else default
  }
  plot_placeholder <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 5) +
      ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::theme_void()
  }
  
  
  # Carpeta de trabajo temporal para toda la sesión
  workdir <- file.path(tempdir(), paste0("gwas_inspector_", as.integer(Sys.time())))
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  output$workdir_txt <- renderText(workdir)
  
  # ---------- Etapa 1: Lectura & Manhattan ----------
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
    # posiciones Manhattan
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
      labs(x = "Cromosoma", y = "-log10(P)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text")
  })
  
  # Tabla de hits por umbral y selección de filas
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
  
  # --- helper: número PLINK -> etiqueta (1..22, X, Y, MT)
  chr_label_plink <- function(chr_num) {
    out <- as.character(chr_num)
    out[out == "23"] <- "X"
    out[out == "24"] <- "Y"
    out[out == "26"] <- "MT"
    out
  }
  
  # ---------- Etapa 1 -> generar selected_intervals (4 columnas) ----------
  selected_ranges_path <- reactiveVal(NULL)
  
  observeEvent(input$build_ranges, {
    h <- hits_df(); req(nrow(h) > 0)
    
    sel_idx <- input$hits_tbl_rows_selected
    picks <- if (length(sel_idx) > 0) {
      h[sel_idx, , drop = FALSE]
    } else {
      showNotification("No seleccionaste filas: usaré TODOS los SNPs ≥ umbral.", type = "warning")
      h
    }
    
    flank <- as.integer(if (!is.null(input$flank) && !is.na(input$flank)) input$flank else 0)
    
    # columnas seguras
    CHR <- as.integer(picks$CHR)
    BP  <- as.integer(picks$BP)
    
    ranges <- data.frame(
      CHR_lab = chr_label_plink(CHR),
      BP1     = pmax(1L, BP - flank),
      BP2     = BP + flank,
      stringsAsFactors = FALSE
    )
    
    # LABEL: rsID si hay; si no, "C{CHR}:{BP}±{flanco}"
    base_label <- ifelse(!is.na(picks$snp) & nzchar(picks$snp) & grepl("^rs", picks$snp, ignore.case = TRUE),
                         picks$snp,
                         paste0("C", ranges$CHR_lab, ":", BP))
    label <- if (flank > 0) paste0(base_label, "±", flank) else base_label
    
    # ordenar y deduplicar etiquetas si se repiten
    ord <- order(ranges$CHR_lab, ranges$BP1, ranges$BP2, na.last = NA)
    ranges <- ranges[ord, , drop = FALSE]
    label  <- label[ord]
    label  <- make.unique(label, sep = "_")
    
    # quitar filas con NA
    keep <- !is.na(ranges$CHR_lab) & !is.na(ranges$BP1) & !is.na(ranges$BP2)
    ranges <- ranges[keep, , drop = FALSE]
    label  <- label[keep]
    
    # construir líneas con ESPACIOS (no tabs) → "CHR BP1 BP2 LABEL"
    lines <- sprintf("%s %d %d %s", ranges$CHR_lab, ranges$BP1, ranges$BP2, label)
    lines <- unique(lines)
    lines <- lines[nchar(lines) > 0]
    
    # escribir archivo
    rng_path <- file.path(workdir, "selected_intervals.range")
    writeLines(lines, rng_path, useBytes = TRUE)
    
    # chequeo de las primeras líneas (≥4 tokens)
    head5 <- readLines(rng_path, n = 5, warn = FALSE)
    ok_tokens <- length(head5) == 0 || all(sapply(strsplit(head5, "\\s+"), function(x) length(x) >= 4))
    validate(need(ok_tokens, "El archivo de intervalos no cumple el formato CHR BP1 BP2 LABEL."))
    
    selected_ranges_path(rng_path)
    
    # Vista previa en UI + consola
    # Mostrar preview compacto en la UI (solo los 10 primeros) + resumen total
    output$ranges_preview <- renderText({
      n_total <- length(lines)
      n_show  <- min(10, n_total)
      head10  <- utils::head(lines, n_show)
      
      paste0(
        "Total de intervalos: ", n_total, "\n",
        if (n_total > 10) "Mostrando los 10 primeros:\n" else "Mostrando todos:\n",
        paste(head10, collapse = "\n")
      )
    })
  })
  
  ## ======================
  ##   ## Etapa 2 ====== 
  ## ======================
  
  ##  ====== Helpers específicos mac
  
  ensure_exec_mac <- function(path) {
    if (!file.exists(path)) return(FALSE)
    # permisos + quitar quarantine
    try(suppressWarnings(Sys.chmod(path, "0755")), silent = TRUE)
    if (Sys.info()[["sysname"]] == "Darwin") {
      try(system2("xattr", c("-d", "com.apple.quarantine", path),
                  stdout = FALSE, stderr = FALSE), silent = TRUE)
    }
    (.Platform$OS.type == "windows") || file.access(path, 1) == 0
  }
  
  # Devuelve comando + prefijo de argumentos si hace falta Rosetta
  plink_cmd_mac <- function(plink_path) {
    if (Sys.info()[["sysname"]] != "Darwin") return(c(normalizePath(plink_path)))
    os_arch <- tryCatch(system2("uname", "-m", stdout = TRUE), error = function(e) "x86_64")
    os_arch <- trimws(os_arch[1])
    finfo <- tryCatch(system2("file", normalizePath(plink_path), stdout = TRUE), error = function(e) "")
    is_x86_bin <- grepl("x86_64", paste(finfo, collapse = " "))
    is_arm_bin <- grepl("arm64",  paste(finfo, collapse = " "))
    if (identical(os_arch, "arm64") && is_x86_bin && !is_arm_bin && nzchar(Sys.which("arch"))) {
      return(c("arch", "-x86_64", normalizePath(plink_path)))
    }
    c(normalizePath(plink_path))
  }
  
  # Autodetecta bfiles completos (bed+bim+fam) en www/
  detect_bfiles <- function(dir = "www") {
    if (!dir.exists(dir)) return(character(0))
    fs <- list.files(dir, pattern = "\\.(bed|bim|fam)$", full.names = TRUE)
    if (!length(fs)) return(character(0))
    bases <- unique(sub("\\.(bed|bim|fam)$", "", fs))
    keep <- bases[sapply(bases, function(p) all(file.exists(paste0(p, c(".bed",".bim",".fam")))))]
    if (!length(keep)) return(character(0))
    ord <- order(file.info(paste0(keep, ".bed"))$size, decreasing = TRUE)
    keep[ord]
  }
  
  ## Autoselección de bfile en www si el campo está vacío
  observe({
    if (!nzchar(input$ref_bfile)) {
      cands <- detect_bfiles("www")
      if (length(cands)) updateTextInput(session, "ref_bfile", value = cands[1])
    }
  })
  
  # --- helpers mínimos ---
  resolve_in_www <- function(path_or_basename) {
    # si existe tal cual
    if (file.exists(path_or_basename)) return(normalizePath(path_or_basename, mustWork = TRUE))
    # prueba en ./www
    ww <- try(normalizePath("www", mustWork = TRUE), silent = TRUE)
    if (!inherits(ww, "try-error")) {
      cand <- file.path(ww, basename(path_or_basename))
      if (file.exists(cand)) return(normalizePath(cand, mustWork = TRUE))
    }
    NULL
  }
  
  detect_bfile_prefix <- function(dir = "www") {
    if (!dir.exists(dir)) return(character(0))
    fs <- list.files(dir, pattern = "\\.(bed|bim|fam)$", full.names = TRUE)
    if (!length(fs)) return(character(0))
    bases <- unique(sub("\\.(bed|bim|fam)$", "", fs))
    # elige el .bed más grande
    beds <- paste0(bases, ".bed")
    sizes <- suppressWarnings(file.info(beds)$size)
    bases[order(sizes, decreasing = TRUE)]
  }
  
  # --- estado compartido ---
  vcf_out_path <- reactiveVal(NULL)
  
  observeEvent(input$run_plink_simple, {
    # 1) paths
    rng_rel <- selected_ranges_path(); req(!is.null(rng_rel), file.exists(rng_rel))
    rng <- normalizePath(rng_rel, mustWork = TRUE)
    nr <- length(readLines(rng, warn = FALSE))
    validate(need(nr > 0, "selected_intervals.range está vacío."))
    
    # plink bin (NO forzamos normalizePath si te funcionaba con relativo)
    plink_bin_abs <- resolve_in_www(input$plink_path)
    validate(need(!is.null(plink_bin_abs), paste0("No encuentro el binario: ", input$plink_path)))
    # permisos + quitar quarantine por si acaso (mac)
    try(suppressWarnings(Sys.chmod(plink_bin_abs, "0755")), silent = TRUE)
    if (Sys.info()[["sysname"]] == "Darwin") {
      try(system2("xattr", c("-d", "com.apple.quarantine", plink_bin_abs),
                  stdout = FALSE, stderr = FALSE), silent = TRUE)
    }
    
    # bfile: usa el input si viene; si no, autodetecta en www
    ref_input <- input$ref_bfile
    pref <- NULL
    if (nzchar(ref_input)) {
      # acepta relativo (www/...) o absoluto; no normalizamos el prefijo (no existe como archivo)
      # sólo comprobaremos la existencia de .bed/.bim/.fam
      bed <- resolve_in_www(paste0(ref_input, ".bed"))
      bim <- resolve_in_www(paste0(ref_input, ".bim"))
      fam <- resolve_in_www(paste0(ref_input, ".fam"))
      validate(need(!is.null(bed) && !is.null(bim) && !is.null(fam),
                    "No encuentro los 3 ficheros bed/bim/fam con ese prefijo."))
      pref <- sub("\\.(bed|bim|fam)$", "", bed)  # prefijo ABS sin extensión
    } else {
      cands <- detect_bfile_prefix("www")
      validate(need(length(cands) > 0, "No hay BED/BIM/FAM en www/."))
      # primer candidato con los 3 archivos
      for (b in cands) {
        bed <- resolve_in_www(paste0(b, ".bed"))
        bim <- resolve_in_www(paste0(b, ".bim"))
        fam <- resolve_in_www(paste0(b, ".fam"))
        if (!is.null(bed) && !is.null(bim) && !is.null(fam)) {
          pref <- sub("\\.(bed|bim|fam)$", "", bed)
          break
        }
      }
      validate(need(!is.null(pref), "No se encontró un prefijo válido de bfile en www/."))
    }
    
    # 2) rutas de salida
    work <- normalizePath(workdir, mustWork = TRUE)
    candi_tmp  <- file.path(work, "candi_tmp")         # prefijo temporal BED/BIM/FAM
    out_prefix <- file.path(work, "selected_intervals") # prefijo VCF final
    vcf_gz     <- paste0(out_prefix, ".vcf.gz")
    
    # 3) construir args (como tu ejemplo; shQuote SOLO en valores)
    # Paso A: extraer por rango a BED
    argsA <- c(
      "--bfile", shQuote(pref),
      "--extract", "range", shQuote(rng),
      "--make-bed",
      "--out", shQuote(candi_tmp),
      "--allow-extra-chr",
      "--allow-no-sex",
      "--threads", as.character(input$threads19),
      "--chr-set", "26"
    )
    
    # Paso B: exportar a VCF bgzip
    argsB <- c(
      "--bfile", shQuote(candi_tmp),
      if (isTRUE(input$out_chr_mt)) c("--output-chr", "MT") else NULL,
      "--recode", "vcf", "bgz",
      "--out", shQuote(out_prefix),
      "--allow-extra-chr",
      "--allow-no-sex",
      "--threads", as.character(input$threads19)
    )
    
    # 4) Ejecutar (opcional Rosetta simple)
    run_cmd <- function(bin, args) {
      if (isTRUE(input$use_rosetta) &&
          Sys.info()[["sysname"]] == "Darwin" &&
          nzchar(Sys.which("arch"))) {
        system2("arch", args = c("-x86_64", bin, args), stdout = TRUE, stderr = TRUE)
      } else {
        system2(bin, args = args, stdout = TRUE, stderr = TRUE)
      }
    }
    
    # 5) Progreso + logs
    withProgress(message = "Ejecutando PLINK", value = 0, {
      setProgress(0.2, detail = "Paso A: --extract range → BED")
      outA <- run_cmd(plink_bin_abs, argsA)
      # muestra el tail del log A en UI
      output$plink_log <- renderText(paste(outA, collapse = "\n"))
      
      # comprueba que candi_tmp.bed existe
      validate(need(file.exists(paste0(candi_tmp, ".bed")), "PLINK fallo en el Paso A (no se creó candi_tmp.bed)."))
      
      setProgress(0.7, detail = "Paso B: --recode vcf bgz → VCF")
      outB <- run_cmd(plink_bin_abs, argsB)
      output$plink_log <- renderText(paste(c(outA, "\n---\n", outB), collapse = "\n"))
      
      validate(need(file.exists(vcf_gz), "PLINK terminó, pero no apareció el .vcf.gz"))
      
      vcf_out_path(vcf_gz)
      showNotification(paste0("VCF generado: ", basename(vcf_gz)), type = "message")
      setProgress(1)
    })
  })
  
  # botón de descarga (ya lo tenías)
  output$dl_vcf <- downloadHandler(
    filename = function() "selected_intervals.vcf.gz",
    content  = function(file) file.copy(vcf_out_path(), file, overwrite = TRUE)
  )
  
  ## ======================
  ##   Etapa 3 · dbNSFP
  ## ======================
  
  # ---- salida anotación
  ann_path <- reactiveVal(NULL)
  
  
  ensure_plain_vcf <- function(vcf_path, workdir) {
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
    validate(need(length(lines) > 0, "El VCF está vacío."))
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
  
  # ---- envío de email al terminar (opcional)
  send_completion_email <- function(to, subject = "dbNSFP: anotación completada", body = NULL, attach_path = NULL) {
    if (!isTRUE(input$email_notify) || !nzchar(to)) return(invisible(FALSE))
    if (!requireNamespace("blastula", quietly = TRUE)) {
      showNotification("No se pudo enviar email: paquete 'blastula' no instalado.", type = "warning")
      return(invisible(FALSE))
    }
    host <- Sys.getenv("SMTP_HOST")
    port <- as.integer(Sys.getenv("SMTP_PORT"))
    user <- Sys.getenv("SMTP_USER")
    pass <- Sys.getenv("SMTP_PASS")
    from <- Sys.getenv("SMTP_SENDER")
    if (!nzchar(host) || !nzchar(user) || !nzchar(pass) || !nzchar(from)) {
      showNotification("SMTP no configurado (variables de entorno). No se envía email.", type = "warning")
      return(invisible(FALSE))
    }
    
    body <- body %||% paste0(
      "Hola,\n\nTu anotación con dbNSFP ha finalizado.\n\n",
      if (!is.null(attach_path)) paste0("Adjunto: ", basename(attach_path), "\n") else "",
      "Carpeta temporal: ", workdir, "\n",
      "Fecha: ", format(Sys.time()), "\n\nSaludos."
    )
    email <- blastula::compose_email(
      body = blastula::md(paste0("**dbNSFP completado**\n\n", gsub("\n", "  \n", body)))
    )
    creds <- blastula::creds(user = user, pass = pass, host = host, port = port, use_ssl = FALSE, use_tls = TRUE)
    try({
      if (!is.null(attach_path) && file.exists(attach_path)) {
        email <- blastula::add_attachment(email, file = attach_path)
      }
      blastula::smtp_send(email = email, to = to, from = from, subject = subject, credentials = creds)
      showNotification(paste0("Email enviado a ", to), type = "message")
      TRUE
    }, silent = TRUE)
  }
  
  # ---- función que ejecuta la anotación (lo “gordo”)
  run_dbnsfp_real <- function() {
    req(vcf_out_path())
    vcf_in <- normalizePath(vcf_out_path(), mustWork = TRUE)
    
    db_dir  <- input$dbnsfp_dir
    tool_in <- input$dbnsfp_tool  # "search_dbNSFP50a" o ruta al .class
    
    validate(need(dir.exists(db_dir), paste0("No existe la carpeta: ", db_dir)))
    chr_files <- list.files(db_dir, pattern="^dbNSFP5\\.0a_variant\\.chr(\\d+|[XYM])\\.gz$", full.names = TRUE)
    validate(need(length(chr_files) >= 24 && file.exists(file.path(db_dir, "dbNSFP5.0_gene.gz")),
                  "Faltan ficheros en la base (chr*.gz / dbNSFP5.0_gene.gz)."))
    
    java_mem <- as.integer(input$java_gb %||% 6)
    main_class <- basename(tool_in)
    main_class <- sub("\\.class$", "", main_class, ignore.case = TRUE)
    
    # preparar VCF v4.0
    vcf_plain <- ensure_plain_vcf(vcf_in, workdir)
    vcf_v40   <- file.path(workdir, "selected_intervals.v40.vcf")
    vcf_v40   <- force_vcf_version(vcf_plain, vcf_v40, "VCFv4.0")
    
    # contar variantes
    n_var <- sum(!grepl("^#", readLines(vcf_v40, warn = FALSE)))
    validate(need(n_var > 0, "El VCF (v4.0) no contiene variantes que anotar."))
    
    out_file <- file.path(workdir, "dbnsfp.out")
    cmd <- list(
      cmd  = "java",
      args = c(paste0("-Xmx", java_mem, "g"),
               main_class,
               "-i", normalizePath(vcf_v40, mustWork = TRUE),
               "-o", normalizePath(out_file, mustWork = FALSE))
    )
    
    withProgress(message = "Ejecutando dbNSFP", value = 0.5, {
      res <- tryCatch({
        old <- getwd(); on.exit(setwd(old), add = TRUE)
        setwd(db_dir)  # ejecutar en la carpeta de la base
        system2(cmd$cmd, args = cmd$args, stdout = TRUE, stderr = TRUE, wait = TRUE)
      }, error = function(e) { attr(e, "status") <- 1L; e })
      
      status <- attr(res, "status"); if (is.null(status)) status <- 0L
      out_exists <- file.exists(out_file) && file.info(out_file)$size > 0
      
      header <- paste0(
        "DB dir: ", db_dir, "\n",
        "Tool:   ", main_class, "  (ejecución en wd=db_dir)\n",
        "VCF:    ", vcf_v40, "\n\n",
        "CMD:\n", cmd$cmd, " ", paste(cmd$args, collapse = " "), "\n\n",
        "exit=", status, "\n\n"
      )
      output$dbnsfp_log <- renderText(paste0(
        header,
        paste(res, collapse = "\n"),
        if (out_exists) paste0("\n\n✔ Anotación completada.\nPreview (20 líneas):\n",
                               paste(utils::head(readLines(out_file, warn = FALSE), 20), collapse = "\n"))
        else "\n✖ No se generó la salida."
      ))
      
      validate(need(out_exists, "dbNSFP no generó la salida; revisa el log de arriba."))
      ann_path(out_file)
      showNotification(paste0("dbNSFP OK: ", basename(out_file)), type = "message")
      setProgress(1)
    })
    
    # Email opcional
    if (isTRUE(input$email_notify) && nzchar(input$notify_email)) {
      send_completion_email(
        to = input$notify_email,
        body = paste0(
          "El archivo de salida es: ", basename(ann_path()), "\n",
          "Vista previa (primeras 20 líneas):\n\n",
          paste(utils::head(readLines(ann_path(), warn = FALSE), 20), collapse = "\n")
        ),
        attach_path = ann_path()
      )
    }
  }
  
  # ---- aviso/confirmación previa
  observeEvent(input$run_dbnsfp, {
    if (isTRUE(input$show_dbnsfp_warning)) {
      showModal(modalDialog(
        title = "Este proceso puede tardar (minutos u horas)",
        tagList(
          p("La anotación con dbNSFP puede ser lenta, especialmente desde disco externo."),
          p("Puedes dejar la app trabajando; verás el log de progreso aquí.")
        ),
        footer = tagList(
          modalButton("Cancelar"),
          actionButton("confirm_run_dbnsfp", "Entendido, continuar", class = "btn-primary")
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
  
  # ---- descarga
  output$dl_ann <- downloadHandler(
    filename = function() "dbnsfp_annotation.out",
    content  = function(file) file.copy(ann_path(), file, overwrite = TRUE)
  )

  #  ============= Etapa 3.5
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
      paste0("dbNSFP normalizado: ", basename(res$csv), " (", res$n_rows, " x ", res$n_cols, ")"),
      type = "message"
    )
  })
  
 
  # 4.0) Dataset normalizado para la Visualización (independiente)
  dbnsfp_norm_df <- reactive({
    src <- input$viz_source %||% "upload"
    
    # a) Modo subir archivo
    if (identical(src, "upload")) {
      req(input$viz_file)
      f <- input$viz_file$datapath
      validate(need(file.exists(f), "El archivo subido no existe."))
      df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
      validate(need(!is.null(df) && nrow(df) > 0, "No pude leer el CSV subido."))
      attr(df, "source") <- "upload"
      return(df)
    }
    
    # b) Modo usar el generado por la app (Etapa 3.5)
    if (identical(src, "generated")) {
      path_csv <- dbnsfp_norm_path_csv()
      validate(need(!is.null(path_csv) && file.exists(path_csv),
                    "No existe un dbnsfp_normalized.csv generado en esta sesión."))
      df <- tryCatch(read.csv(path_csv, check.names = FALSE), error = function(e) NULL)
      validate(need(!is.null(df) && nrow(df) > 0, "El CSV generado está vacío o ilegible."))
      attr(df, "source") <- "generated"
      return(df)
    }
    
    # fallback
    validate(need(FALSE, "Selecciona una fuente de datos para la visualización."))
  })
  
  # === Base visualización: aplica filtro chr y garantiza rsid_final ===
  viz_base_df <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    snp_col <- session$userData$snp_col %||% "SNP"
    chr_col <- session$userData$chr_col
    bp_col  <- session$userData$bp_col
    
    # normalizador chr: "chr1" -> "1"
    norm_chr <- function(x) sub("^chr", "", tolower(trimws(as.character(x))))
    
    # filtro por cromosoma (si procede)
    if (!is.null(chr_col) && chr_col %in% names(df) &&
        !is.null(input$hm_chr) && input$hm_chr != "all") {
      key <- norm_chr(input$hm_chr)
      df  <- df[norm_chr(df[[chr_col]]) == key, , drop = FALSE]
    }
    
    # construir columna rsid_final 100% presente
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
    
    # anota cuál es la columna efectiva de ID para el resto de módulos
    attr(df, "snp_col_effective") <- "rsid_final"
    attr(df, "chr_col") <- chr_col
    attr(df, "bp_col")  <- bp_col
    df
  })

  # cuando haya datos normalizados disponibles para la visualización
  observeEvent(dbnsfp_norm_df(), {
    df <- dbnsfp_norm_df(); req(nrow(df) > 0)
    
    # Usa la columna de cromosoma detectada previamente
    chr_col <- session$userData$chr_col
    if (is.null(chr_col) || !(chr_col %in% names(df))) {
      updateSelectInput(session, "hm_chr", choices = c("Todos" = "all"), selected = "all")
      return()
    }
    
    # Normaliza etiquetas tipo "chr1" -> "1"
    norm_chr <- function(x) sub("^chr", "", tolower(trimws(as.character(x))))
    vals <- unique(df[[chr_col]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    vals_n <- norm_chr(vals)
    
    # Presenta como etiquetas "chrN" pero guarda el valor normalizado
    # separa numéricos y especiales para ordenar mejor
    is_num <- suppressWarnings(!is.na(as.integer(vals_n)))
    nums   <- sort(unique(as.integer(vals_n[is_num])))
    specs  <- setdiff(unique(toupper(vals_n[!is_num])), character(0))
    # prioriza X,Y,MT si existen
    spec_order <- c("X","Y","MT","M","MTDNA")
    specs <- c(intersect(spec_order, specs), setdiff(specs, spec_order))
    
    choices <- character(0)
    if (length(nums))  choices <- c(choices, stats::setNames(as.character(nums), paste0("chr", nums)))
    if (length(specs)) choices <- c(choices, stats::setNames(specs, paste0("chr", specs)))
    
    if (!length(choices)) {
      updateSelectInput(session, "hm_chr", choices = c("Todos"="all"), selected="all")
    } else {
      updateSelectInput(session, "hm_chr",
                        choices = c("Todos" = "all", choices),
                        selected = "all")
    }
  }, ignoreInit = FALSE)
  
  # =========================
  #  Etapa 4 · Visualización
  # =========================
  
  viz_df_go <- reactive({
    df <- dbnsfp_norm_df()
    req(df, nrow(df) > 0)
    
    # Determinar columna chr disponible
    chr_col <- if ("chr" %in% names(df)) "chr" else if ("CHR" %in% names(df)) "CHR" else NULL
    
    # Filtro de cromosoma, si procede
    if (!is.null(chr_col) && !is.null(input$hm_chr) && !identical(input$hm_chr, "all")) {
      df <- df[as.character(df[[chr_col]]) == as.character(input$hm_chr), , drop = FALSE]
    }
    
    # Orden chr/BP si existen
    bp_col <- if ("BP" %in% names(df)) "BP" else if ("bp" %in% names(df)) "bp" else NULL
    if (!is.null(chr_col) && !is.null(bp_col)) {
      BPnum <- suppressWarnings(as.numeric(df[[bp_col]]))
      ord <- order(as.character(df[[chr_col]]), BPnum, na.last = TRUE)
      df <- df[ord, , drop = FALSE]
    }
    df
  })
  
  # df para visualización (heatmap/GO/KEGG), respeta hm_chr y hm_topn
  viz_df <- reactive({
    df <- dbnsfp_norm_df(); req(df, nrow(df) > 0)
    
    chr_col <- session$userData$chr_col
    bp_col  <- session$userData$bp_col %||% if ("BP" %in% names(df)) "BP" else NULL
    
    # Filtro por cromosoma si procede
    if (!is.null(chr_col) && chr_col %in% names(df) &&
        !is.null(input$hm_chr) && !identical(input$hm_chr, "all")) {
      norm_chr <- function(x) sub("^chr", "", tolower(trimws(as.character(x))))
      key <- sub("^chr", "", tolower(trimws(as.character(input$hm_chr))))
      df  <- df[norm_chr(df[[chr_col]]) == key, , drop = FALSE]
    }
    
    # Orden por chr/BP si existen
    if (!is.null(chr_col) && !is.null(bp_col) &&
        chr_col %in% names(df) && bp_col %in% names(df)) {
      BPnum <- suppressWarnings(as.numeric(df[[bp_col]]))
      ord <- order(as.character(df[[chr_col]]), BPnum, na.last = TRUE)
      df <- df[ord, , drop = FALSE]
    }
    
    # top-N para heatmap (si quieres sólo aplicarlo al heatmap, puedes quitarlo de aquí)
    topn <- if (!is.null(input$hm_topn) && is.finite(input$hm_topn)) input$hm_topn else 100
    if (nrow(df) > topn) df <- utils::head(df, topn)
    
    df
  })
  
  
  
  # ===== Heatmap data (equiv. score_input) =====
  heatmap_df <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    # columnas clave (intenta respetar tus nombres normalizados)
    snp_col <- session$userData$snp_col %||% "SNP"
    chr_col <- session$userData$chr_col %||% "chr"
    bp_col  <- session$userData$bp_col  %||% "BP"
    
    # candidatos por defecto (tus métricas típicas de ejemplo)
    default_scores <- c(
      "SIFT_score","SIFT4G_score",
      "Polyphen2_HDIV_score","Polyphen2_HVAR_score",
      "MutationTaster_score","MutationAssessor_score",
      "PROVEAN_score","REVEL_score",
      "GERP_91_mammals","phyloP17way_primate","phastCons17way_primate"
    )
    default_ranks  <- c(
      "SIFT_converted_rankscore","SIFT4G_converted_rankscore",
      "Polyphen2_HDIV_rankscore","Polyphen2_HVAR_rankscore",
      "MutationAssessor_rankscore","PROVEAN_converted_rankscore",
      "REVEL_rankscore","GERP_91_mammals_rankscore",
      "phyloP17way_primate_rankscore","phastCons17way_primate_rankscore"
    )
    
    # columnas elegidas por UI (o fallback si vienen vacías)
    if (identical(input$hm_mode, "Rankscore")) {
      cols_metrics <- input$hm_ranks
      if (length(cols_metrics) < 2) cols_metrics <- base::intersect(default_ranks, names(df))
    } else {
      cols_metrics <- input$hm_scores
      if (length(cols_metrics) < 2) cols_metrics <- base::intersect(default_scores, names(df))
    }
    validate(need(length(cols_metrics) >= 2,
                  "No hay al menos 2 columnas de métricas disponibles para el heatmap."))
    
    # columnas base presentes
    cols_base <- base::intersect(c(snp_col, chr_col, bp_col, "ref","alt","genename"), names(df))
    out <- dplyr::bind_cols(df[, cols_base, drop = FALSE], df[, cols_metrics, drop = FALSE])
    
    # filtro cromosoma (si procede)
    if (!is.null(input$hm_chr) && input$hm_chr != "all" && chr_col %in% names(out)) {
      key <- suppressWarnings(as.numeric(input$hm_chr))
      key <- if (!is.na(key)) key else input$hm_chr
      out <- out[out[[chr_col]] == key, , drop = FALSE]
    }
    
    # ordenar y limitar top-N
    if (chr_col %in% names(out) && bp_col %in% names(out)) {
      out <- out %>% dplyr::arrange(.data[[chr_col]], .data[[bp_col]])
    }
    n_show <- input$hm_topn %||% 100
    out <- utils::head(out, n_show)
    
    # renombra columnas métricas con salto de línea para estética
    new_names <- names(out)
    met_idx <- which(names(out) %in% cols_metrics)
    new_names[met_idx] <- gsub("_", "\n", new_names[met_idx])
    names(out) <- new_names
    
    attr(out, "snp_col") <- snp_col
    attr(out, "chr_col") <- chr_col
    attr(out, "bp_col")  <- bp_col
    out
  })
  
  ################
  # Debug del esquema detectado (opcional)
  output$dbnsfp_schema <- renderText({ "" })
  
  # Detección de columnas (idéntico a lo que ya te pasé, pero leyendo de dbnsfp_norm_df())
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
    
    # Métricas
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
    
    # guardar
    session$userData$snp_col    <- snp_col
    session$userData$bp_col     <- bp_col
    session$userData$chr_col    <- chr_col
    session$userData$gene_col   <- gene_col
    session$userData$score_cols <- sort(unique(score_cols))
    session$userData$rank_cols  <- sort(unique(rank_cols))
    session$userData$pred_cols  <- sort(unique(pred_cols))
    session$userData$go_cols    <- sort(unique(go_cols))
    session$userData$kegg_cols  <- sort(unique(kegg_cols))
    
    output$dbnsfp_schema <- renderText(paste(
      sprintf("Fuente: %s", attr(df, "source") %||% "desconocida"),
      sprintf("Filas: %s, Columnas: %s", nrow(df), ncol(df)),
      paste("SNP:",  ifelse(is.null(snp_col),"—", snp_col)),
      paste("CHR:",  ifelse(is.null(chr_col),"—", chr_col)),
      paste("BP :",  ifelse(is.null(bp_col) ,"—", bp_col)),
      paste("GENE:", ifelse(is.null(gene_col),"—", gene_col)),
      sep = "\n"
    ))
  })
  
  # 4.2) Mostrar controles solo si hay dataset
  output$viz_controls <- renderUI({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    
    score_cols <- session$userData$score_cols %||% character(0)
    rank_cols  <- session$userData$rank_cols  %||% character(0)
    chr_col    <- session$userData$chr_col
    
    # Construir selector de cromosoma solo si existe la columna
    chr_ui <- NULL
    if (!is.null(chr_col) && chr_col %in% names(df)) {
      chr_vals <- unique(df[[chr_col]])
      if (length(chr_vals) > 0) {
        if (is.numeric(chr_vals)) {
          labs <- paste0("chr", chr_vals)
          choices <- c("Todos" = "all", stats::setNames(as.character(chr_vals), labs))
        } else {
          u <- as.character(chr_vals)
          choices <- c("Todos" = "all", stats::setNames(u, u))
        }
        chr_ui <- selectInput("hm_chr", "Cromosoma (filtro opcional)", choices = choices, selected = "all")
      } else {
        chr_ui <- tags$em("Columna de cromosoma vacía; se omite filtro.")
      }
    } else {
      chr_ui <- tags$em("No se encontró columna de cromosoma; no se aplicará filtro por chr.")
    }
    
    tagList(
      tags$div(style="margin-bottom:8px; font-size: 12px; white-space: pre-wrap; color:#666; border:1px dashed #ccc; padding:6px;",
               strong("Esquema detectado"), tags$br(), verbatimTextOutput("dbnsfp_schema")),
      h4("Heatmap"),
     # radioButtons("hm_mode", "Tipo de métrica", choices = c("Score", "Rankscore"), inline = TRUE),
      conditionalPanel(
        "input.hm_mode == 'Score'",
        selectizeInput("hm_scores", "Selecciona scores (2–15)", choices = sort(score_cols),
                       multiple = TRUE, options = list(maxItems = 15, placeholder = "elige columnas score…"))
      ),
      conditionalPanel(
        "input.hm_mode == 'Rankscore'",
        selectizeInput("hm_ranks", "Selecciona rank-scores (2–15)", choices = sort(union(rank_cols, character(0))),
                       multiple = TRUE, options = list(maxItems = 15, placeholder = "elige columnas rankscore…"))
      ),
      chr_ui,
      numericInput("hm_topn", "Máx. filas (SNPs) a mostrar", value = 100, min = 10, max = 5000, step = 10),
      checkboxInput("hm_scale_cols", "Escalar 0–1 por columna", TRUE),
      
      tags$hr(),
      h4("GO / KEGG"),
      sliderInput("go_kegg_pcut", "Umbral FDR (BH)", min = 0, max = 0.25, value = 0.05, step = 0.005),
      numericInput("kegg_top", "Top KEGG a mostrar", value = 15, min = 5, max = 50, step = 1)
    )
  })
  
  
  # 4.3) Tabs de visualización
  output$viz_tabs <- renderUI({
    req(input$enable_viz)
    req(dbnsfp_norm_df())
    
    tagList(
      h3("Visualización"),
      p(em(sprintf("CSV: %s", attr(dbnsfp_norm_df(), "source") %||% "—"))),
      tabsetPanel(
        id = "viz_tabs",
        tabPanel("Heatmap", uiOutput("dynamic_heatmap_plot"), DT::DTOutput("heatmap_table")),
        tabPanel("GO (BP/CC/MF)", plotOutput("go_bar")),
        tabPanel("KEGG", plotOutput("kegg_bar"), tags$hr(), h4("Tabla KEGG"), DT::DTOutput("kegg_table"))
      )
    )
  })
  
  
  # ==== Fuente del heatmap (equivale a tu score_input()) ====
  score_input <- reactive({
    df <- viz_base_df()
    validate(need(nrow(df) > 0, "No hay datos tras el filtro de cromosoma."))
    
    snp_eff <- attr(df, "snp_col_effective") %||% "rsid_final"
    chr_col <- attr(df, "chr_col")
    bp_col  <- attr(df, "bp_col")
    
    wanted <- c(
      snp_eff, chr_col, bp_col, "ref", "alt", "genename",
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
    validate(need(length(keep) >= 5, "El CSV no tiene suficientes columnas esperadas para el heatmap."))
    
    out <- df[, keep, drop = FALSE]
    
    # ordenar y truncar
    if (!is.null(chr_col) && !is.null(bp_col) && chr_col %in% names(out) && bp_col %in% names(out)) {
      out <- dplyr::arrange(out, .data[[chr_col]], .data[[bp_col]])
    }
    n_show <- input$hm_topn %||% 100
    out <- utils::head(out, n_show)
    
    attr(out, "snp_col") <- snp_eff
    attr(out, "chr_col") <- chr_col
    attr(out, "bp_col")  <- bp_col
    out
  })
  
  # ====== Núcleo común para Heatmap (reutilizable) ======
  heatmap_data <- reactive({
    df <- dbnsfp_norm_df(); req(nrow(df) > 0)
    
    chr_col <- session$userData$chr_col
    bp_col  <- session$userData$bp_col %||% if ("BP" %in% names(df)) "BP" else NULL
    
    # Filtro por cromosoma si procede
    if (!is.null(chr_col) && chr_col %in% names(df) &&
        !is.null(input$hm_chr) && !identical(input$hm_chr, "all")) {
      norm_chr <- function(x) sub("^chr", "", tolower(trimws(as.character(x))))
      key <- sub("^chr", "", tolower(trimws(as.character(input$hm_chr))))
      df  <- df[norm_chr(df[[chr_col]]) == key, , drop = FALSE]
    }
    validate(need(nrow(df) > 0, "No hay filas tras filtrar por cromosoma."))
    
    # Orden estable por chr/BP si existen
    if (!is.null(chr_col) && !is.null(bp_col) &&
        chr_col %in% names(df) && bp_col %in% names(df)) {
      BPnum <- suppressWarnings(as.numeric(df[[bp_col]]))
      ord <- order(as.character(df[[chr_col]]), BPnum, na.last = TRUE)
      df <- df[ord, , drop = FALSE]
    }
    
    # SNP preferente
    snp_col <- if ("rsid_final" %in% names(df)) "rsid_final" else if ("SNP" %in% names(df)) "SNP" else NULL
    
    nms <- names(df)
    if (identical(input$scoreorrank, "Score")) {
      metric_cols <- c(
        grep("_score$", nms, value = TRUE),
        intersect(c("GERP_91_mammals","phyloP17way_primate","phastCons17way_primate"), nms)
      )
      title_txt <- "Heatmap Scores"
    } else {
      metric_cols <- c(
        grep("_rankscore$", nms, value = TRUE),
        intersect("GERP++_RS_rankscore", nms)
      )
      title_txt <- "Heatmap Rank-scores"
    }
    
    id_cols <- intersect(unique(c(snp_col, "rsid_final", chr_col, bp_col, "ref","alt","genename")), nms)
    metric_cols <- setdiff(unique(metric_cols), id_cols)
    validate(need(length(metric_cols) >= 2, "Elige al menos 2 columnas para el heatmap."))
    
    topn <- if (!is.null(input$hm_topn) && is.finite(input$hm_topn)) as.integer(input$hm_topn) else 100
    if (nrow(df) > topn) df <- head(df, topn)
    
    df_hm <- df[, c(id_cols, metric_cols), drop = FALSE]
    for (cc in metric_cols) df_hm[[cc]] <- suppressWarnings(as.numeric(df_hm[[cc]]))
    
    mat <- as.matrix(df_hm[, metric_cols, drop = FALSE])
    
    rs <- if ("rsid_final" %in% names(df_hm)) {
      df_hm$rsid_final
    } else if (!is.null(snp_col) && snp_col %in% names(df_hm)) {
      df_hm[[snp_col]]
    } else if (!is.null(chr_col) && !is.null(bp_col) &&
               chr_col %in% names(df_hm) && bp_col %in% names(df_hm)) {
      paste0(df_hm[[chr_col]], ":", df_hm[[bp_col]])
    } else {
      paste0("var_", seq_len(nrow(df_hm)))
    }
    rownames(mat) <- as.character(rs)
    
    mat_norm <- apply(mat, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || (rng[2] - rng[1] == 0)) return(rep(0, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    })
    colnames(mat_norm) <- gsub("_", "\n", colnames(mat_norm), fixed = TRUE)
    
    show_numbers <- nrow(mat) <= 60 && ncol(mat) <= 20
    chr_tag <- if (!is.null(input$hm_chr) && !identical(input$hm_chr, "all")) paste0(" — chr ", input$hm_chr) else ""
    list(mat = mat, mat_norm = mat_norm, rs = rs, title = paste0(title_txt, chr_tag), show_numbers = show_numbers)
  })
  
  
  
  # ===== Heatmap plot =====
  output$heatmap_plot <- renderPlot({
    hd <- heatmap_data()  # <- usa el reactive de arriba
    validate(
      need(!is.null(hd$mat_norm) && nrow(hd$mat_norm) > 0 && ncol(hd$mat_norm) > 0,
           "No hay datos para el heatmap.")
    )
    
    heatmap_plot <- ggplotify::as.ggplot(function() {
      pheatmap::pheatmap(
        hd$mat_norm,
        cluster_rows   = FALSE,
        cluster_cols   = FALSE,
        display_numbers = if (hd$show_numbers) round(hd$mat, 2) else FALSE,
        labels_row     = hd$rs,
        show_rownames  = TRUE,
        main           = hd$title,
        fontsize_number = 10,
        fontsize_row    = 8,
        fontsize_col    = 10
      )
    })
    print(heatmap_plot)
  })
  
  
  output$heatmap_plotXXXXX <- renderPlot({
    req(score_input())
    df <- score_input()
    validate(need(nrow(df) > 0, "⚠️ score_input es NULL o vacío"))
    
    # 1) columnas efectivas (siempre preferimos 'rsid_final' si existe)
    snp_col <- if ("rsid_final" %in% names(df)) "rsid_final" else (attr(df, "snp_col") %||% "SNP")
    chr_col <- attr(df, "chr_col") %||% "chr"
    bp_col  <- attr(df, "bp_col")  %||% "BP"
    
    # 2) métricas según modo
    nms <- names(df)
    if (identical(input$scoreorrank, "Score")) {
      metric_cols <- c(
        nms[grepl("_score$", nms, perl = TRUE)],
        base::intersect(c("GERP_91_mammals","phyloP17way_primate","phastCons17way_primate"), nms)
      )
    } else {
      metric_cols <- c(
        nms[grepl("_rankscore$", nms, perl = TRUE)],
        base::intersect("GERP++_RS_rankscore", nms)
      )
    }
    
    # 3) columnas identificadoras (incluye rsid_final explícitamente si existe)
    id_cols <- base::intersect(
      unique(c(snp_col, "rsid_final", chr_col, bp_col, "ref", "alt", "genename")),
      nms
    )
    
    # 4) evitar duplicados/colisiones
    metric_cols <- setdiff(unique(metric_cols), id_cols)
    cols_to_take <- unique(c(id_cols, metric_cols))
    validate(need(length(metric_cols) >= 2, "Elige al menos 2 columnas para el heatmap."))
    
    # 5) construir df para heatmap
    df_hm <- df[, cols_to_take, drop = FALSE]
    
    # 6) renombrar SOLO las columnas de métricas (nunca las ID)
    new_names <- names(df_hm)
    m_idx <- match(metric_cols, new_names)
    new_names[m_idx] <- gsub("_", "\n", new_names[m_idx], fixed = TRUE)
    names(df_hm) <- new_names
    metric_cols_renamed <- names(df_hm)[m_idx]
    
    # 7) convertir métricas a numéricas
    df_hm[metric_cols_renamed] <- lapply(
      df_hm[metric_cols_renamed],
      function(x) suppressWarnings(as.numeric(as.character(x)))
    )
    
    # 8) matriz y etiquetas de fila (rsid_final siempre que esté)
    score_matrix <- as.matrix(df_hm[, metric_cols_renamed, drop = FALSE])
    
    rs_vec <- if ("rsid_final" %in% names(df_hm)) df_hm[["rsid_final"]] else {
      if (snp_col %in% names(df_hm)) df_hm[[snp_col]] else NA
    }
    # último fallback: construir chr:BP si rs_vec viene vacío
    if (is.null(rs_vec) || all(is.na(rs_vec) | rs_vec == "")) {
      if (all(c(chr_col, bp_col) %in% names(df_hm))) {
        rs_vec <- paste0(df_hm[[chr_col]], ":", df_hm[[bp_col]])
      } else {
        rs_vec <- paste0("var_", seq_len(nrow(df_hm)))
      }
    }
    rownames(score_matrix) <- as.character(rs_vec)
    
    # 9) normalización protegida
    score_matrix_normalized <- apply(score_matrix, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || (rng[2] - rng[1] == 0)) return(rep(0, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    })
    
    # 10) si es grande, no numera celdas
    show_numbers <- nrow(score_matrix) <= 60 && ncol(score_matrix) <= 20
    
    heatmap_plot <- ggplotify::as.ggplot(function() {
      pheatmap::pheatmap(
        score_matrix_normalized,
        color = colorRampPalette(c("#003366", "#66CCFF", "#F7F7F7", "#FF9966", "#660000"))(100),
        
        cluster_rows = FALSE, cluster_cols = FALSE,
        display_numbers = if (show_numbers) score_matrix else FALSE,
        labels_row = rs_vec,
        show_rownames = TRUE,
        number_format = "%.2f",
        main = if (identical(input$scoreorrank, "Score")) "Heatmap Scores" else "Heatmap Rank-scores",
        fontsize_number = 10, fontsize_row = 8, fontsize_col = 10
      )
    })
    print(heatmap_plot)
  })
  
  # ===== Heatmap dynamic height wrapper =====
  output$dynamic_heatmap_plot <- renderUI({
    df <- score_input()
    n_rows <- if (!is.null(df)) nrow(df) else 0
    plot_height <- max(400, n_rows * 20)  # 20 px por fila
    plotOutput("heatmap_plot", height = paste0(plot_height, "px"))
  })
  
  output$dl_heatmap_png <- downloadHandler(
    filename = function() sprintf("heatmap_%s.png", tolower(input$scoreorrank)),
    contentType = "image/png",
    content = function(file) {
      # Congelar reactividad y recuperar datos listos para pintar
      hd <- tryCatch(isolate(heatmap_data()), error = function(e) NULL)
      
      # Si no hay datos válidos, genera un PNG con mensaje (no un 404)
      if (is.null(hd) || is.null(hd$mat) || nrow(hd$mat) == 0 || ncol(hd$mat) == 0) {
        png(file, width = 1200L, height = 400L, res = 120)
        par(mar = c(0,0,0,0)); plot.new()
        text(0.5, 0.5, "Sin datos para descargar el heatmap.", cex = 1.4)
        dev.off()
        return(invisible())
      }
      
      # Dimensiones en pulgadas (pheatmap usa inches cuando se pasa 'filename')
      # Ajusta si quieres celdas más grandes/pequeñas.
      height_in <- max(6, nrow(hd$mat) * 0.20)  # ~0.20" por fila, mínimo 6"
      width_in  <- max(9, ncol(hd$mat) * 0.50)  # ~0.50" por columna, mínimo 9"
      
      # Escribir directamente con pheatmap al fichero de salida
      pheatmap::pheatmap(
        hd$mat_norm,
        cluster_rows = FALSE, cluster_cols = FALSE,
        display_numbers = if (hd$show_numbers) round(hd$mat, 2) else FALSE,
        labels_row = hd$rs,
        show_rownames = TRUE,
        main = hd$title,
        fontsize_number = 10, fontsize_row = 8, fontsize_col = 10,
        filename = file,       # <- clave: pheatmap escribe el PNG
        width = width_in,      # pulgadas
        height = height_in     # pulgadas
      )
    }
  )
  
  
  output$heatmap_table <- DT::renderDT({
    df <- score_input()
    id <- if ("rsid_final" %in% names(df)) "rsid_final" else attr(df, "snp_col") %||% "SNP"
    chr <- attr(df, "chr_col") %||% "chr"
    bp  <- attr(df, "bp_col")  %||% "BP"
    show_df <- head(df[, c(id, chr, bp, setdiff(names(df), c(id, chr, bp))), drop = FALSE], 200)
    DT::datatable(
      show_df,
      options = list(pageLength = 10, scrollX = TRUE),
      escape = FALSE
    )
  })
  
  # ==========================
  #  GO Enrichment (Etapa 4)
  # ==========================
  output$go_bar <- renderPlot({
    req(input$enable_viz)
    
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(plot_placeholder("Instala 'clusterProfiler'"))
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))    return(plot_placeholder("Instala 'org.Hs.eg.db'"))
    
    df <- viz_df(); req(nrow(df) > 0)
    
    # Columna de genes
    gene_col <- session$userData$gene_col
    if (is.null(gene_col) || !(gene_col %in% names(df))) return(plot_placeholder("No existe columna 'genename'"))
    
    # ⚠️ Split por ';' (a diferencia de antes)
    genes <- df[[gene_col]]
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genes <- unique(trimws(unlist(strsplit(as.character(genes), ";", fixed = TRUE))))
    if (!length(genes)) return(plot_placeholder("Sin genes válidos en 'genename'"))
    
    suppressMessages(library(clusterProfiler))
    suppressMessages(library(org.Hs.eg.db))
  conv <- try(clusterProfiler::bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db), silent = TRUE)
  if (inherits(conv, "try-error") || is.null(conv) || nrow(conv) == 0) return(plot_placeholder("No se pudieron mapear SYMBOL→ENTREZID"))
  
  entrez <- unique(na.omit(as.character(conv$ENTREZID)))
  if (!length(entrez)) return(plot_placeholder("Mapeo a ENTREZID vacío"))
  
  qcut <- get_qcut(input$go_kegg_qcut, 0.05)
  ont  <- "BP"  # si quieres UI para BP/CC/MF, añade un radioButtons y úsalo aquí
  
  ego <- try(clusterProfiler::enrichGO(
    gene = entrez, OrgDb = org.Hs.eg.db, ont = ont, keyType = "ENTREZID",
    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE
  ), silent = TRUE)
  if (inherits(ego, "try-error") || is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(plot_placeholder(sprintf("Sin términos GO (%s)", ont)))
  }
  
  # aplica FDR del slider
  res <- ego@result
  res <- res[!is.na(res$p.adjust) & res$p.adjust <= qcut, , drop = FALSE]
  if (!nrow(res)) return(plot_placeholder(sprintf("Sin GO (%s) con FDR ≤ %.3f", ont, qcut)))
  
  ego2 <- ego; ego2@result <- res
  clusterProfiler::barplot(ego2, showCategory = min(20, nrow(res)))
  })

  
  
  # 4.6) KEGG enrichment con clusterProfiler
  # Necesita: org.Hs.eg.db, clusterProfiler
  kegg_result <- reactive({
    req(dbnsfp_norm_df())
    df <- dbnsfp_norm_df()
    gene_col <- session$userData$gene_col
    validate(need(!is.null(gene_col) && gene_col %in% names(df),
                  "No hay columna de símbolos de gen (genename)."))
    
    # desdoblar genename (separado por ;)
    genes <- df[[gene_col]]
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(NULL)
    genes <- unique(unlist(strsplit(genes, ";", fixed = TRUE)))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    validate(need(length(genes) > 0, "No hay símbolos de gen válidos."))
    
    # mapear a ENTREZ
    suppressWarnings({
      gmap <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db)
    })
    validate(need(!is.null(gmap) && nrow(gmap) > 0, "No se pudieron mapear genes a ENTREZ."))
    
    ek <- clusterProfiler::enrichKEGG(
      gene = unique(gmap$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 1  # cortamos por FDR luego en el plot
    )
    ek
  })
  
  output$kegg_bar <- renderPlot({
    req(input$enable_viz)
    
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      return(plot_placeholder("Instala 'clusterProfiler' (CRAN/Bioc)"))
    }
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      return(plot_placeholder("Instala 'org.Hs.eg.db' (Bioconductor)"))
    }
    
    df <- viz_df_go(); req(nrow(df) > 0)
    
    gene_col <- NULL
    if ("genename" %in% names(df)) gene_col <- "genename"
    if (is.null(gene_col)) return(plot_placeholder("No existe columna 'genename'"))
    
    genes <- unique(na.omit(as.character(df[[gene_col]])))
    if (!length(genes)) return(plot_placeholder("Sin genes válidos en 'genename'"))
    
    suppressMessages(library(clusterProfiler))
    suppressMessages(library(org.Hs.eg.db))
    conv <- try(clusterProfiler::bitr(genes,
                                      fromType = "SYMBOL",
                                      toType   = "ENTREZID",
                                      OrgDb    = org.Hs.eg.db),
                silent = TRUE)
    if (inherits(conv, "try-error") || is.null(conv) || nrow(conv) == 0) {
      return(plot_placeholder("No se pudieron mapear SYMBOL→ENTREZID"))
    }
    entrez <- unique(na.omit(as.character(conv$ENTREZID)))
    if (!length(entrez)) return(plot_placeholder("Mapeo a ENTREZID vacío"))
    
    qcut <- get_qcut(input$go_kegg_qcut, 0.05)
    
    kk <- try(clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa",
      pvalueCutoff  = 1, qvalueCutoff = 1, pAdjustMethod = "BH"
    ), silent = TRUE)
    if (inherits(kk, "try-error") || is.null(kk) || nrow(as.data.frame(kk)) == 0) {
      return(plot_placeholder("Sin términos KEGG"))
    }
    
    # Aplica FDR por p.adjust
    res <- kk@result
    res <- res[!is.na(res$p.adjust) & res$p.adjust <= qcut, , drop = FALSE]
    if (!nrow(res)) return(plot_placeholder(sprintf("Sin KEGG con FDR ≤ %.3f", qcut)))
    
    kk2 <- kk; kk2@result <- res
    clusterProfiler::barplot(kk2, showCategory = min(20, nrow(res)))
  })
  
  
  
  output$kegg_table <- DT::renderDT({
    ek <- kegg_result(); validate(need(!is.null(ek) && nrow(ek@result) > 0, ""))
    df <- ek@result %>%
      dplyr::mutate(
        p.adjust = p.adjust(pvalue, method = "BH"),
        KEGG = paste0("<a href='https://www.kegg.jp/pathway/", ID, "' target='_blank'>", ID, "</a>")
      ) %>%
      dplyr::arrange(p.adjust)
    DT::datatable(df, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
}

shinyApp(ui, server)

