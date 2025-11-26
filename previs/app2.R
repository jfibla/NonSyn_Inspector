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

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("GWAS Inspector — Etapas 1–3"),
  sidebarLayout(
    sidebarPanel(width = 4,
                 h3("Etapa 1 · Lectura & selección"),
                 fileInput("file", "Carga GWAS (TSV/CSV)", accept = c(".tsv",".txt",".csv")),
                 checkboxInput("header", "Primera fila con nombres", TRUE),
                 radioButtons("sep", "Separador",
                              c("Tab \\t" = "\t", "Coma ," = ",", "Punto y coma ;" = ";"),
                              selected = "\t"),
                 sliderInput("pthr", "Umbral -log10(P)", min=2, max=20, value=7.3, step=0.1),
                 numericInput("flank", "Flanco (+/- bp)", value = 10000, min = 0, step = 1000),
                 actionButton("build_ranges", "➊ Generar selected_intervals"),
                 h4("Preview selected_intervals.range"),
                 verbatimTextOutput("ranges_preview"),
                 tags$hr(),
                 ## --- UI: Etapa 2 · VCF con PLINK 1.9 (simple system2) ---
                 h3("Etapa 2 · VCF con PLINK 1.9"),
                 helpText("Coloca plink y los ref.* (bed/bim/fam) en la carpeta www/"),
                 textInput("plink_path", "Binario PLINK", value = "www/plink19"),
                 textInput("ref_bfile", "Prefijo --bfile (auto si vacío)", value = ""),
                 numericInput("threads19", "Hilos", 4, min = 1),
                 checkboxInput("out_chr_mt", "Etiquetas chr con MT (1..22,X,Y,MT)", TRUE),
                 checkboxInput("use_rosetta", "Usar Rosetta (Mac ARM con PLINK x86_64)", FALSE),
                 
                 div(
                   actionButton("run_plink_simple", "➋ Extraer VCF (PLINK)"),
                   downloadButton("dl_vcf", "Descargar selected_intervals.vcf.gz")
                 ),
                 h4("Log PLINK"),
                 verbatimTextOutput("plink_log"),
                 
                 tags$hr(),
                 # En tu sidebar (reemplaza las entradas de dbNSFP si ya estaban):
                 
                 # ---- UI: Etapa 3 (simple) ----
                 h3("Etapa 3 · Anotar con dbNSFP"),
                 textInput("dbnsfp_dir",  "Carpeta dbNSFP (contiene *.gz y el buscador)",
                           value = "/Volumes/LaCie_Drive/dbNSFP5.0a"),
                 # Acepta ruta al .class o directamente el nombre de la clase (recomendado: 'search_dbNSFP50a')
                 textInput("dbnsfp_tool", "Buscador (clase Java)", value = "search_dbNSFP50a"),
                 numericInput("java_gb",   "Memoria Java (GB)", value = 6, min = 2, step = 1),
                 actionButton("run_dbnsfp", "➌ Ejecutar dbNSFP"),
                 verbatimTextOutput("dbnsfp_log"),
                 downloadButton("dl_ann",  "Descargar dbNSFP.out")
                 
                 
    ),
    mainPanel(
      h3("Manhattan (Etapa 1)"),
      plotlyOutput("manhattan", height = 520),
      tags$hr(),
      h4("SNPs ≥ umbral (selecciona filas)"),
      DTOutput("hits_tbl"),
      tags$hr(),
      h4("Logs"),
      strong("PLINK 1.9"), verbatimTextOutput("plink_log"),
      strong("dbNSFP"), verbatimTextOutput("dbnsfp_log")
    )
  )
)

# ---------- SERVER ----------
server <- function(input, output, session){
  
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
    chr_lengths <- df %>% group_by(CHR) %>%
      summarise(chr_len = max(BP, na.rm=TRUE), .groups="drop") %>%
      arrange(CHR) %>% mutate(tot = cumsum(chr_len) - chr_len)
    
    dfp <- df %>% inner_join(chr_lengths, by = "CHR") %>%
      arrange(CHR, BP) %>% mutate(pos = BP + tot, colgrp = CHR %% 2)
    
    axis_df <- dfp %>% group_by(CHR) %>% summarise(center = median(pos), .groups="drop")
    
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
      filter(logp >= input$pthr) %>%
      arrange(desc(logp)) %>%
      select(CHR, BP, snp, p = Pval, logp, any_of("a1"))
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
    output$ranges_preview <- renderText(paste(readLines(rng_path, warn = FALSE), collapse = "\n"))
    showNotification(paste0("Archivo de intervalos (4 columnas) generado: ", basename(rng_path)), type = "message")
    cat("\n=== selected_intervals.range ===\n",
        "Ruta: ", rng_path, "\n",
        paste(utils::head(readLines(rng_path, warn = FALSE), 20), collapse = "\n"),
        if (length(readLines(rng_path)) > 20) "\n..." else "",
        "\n=== fin ===\n", sep = "")
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
  
  # ======================================
  # ---------- Etapa 3: dbNSFP ----------
  # ======================================
  ## ======================
  ##   Etapa 3 · dbNSFP (simple; sin pruebas chr)
  ## ======================
  
  ann_path <- reactiveVal(NULL)
  
  observeEvent(input$run_dbnsfp, {
    req(vcf_out_path())
    vcf_in <- normalizePath(vcf_out_path(), mustWork = TRUE)
    
    db_dir  <- input$dbnsfp_dir
    tool_in <- input$dbnsfp_tool   # puede ser "search_dbNSFP50a" o "/.../search_dbNSFP50a.class"
    
    validate(need(dir.exists(db_dir), paste0("No existe la carpeta: ", db_dir)))
    
    # --- existencia mínima de la base ---
    chr_files <- list.files(db_dir, pattern="^dbNSFP5\\.0a_variant\\.chr(\\d+|[XYM])\\.gz$", full.names = TRUE)
    validate(need(length(chr_files) >= 24 &&
                    file.exists(file.path(db_dir, "dbNSFP5.0_gene.gz")),
                  "Faltan ficheros en la base (chr*.gz / dbNSFP5.0_gene.gz)."))
    
    # --- helpers ---
    `%||%` <- function(a,b) if (!is.null(a) && length(a)>0 && !is.na(a)) a else b
    java_mem <- as.integer(input$java_gb %||% 6)
    
    # Descomprimir si viene .vcf.gz
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
    
    # Forzar cabecera VCFv4.0 y añadir reference=hg38 si falta
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
    
    # Construir comando EXACTO que te funciona: java -Xmx{g} search_dbNSFP50a -i ... -o ...
    # Si nos pasaron una ruta al .class, extraemos el nombre de la clase
    main_class <- basename(tool_in)
    main_class <- sub("\\.class$", "", main_class, ignore.case = TRUE)
    build_args <- function(vcf_path, out_path) {
      list(cmd = "java",
           args = c(paste0("-Xmx", java_mem, "g"),
                    main_class,
                    "-i", vcf_path,
                    "-o", out_path))
    }
    
    # Preparación del VCF v4.0
    vcf_plain <- ensure_plain_vcf(vcf_in, workdir)
    vcf_v40   <- file.path(workdir, "selected_intervals.v40.vcf")
    vcf_v40   <- force_vcf_version(vcf_plain, vcf_v40, "VCFv4.0")
    
    # Conteo rápido de variantes
    n_var <- sum(!grepl("^#", readLines(vcf_v40, warn = FALSE)))
    validate(need(n_var > 0, "El VCF (v4.0) no contiene variantes que anotar."))
    
    # Salida
    out_file <- file.path(workdir, "dbnsfp.out")
    
    # Ejecutar dentro de db_dir (así Java encuentra la clase y los chr*.gz)
    # Ejecutar dentro de db_dir (cambiando temporalmente el WD)
    cmd <- build_args(normalizePath(vcf_v40, mustWork = TRUE),
                      normalizePath(out_file, mustWork = FALSE))
    
    withProgress(message = "Ejecutando dbNSFP", value = 0.5, {
      res <- tryCatch({
        old <- getwd()
        on.exit(setwd(old), add = TRUE)
        setwd(db_dir)  # <- MUY IMPORTANTE: ejecutamos en la carpeta de la base
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
    
  })
  
  # Descarga
  output$dl_ann <- downloadHandler(
    filename = function() "dbnsfp_annotation.out",
    content  = function(file) file.copy(ann_path(), file, overwrite = TRUE)
  )
  

  
}

shinyApp(ui, server)

