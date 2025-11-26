# app.R
library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(plotly)

ui <- fluidPage(
  titlePanel("Manhattan: -log10(P)"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 h4("Archivo"),
                 fileInput("file", "Carga el fichero (TSV/CSV)", accept = c(".tsv", ".txt", ".csv")),
                 checkboxInput("header", "Primera fila con nombres", TRUE),
                 radioButtons("sep", "Separador",
                              c("Tab \\t" = "\t", "Coma ," = ",", "Punto y coma ;" = ";"),
                              selected = "\t"),
                 tags$hr(),
                 h4("Gráfico"),
                 sliderInput("pthr", "Umbral p (línea horizontal)", min = 1e-12, max = 1e-2,
                             value = 5e-8, step = 1e-12, ticks = FALSE),
                 numericInput("maxn", "Máx. puntos a mostrar (submuestreo)", value = 50000, min = 1000, step = 1000),

                 tags$hr(),
                 h4("Referencia PLINK 1.9"),
                 textInput("ref_bfile", "Prefijo --bfile (sin extensión)", "/data/TopMed_Imputed_GENPIR"),
                 checkboxInput("chr_mt", "Etiquetas chr con MT (1-22,X,Y,MT)", TRUE),
                 numericInput("threads19", "Hilos (PLINK 1.9)", 4, min = 1),
                 actionButton("make_vcf19", "Generar VCF candidato (PLINK 1.9)"),
                 verbatimTextOutput("vcf19_log"),
                 tags$hr(),
                 h4("dbNSFP"),
                 actionButton("run_dbnsfp", "Anotar con dbNSFP"),
                 br(), br(),
                 verbatimTextOutput("dbnsfp_log"),
                 downloadButton("dl_dbnsfp", "Descargar salida")
    ),
    mainPanel(
      plotlyOutput("manhattan", height = 520)
    )
  )
)


server <- function(input, output, session) {
  
  # al inicio del server (o global.R)
  plink_bin <- reactiveVal(NULL)
  
  ensure_plink <- function() {
    # añade www al PATH y detecta el binario
    ww <- normalizePath("www", mustWork = FALSE)
    if (dir.exists(ww)) Sys.setenv(PATH = paste(ww, Sys.getenv("PATH"), sep = ":"))
    cand <- c(
      file.path(ww, "plink"),          # mac/linux
      file.path(ww, "plink1.9"),       # por si lo nombraste así
      file.path(ww, "plink.exe"),      # windows
      Sys.which("plink")
    )
    cand <- unique(path.expand(cand))
    cand <- cand[file.exists(cand)]
    if (length(cand)) {
      # intenta activar permisos de ejecución (unix)
      try(Sys.chmod(cand[1], mode = "0755"), silent = TRUE)
      # comprueba permiso de ejecución
      if (.Platform$OS.type == "windows" || file.access(cand[1], 1) == 0) {
        plink_bin(cand[1])
      } else {
        plink_bin(NULL)
      }
    } else {
      plink_bin(NULL)
    }
  }
  
  # llama una vez al arrancar
  observe({
    ensure_plink()
  })
  
  # botón opcional para forzar detección
  observeEvent(input$detect_plink, { ensure_plink()
    output$plink_check <- renderText({
      if (is.null(plink_bin())) "PLINK no encontrado en www/ ni en PATH."
      else paste("Usando:", plink_bin())
    })
  })
  
  
  # Lectura robusta (soporta comas decimales)
  dat_raw <- reactive({
    req(input$file)
    read_delim(
      file   = input$file$datapath,
      delim  = input$sep,
      col_names = input$header,
      locale = locale(decimal_mark = ",", grouping_mark = "."),
      show_col_types = FALSE,
      na = c("", "NA")
    )
  })
  
  # Normalización de columnas esperadas y cálculo de -log10(P)
  dat <- reactive({
    df <- dat_raw()
    names(df) <- tolower(names(df))
    req(all(c("chr","bp","p") %in% names(df)))
    if (!"snp" %in% names(df)) df$snp <- NA_character_
    
    # --- CHR a num (X=23, Y=24, MT=25)
    chr_map <- function(x){
      x <- toupper(as.character(x))
      x <- sub("^CHR", "", x)
      x[x == "X"]  <- "23"
      x[x == "Y"]  <- "24"
      x[x %in% c("MT","M","MTDNA")] <- "25"
      suppressWarnings(as.integer(x))
    }
    
    # --- BP: números seguros
    BP <- if (is.numeric(df$bp)) df$bp else {
      readr::parse_number(as.character(df$bp), locale = readr::locale())
    }
    
    # --- P: admite punto o coma, incluida notación científica con coma
    if (is.numeric(df$p)) {
      Pval <- df$p
    } else {
      p_chr <- trimws(as.character(df$p))
      p_chr <- gsub("\\s+", "", p_chr)
      # si NO hay punto y sí hay coma -> cambia coma por punto (1,30E-06 -> 1.30E-06)
      needs_swap <- !grepl("\\.", p_chr) & grepl(",", p_chr)
      p_chr[needs_swap] <- gsub(",", ".", p_chr[needs_swap], fixed = TRUE)
      # quitar separadores extraños por si acaso
      Pval <- suppressWarnings(as.numeric(p_chr))
    }
    
    df <- df %>%
      mutate(
        CHR  = chr_map(chr),
        BP   = as.numeric(BP),
        Pval = as.numeric(Pval),
        logp = -log10(Pval)
      ) %>%
      filter(!is.na(CHR), !is.na(BP), !is.na(Pval), Pval > 0)
    
    # posiciones acumuladas por cromosoma
    chr_lengths <- df %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP, na.rm = TRUE), .groups="drop") %>%
      arrange(CHR) %>%
      mutate(tot = cumsum(chr_len) - chr_len)
    
    df <- df %>%
      inner_join(chr_lengths, by = "CHR") %>%
      arrange(CHR, BP) %>%
      mutate(pos = BP + tot)
    
    if (nrow(df) > input$maxn) df <- dplyr::slice_sample(df, n = input$maxn)
    df
  })
  
  
  output$manhattan <- renderPlotly({
    df <- dat(); req(nrow(df) > 0)
    
    # etiquetas de eje X en el centro de cada cromosoma
    axis_df <- df %>% group_by(CHR) %>% summarise(center = median(pos), .groups="drop")
    
    # colores alternos por cromosoma
    df <- df %>% mutate(colgrp = CHR %% 2)
    
    p_thr <- input$pthr
    thr_y <- -log10(p_thr)
    
    g <- ggplot(df, aes(x = pos, y = logp,
                        text = paste0(
                          "SNP: ", ifelse(is.na(snp), "(sin rsID)", snp),
                          "<br>CHR: ", CHR,
                          "<br>BP: ", BP,
                          "<br>P: ", signif(Pval, 3)
                        ))) +
      geom_point(aes(color = factor(colgrp)), size = 0.6) +
      geom_hline(yintercept = thr_y, linetype = "dashed") +
      scale_color_manual(values = c("#377eb8", "#ff7f00")) +
      scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
      labs(x = "Cromosoma", y = expression(-log[10](P))) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 9))
    
    ggplotly(g, tooltip = "text") %>%
      layout(xaxis = list(title = "Cromosoma"),
             yaxis = list(title = "-log10(P)"))
  })
  
  observeEvent(input$make_vcf19, {
    req(input$file, nzchar(input$ref_bfile))
    
    # 1) Toma tu tabla subida y construye lista CHR:BP (una línea por SNP)
    df <- dat_raw()
    names(df) <- tolower(names(df))
    validate(need(("chr" %in% names(df) & "bp" %in% names(df)) | "snp" %in% names(df),
                  "Necesito columnas CHR y BP o una columna SNP tipo 'chr:pos(:ref:alt)'"))
    
    # Mapear CHR a enteros PLINK 1.9 (X=23, Y=24, MT=26)
    chr_map <- function(x){
      x <- toupper(as.character(x))
      x <- sub("^CHR", "", x)
      x[x=="X"] <- "23"; x[x=="Y"] <- "24"
      x[x %in% c("MT","M","MTDNA")] <- "26"   # PLINK1.9 usa 26 para MT
      suppressWarnings(as.integer(x))
    }
    BP <- if (is.numeric(df$bp)) df$bp else readr::parse_number(as.character(df$bp))
    CHR <- if ("chr" %in% names(df)) chr_map(df$chr) else NA_integer_
    
    # Si no hay CHR/BP pero sí SNP con “chr:pos”, extrae de ahí
    if (all(is.na(CHR)) && "snp" %in% names(df)) {
      s <- toupper(as.character(df$snp))
      m <- regexec("^([0-9XYMT]+):([0-9]+)(?::[ACGTN]+:[ACGTN]+)?$", s, perl=TRUE)
      mm <- regmatches(s, m)
      chr_from <- sapply(mm, function(x) if (length(x)>=3) x[2] else NA)
      bp_from  <- sapply(mm, function(x) if (length(x)>=3) x[3] else NA)
      CHR <- chr_map(chr_from)
      BP  <- suppressWarnings(as.numeric(bp_from))
    }
    
    keep <- !is.na(CHR) & !is.na(BP)
    validate(need(any(keep), "No se pudieron derivar posiciones CHR/BP válidas."))
    
    # 2) Rango por SNP (start=end) para --extract range
    work <- file.path(tempdir(), "gwas_inspector_plink19"); dir.create(work, showWarnings = FALSE)
    ranges <- data.frame(chr=CHR[keep], start=as.integer(BP[keep]), end=as.integer(BP[keep]))
    # PLINK 1.9 acepta 3 o 4 columnas: CHR START END [LABEL]
    ranges_path <- file.path(work, "targets.range")
    write.table(ranges, file=ranges_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    
    # 3) Ejecuta PLINK 1.9: extrae por rango y exporta a VCF bgzip
    out_prefix <- file.path(work, "TopMed_Imputed_GENPIR_candi")
    plog <- c()
    
    # A) Extraer a BED temporal (más rápido/seguro antes de recode)
    cmdA <- c("--bfile", normalizePath(input$ref_bfile, mustWork = FALSE),
              "--threads", input$threads19,
              "--chr-set", "26",
              "--extract", "range", ranges_path,
              "--make-bed",
              "--out", file.path(work, "candi_tmp"))
    req(!is.null(plink_bin()), file.exists(plink_bin()))
    pxA <- processx::run(plink_bin(), cmdA, error_on_status = FALSE)
    plog <- c(plog, "CMD A: plink ", paste(cmdA, collapse=" "),
              "\nSTATUS:", pxA$status, "\nSTDERR:\n", pxA$stderr, "\n")
    
    validate(need(pxA$status==0 && file.exists(file.path(work,"candi_tmp.bed")),
                  "PLINK 1.9 (A) falló al extraer por rango."))
    
    # B) Exportar a VCF bgzip (con etiquetas chr estándar)
    cmdB <- c("--bfile", file.path(work, "candi_tmp"),
              "--threads", input$threads19,
              if (isTRUE(input$chr_mt)) c("--output-chr","MT") else NULL,
              "--recode", "vcf", "bgz",
              "--out", out_prefix)
    pxB <- processx::run(plink_bin(), cmdB, error_on_status = FALSE)
    plog <- c(plog, "CMD B: plink ", paste(cmdB, collapse=" "),
              "\nSTATUS:", pxB$status, "\nSTDERR:\n", pxB$stderr, "\n")
    
    output$vcf19_log <- renderText(paste(plog, collapse="\n"))
    
    validate(need(pxB$status==0 && file.exists(paste0(out_prefix, ".vcf.gz")),
                  "PLINK 1.9 no generó el VCF .vcf.gz"))
    
    # Deja el VCF listo para dbNSFP
    vcf_path <<- paste0(out_prefix, ".vcf.gz")
    showNotification("VCF candidato generado (PLINK 1.9).", type="message")
  })
  
  observeEvent(input$run_dbnsfp, {
    req(exists("vcf_path"), file.exists(vcf_path))
    out <- sub("\\.vcf\\.gz$", ".out", vcf_path)
    tool <- Sys.getenv("DBNSFP_TOOL", unset = "search_dbNSFP50a") # o ruta .jar
    
    validate(need(nzchar(Sys.which("java")), "No se encontró Java en PATH."))
    args <- c("-Xmx5g", tool, "-i", vcf_path, "-o", out) # si es .jar: c("-Xmx5g","-jar", tool, ...)
    px <- processx::run("java", args=args, error_on_status = FALSE)
    
    output$dbnsfp_log <- renderText(paste0(
      "CMD: java ", paste(args, collapse=" "),
      "\nExit status: ", px$status, "\n\nSTDOUT:\n", px$stdout, "\n\nSTDERR:\n", px$stderr
    ))
    validate(need(px$status==0 && file.exists(out), "dbNSFP falló; revisa el log."))
    
    output$dl_dbnsfp <- downloadHandler(
      filename = basename(out),
      content = function(file) file.copy(out, file, overwrite = TRUE)
    )
  })
  
}

shinyApp(ui, server)
