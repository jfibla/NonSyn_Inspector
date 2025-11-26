# setwd("~/Library/CloudStorage/Dropbox/Signature_Selection_by_Shiny/SIG_SELEC_REHH_CANDI_REGIONS/SHINY_SELECTION_REHH/Shiny_selection_REHH_dbNSFP50")
# 
library(shiny)
library(plotly) # For interactive plots
library(readr)
library(data.table)
library(ggplot2)
library(shinyjs)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(rehh)
library(conflicted)
library(ggplotify)
library(DT)
library(grid)
library(ComplexHeatmap)
library(shinycssloaders)
library(Hmisc)
library(ggpubr)
library(stringr)
library(tibble)
library(AnnotationDbi)
library(GO.db)
library(wordcloud)
library(tm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(patchwork)

conflicts_prefer(
  plotly::subplot,
  dplyr::filter,
  plotly::layout,
  AnnotationDbi::select
)

# Load your datasets
load(file = "www/dbNSFP50_candi_unrel_1KG_hg38.rda") # input file plot scores candidates hg38 regions
load(file = "www/dbNSFP50_IMP.rda") # input file plot scores candidates Imputed GENPIR regions
load(file = "www/candi_gene_regions.rda") # input file plot genomic structure
load(file = "www/candi_freqs_95CI.all_long.rda") 
load(file = "www/sorted_df.rda") # input file for haplotype plots
load(file = "www/f.3.rs4000279.rda") # input file furcation
load(file = "www/f.5.rs55740022.rda") # input file furcation
load(file = "www/f.9.rs75427350.rda") # input file furcation
load(file = "www/f.11.rs72923901.rda") # input file furcation
load(file = "www/f.17.rs62050498.rda") # input file furcation
load(file = "www/f.20.rs2144899.rda") # input file furcation

load(file = "www/ihs_CEU.rda")
load(file = "www/ihs_IBS.rda")
load(file = "www/ihs_PIR.rda")

#head(dbNSFP50.candi_scores)
#str(dbNSFP50.candi_scores)


# Define the UI
ui <- fluidPage(
  titlePanel("Selection signatures at Catalan Pyrenean region based on REHH"),
  useShinyjs(),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Selector principal
      selectInput("study", "Select study:", 
                  choices = c("Candidate region compared" = "ihs",
                              "Haplotype furcation" = "rehh",
                              "Non Synonymous candidates" = "ns",
                              "Annotations GO and KEGG" = "gokegg")),
      
      # Archivos de entrada para todo excepto 'ihs'
      conditionalPanel(
        condition = "input.study != 'ihs'",
        selectInput("inputfile", "Select input file:",
                    choices = c("1000KG candidate regions" = "kg",
                                "GENPIR imputed" = "imputed",
                                "Common candidate regions" = "common"))
      ),
      
      # Chromosome selector (condicional según estudio)
      conditionalPanel(
        condition = "input.study == 'gokegg'",
        selectInput("chromosome", "Select chromosome:", 
                    choices = c("chr3" = 3, "chr9" = 9, "chr11" = 11,
                                "chr17" = 17, "chr20" = 20, "all set" = 0))
      ),
      conditionalPanel(
        condition = "input.study != 'gokegg'",
        selectInput("chromosome", "Select chromosome:",
                    choices = c("chr3" = 3, "chr9" = 9, "chr11" = 11,
                                "chr17" = 17, "chr20" = 20))
      ),
      
      # Parámetros para KEGG
      conditionalPanel(
        condition = "input.study == 'gokegg'",
        numericInput("kegg_pval", "KEGG p-value cutoff", 1, min = 0, max = 1, step = 0.01),
        verbatimTextOutput("kegg_text")
      ),
      
      # Parámetros específicos de 'rehh'
      conditionalPanel(
        condition = "input.study == 'rehh'",
        selectInput("limits", "Select x limits:",
                    choices = c("Extended plot" = "extended", "Reduced plot" = "reduced")),
        
        selectInput("scoreorrank", "Select score/rank type:", choices = c("Score", "Score rank")),
        
        conditionalPanel(
          condition = "input.scoreorrank == 'Score'",
          selectInput("scores", "Select score:", choices = c(
            "SIFT_score", "SIFT4G_score", "Polyphen2_HDIV_score",
            "Polyphen2_HVAR_score", "MutationTaster_score", "MutationAssessor_score",
            "PROVEAN_score", "REVEL_score", "GERP_91_mammals",
            "phyloP17way_primate", "phastCons17way_primate"
          ))
        ),
        
        conditionalPanel(
          condition = "input.scoreorrank == 'Score rank'",
          selectInput("rank", "Select rank score:", choices = c(
            "SIFT_converted_rankscore", "SIFT4G_converted_rankscore", "Polyphen2_HDIV_rankscore",
            "Polyphen2_HVAR_rankscore", "MutationTaster_converted_rankscore",
            "MutationAssessor_rankscore", "PROVEAN_converted_rankscore",
            "GERP++_RS_rankscore", "GERP_91_mammals_rankscore",
            "phyloP470way_mammalian_rankscore", "phyloP17way_primate_rankscore",
            "phastCons17way_primate_rankscore"
          ))
        ),
        
        selectInput("decorate", "Select predictor:",
                    choices = c("Score predictors", "Gene demage predictors")),
        
        conditionalPanel(
          condition = "input.decorate == 'Gene demage'",
          selectInput("genedemage", "Gene damage predictors:", choices = c(
            "all_disease" = "Gene_damage_prediction_all_disease_causing_genes",
            "all_Mendelian" = "Gene_damage_prediction_all_Mendelian_disease_causing_genes",
            "Mendelian_AD" = "Gene_damage_prediction_Mendelian_AD_disease_causing_genes",
            "Mendelian_AR" = "Gene_damage_prediction_Mendelian_AR_disease_causing_genes",
            "all_PID" = "Gene_damage_prediction_all_PID_disease_causing_genes",
            "PID_AD" = "Gene_damage_prediction_PID_AD_disease_causing_genes",
            "PID_AR" = "Gene_damage_prediction_PID_AR_disease_causing_genes",
            "all_cancer" = "Gene_damage_prediction_all_cancer_disease_causing_genes",
            "cancer_recessive" = "Gene_damage_prediction_cancer_recessive_disease_causing_genes",
            "cancer_dominant" = "Gene_damage_prediction_cancer_dominant_disease_causing_genes"
          ))
        ),
        
        conditionalPanel(
          condition = "input.decorate == 'Score predictors'",
          selectInput("predictors", "Score predictors:", choices = c(
            "SIFT_pred", "SIFT4G_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred",
            "MutationTaster_pred", "MutationAssessor_pred", "PROVEAN_pred",
            "MetaSVM_pred", "MetaLR_pred", "MetaRNN_pred", "M-CAP_pred", "PrimateAI_pred",
            "DEOGEN2_pred", "BayesDel_addAF_pred", "BayesDel_noAF_pred", "ClinPred_pred",
            "LIST-S2_pred", "ESM1b_pred", "AlphaMissense_pred", "Aloft_pred",
            "fathmm-XF_coding_pred", "HIPred"
          ))
        ),
        
        selectInput("genome", "Select genome entity:", choices = c(
          "Coding region" = "CDS", "Transcript" = "transcript", "Exons" = "exon",
          "5' UTR" = "5UTR", "3' UTR" = "3UTR", "Start codon" = "start_codon",
          "Stop codon" = "stop_codon"
        ))
      )
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "input.study == 'ihs'",
        withSpinner(plotOutput("manhat_plot", height = 1200), type = 1, color = "#007BFF"),
        #plotOutput("manhat_plot", height = 800)
      ),
      
      conditionalPanel(
        condition = "input.study == 'ns'",
        div(style = "height:800px; overflow-y: scroll; border: 1px solid #ccc;",
            uiOutput("dynamic_heatmap_plot")),
        DTOutput("snp_table")
      ),
      
      conditionalPanel(
        condition = "input.study == 'rehh'",
        plotOutput("focal_plot", height = 400, width = 900),
        plotlyOutput("combined_plot")
      ),
      
      conditionalPanel(
        condition = "input.study == 'gokegg'",
        withSpinner(plotOutput("go_plot", height = 400), type = 1, color = "#007BFF"),
        plotOutput("kegg_plot"),
        tags$h4(strong("GO and KEGG terms")),  # Título en negrita
        DTOutput("kegg_table"),
        tags$h4(strong("KEGG enrichment results")),  # Título en negrita
        DTOutput("kegg_results_table"),
        plotOutput("kegg_barplot"),
        plotOutput("kegg_wordcloud")
      )
    )
  )
)



# Define the server logic
server <- function(input, output, session) {

##################################################################  
  output$focal_plot <- renderPlot({
    # assign global boundaries
    if (input$limits == "extended") {
      filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
      filtered.sorted_df <- filtered.sorted_df %>%
        mutate(
          haplen.MIN = as.numeric(unlist(haplen.MIN)),
          haplen.MAX = as.numeric(unlist(haplen.MAX)),
          row_number = as.numeric(unlist(row_number))
        ) 
      global_x_min <- min(filtered.sorted_df$haplen.MIN, na.rm = TRUE)
      global_x_max <- max(filtered.sorted_df$haplen.MAX, na.rm = TRUE)
    } else {
    if (input$chromosome == 3) {
      global_x_min <- 93400000
      global_x_max <- 95200000
    } else if (input$chromosome == 5) {
      global_x_min <- 23500000
      global_x_max <- 25400000
    } else if (input$chromosome == 9) {
      global_x_min <- 67700000
      global_x_max <- 69400000
    } else if (input$chromosome == 11) {
      global_x_min <- 55800000
      global_x_max <- 57700000
    } else if (input$chromosome == 17) {
      global_x_min <- 20600000
      global_x_max <- 22400000
    } else if (input$chromosome == 20) {
      global_x_min <- 15300000
      global_x_max <- 17200000
    }
    }
# focal plot by chromosome
    if (input$chromosome == 3) {
      focal_position <- f.3.rs4000279@position
      buffer_min <- focal_position - global_x_min  # Extend the xlim by 500kb on each side
      buffer_max <-  global_x_max -focal_position  # Extend the xlim by 500kb on each side
      extended_xlim <- c(focal_position - buffer_min, focal_position + buffer_max)
      plot(f.3.rs4000279, xlim = extended_xlim)  
      
    } else if (input$chromosome == 5) { 
      focal_position <- f.5.rs55740022@position
      buffer_min <- focal_position - global_x_min  # Extend the xlim by 500kb on each side
      buffer_max <-  global_x_max -focal_position  # Extend the xlim by 500kb on each side
      extended_xlim <- c(focal_position - buffer_min, focal_position + buffer_max)
      plot(f.5.rs55740022, xlim = extended_xlim) 
      
    } else if (input$chromosome == 9) {
      focal_position <- f.9.rs75427350@position
      buffer_min <- focal_position - global_x_min  # Extend the xlim by 500kb on each side
      buffer_max <-  global_x_max -focal_position  # Extend the xlim by 500kb on each side
      extended_xlim <- c(focal_position - buffer_min, focal_position + buffer_max)
      plot(f.9.rs75427350, xlim = extended_xlim)
      
    } else if (input$chromosome == 11) {
      focal_position <- f.11.rs72923901@position
      buffer_min <- focal_position - global_x_min  # Extend the xlim by 500kb on each side
      buffer_max <-  global_x_max -focal_position  # Extend the xlim by 500kb on each side
      extended_xlim <- c(focal_position - buffer_min, focal_position + buffer_max)
      plot(f.11.rs72923901, xlim = extended_xlim) 
      
    } else if (input$chromosome == 17) {
      focal_position <- f.17.rs62050498@position
      buffer_min <- focal_position - global_x_min  # Extend the xlim by 500kb on each side
      buffer_max <-  global_x_max -focal_position  # Extend the xlim by 500kb on each side
      extended_xlim <- c(focal_position - buffer_min, focal_position + buffer_max)
      plot(f.17.rs62050498, xlim = extended_xlim) 
        
    } else if (input$chromosome == 20) {
      focal_position <- f.20.rs2144899@position
      buffer_min <- focal_position - global_x_min  # Extend the xlim by 500kb on each side
      buffer_max <-  global_x_max -focal_position  # Extend the xlim by 500kb on each side
      extended_xlim <- c(focal_position - buffer_min, focal_position + buffer_max)
      plot(f.20.rs2144899, xlim = extended_xlim) 
    }
  
  })
  ##################################################################  
  
  output$combined_plot <- renderPlotly({
 # Filter and prepare input file for haplotype plot
    if (input$limits == "extended") {
      filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
      filtered.sorted_df <- filtered.sorted_df %>%
        mutate(
          haplen.MIN = as.numeric(unlist(haplen.MIN)),
          haplen.MAX = as.numeric(unlist(haplen.MAX)),
          row_number = as.numeric(unlist(row_number))
        ) 
      global_x_min <- min(filtered.sorted_df$haplen.MIN, na.rm = TRUE)
      global_x_max <- max(filtered.sorted_df$haplen.MAX, na.rm = TRUE)
    } else if (input$limits == "reduced") {
      if (input$chromosome == 3) {
        global_x_min <- 93400000
        global_x_max <- 95200000
      } else if (input$chromosome == 5) {
        global_x_min <- 23500000
        global_x_max <- 25400000
      } else if (input$chromosome == 9) {
        global_x_min <- 67700000
        global_x_max <- 69400000
      } else if (input$chromosome == 11) {
        global_x_min <- 55800000
        global_x_max <- 57700000
      } else if (input$chromosome == 17) {
        global_x_min <- 20600000
        global_x_max <- 22400000
      } else if (input$chromosome == 20) {
        global_x_min <- 15300000
        global_x_max <- 17200000
      }
    }

    # Filter and prepare data for plots of score candidates
    # generate common input score file
    filter <- as.data.frame(dbNSFP50_IMP$SNP)
    colnames(filter) <- "SNP"
    common <- merge(dbNSFP50_candi_unrel_1KG_hg38, filter, by="SNP")
    
    # switch between input score
    score_input <- switch(
      input$inputfile,
      "kg" = dbNSFP50_candi_unrel_1KG_hg38,
      "imputed" = dbNSFP50_IMP,
      "common" = common
    )
    
    # Filter and prepare data for frequency plot
    filtered_freqs <- candi_freqs_95CI.all_long %>% filter(CHR == as.numeric(input$chromosome))

    # filter freqs by group (=input$inputfile)
    if (input$inputfile == "kg") {
    filtered_freqs <-   filtered_freqs  %>% filter(group %in% c("kg")) 
    filtered_freqs <-   na.omit(filtered_freqs)
    } else if (input$inputfile == "imputed") {
      filtered_freqs <-   filtered_freqs  %>% filter(group %in% c("imputed")) 
      filtered_freqs <-   na.omit(filtered_freqs)
    } else {
      filtered_freqs <-   filtered_freqs  %>% filter(group == c("common"))
    }
    
    
    if (input$chromosome != 0) {
      score_input <- score_input %>% filter(chr == as.numeric(input$chromosome))
    }
    ##################################################################  
  
    if (input$limits == "extended") {
      filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
      filtered.sorted_df <- filtered.sorted_df %>%
        mutate(
          haplen.MIN = as.numeric(unlist(haplen.MIN)),
          haplen.MAX = as.numeric(unlist(haplen.MAX)),
          row_number = as.numeric(unlist(row_number))
        ) 
      
      global_x_min <- min(filtered.sorted_df$haplen.MIN, na.rm = TRUE)
      global_x_max <- max(filtered.sorted_df$haplen.MAX, na.rm = TRUE)
    } else if (input$limits == "reduced") {
      if (input$chromosome == 3) {
        filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
        global_x_min <- 93400000
        global_x_max <- 95200000
      } else if (input$chromosome == 5) {
        filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
        global_x_min <- 23500000
        global_x_max <- 25400000
      }else if (input$chromosome == 9) {
        filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
        global_x_min <- 67700000
        global_x_max <- 69400000
      }else if (input$chromosome == 11) {
        filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
        global_x_min <- 55800000
        global_x_max <- 57700000
      }else if (input$chromosome == 17) {
        filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
        global_x_min <- 20600000
        global_x_max <- 22400000
      }else if (input$chromosome == 20) {
        filtered.sorted_df <- sorted_df %>% filter(chr == input$chromosome)
        global_x_min <- 15300000
        global_x_max <- 17200000
      }
    }
 
   if (input$chromosome == 3) {
     global_min <- 93400000
     global_max <- 95200000
   } else if (input$chromosome == 5) {
     global_min <- 23500000
     global_max <- 25400000
   }else if (input$chromosome == 9) {
     global_min <- 67700000
     global_max <- 69400000
   }else if (input$chromosome == 11) {
     global_min <- 55800000
     global_max <- 57700000
   }else if (input$chromosome == 17) {
     global_min <- 20600000
     global_max <- 22400000
   }else if (input$chromosome == 20) {
     global_min <- 15300000
     global_max <- 17200000
   }
 
    # Haplotype plot
    hap_plot <- ggplot(filtered.sorted_df, aes()) +
      geom_segment(
        aes(
          x = haplen.MIN,
          xend = haplen.MAX,
          y = row_number,
          yend = row_number,    # Explicitly set yend to the same as y
          color = haplen.DESCRIPTION
        ),
        size = 1.2
      ) +
      geom_vline(
        aes(xintercept = global_min),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      geom_vline(
        aes(xintercept = global_max),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      labs(
        title = "Haplotype Lengths by Coordinates",
        x = "Genomic Position",
        y = "Number of haplotypes",
        color = "Haplotype Type"
      ) +
      scale_x_continuous(limits = c(global_x_min, global_x_max)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top") +
      scale_color_manual(values = c("Ancestral" = "blue", "Derived" = "red"))
    
 #   })
################################################################################
    # Frequency Plot
    p_freq <- ggplot(filtered_freqs, aes(x = BP, y = MAF, color = population)) +
      geom_point(aes(
        text = paste(
          "SNP:", SNP,
          "\nPosition :", paste0("chr",input$chromosome,":",BP),
          "\nPopulation:", population,
          "\nMAF:", MAF
        )
      ), size = 3) +
      scale_x_continuous(limits = c(global_x_min, global_x_max)) +
      facet_grid(~CHR, scales = "free_x") + 
      geom_errorbar(aes(ymin = lw, ymax = up), width = 0.2, size = 0.7) +
      geom_vline(
        aes(xintercept = global_min),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      geom_vline(
        aes(xintercept = global_max),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      labs(
        title = "Frequency Plot",
        x = "Position",
        y = "Minor Allele Frequency (MAF)"
      ) 
################################################################################
    # Score Plot
    p_score <- ggplot() +
      geom_point(
        data = score_input,
        aes(
          x = BP,
          y = if (input$scoreorrank == "Score") .data[[input$scores]] else .data[[input$rank]],
          color = if (input$decorate == "Score predictors") .data[[input$predictors]] else .data[[input$genedemage]],
          text = paste(
            "SNP:", SNP,
            "\nPosition :", paste0("chr",input$chromosome,":",BP),
            if (input$scoreorrank == "Score") paste("<br>Score:", .data[[input$scores]]) else paste("<br>Rank:", .data[[input$rank]]),
            if (input$decorate == "Score predictors") paste("<br>Prediction:", .data[[input$predictors]]) else paste("<br>Gene Damage:", .data[[input$genedemage]])
          )
        ),
        size = 3
      ) +
     scale_x_continuous(limits = c(global_x_min, global_x_max)) +
      scale_color_manual(
        values = c("D" = "red", "T" = "orange", "grey" = "grey", "High" = "red", "Medium" = "orange", "Low" = "green"),
        guide = "none"
      ) +
      facet_grid(~chr, scales = "free_x") + 
      geom_vline(
        aes(xintercept = global_min),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      geom_vline(
        aes(xintercept = global_max),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      labs(
        title = "Score Plot",
        x = "Position",
        y = if (input$scoreorrank == "Score") input$scores else input$rank
      ) 

################################################################################
    # Gene annotations Plot
    if (input$chromosome != "all") {
      chromosome <- paste0("chr",as.numeric(input$chromosome))
      filtered_regions <- candi_gene_regions %>% filter(seqnames == chromosome)
      } else {filtered_regions <- candi_gene_regions}
    
    # Separate exons and group by genes transcript exon 5UTR CDS 3UTR start_codon stop_codon
    exons <- filtered_regions %>%
      filter(type == input$genome) %>%
      group_by(gene_id) %>%
      arrange(start)
    
    # Assign a numeric index to each unique gene for plotting
    gene_levels <- unique(exons$gene_id)
    exons <- exons %>%
      mutate(gene_index = as.numeric(factor(gene_id, levels = gene_levels)))
    
    # Create a plot
    p_gene <-  ggplot() +
      # Plot exons as rectangles
      geom_rect(
        data = exons,
        aes(
          xmin = start,
          xmax = end,
          ymin = gene_index - 0.2,
          ymax = gene_index + 0.2,
        #  fill = gene_id,
          text = paste(
            "<br>Gene name:", gene_name,
            "<br>Transcript ID:", transcript_id,
            "\nStart: ", start,
            "\nEnd: ", end
          )
        ),
        color = "black"
      ) +
      # Plot introns as lines connecting exons
      geom_segment(
        data = exons %>% group_by(gene_id) %>% mutate(next_start = lead(start)),
        aes(
          x = end,
          xend = next_start,
          y = gene_index,
          yend = gene_index,
          text = paste(
            "<br>Gene name:", gene_name,
            "<br>Transcript ID:", transcript_id,
            "\nStart: ", start,
            "\nEnd: ", end
          )
        ),
        color = "black",
        linetype = "dotted"
      ) +
      scale_x_continuous(limits = c(global_x_min, global_x_max)) +
      # Customize the plot
      scale_y_continuous(
        breaks = 1:length(gene_levels),  # Match breaks to gene indices
        labels = gene_levels            # Use gene IDs as labels
      ) +
      geom_vline(
        aes(xintercept = global_min),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      geom_vline(
        aes(xintercept = global_max),
        color = "black", linetype = "dashed", size = 0.3
      ) +
      labs(
        title = paste("Score distribution at candidate region of chromosome ", input$chromosome),
        x = "Genomic Position",
        y = "Gene"
      ) +
      facet_grid(~seqnames, scales = "free_x") + 
   #   theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        legend.position = "none"
      )
    
    ##################################################
 
    p_combined <- subplot(
      ggplotly(hap_plot),
      ggplotly(p_freq, tooltip = "text"),
      ggplotly(p_score, tooltip = "text"),
      ggplotly(p_gene, tooltip = "text"),
      nrows = 4,
      heights = c(0.3,0.2,0.2,0.3),
      shareX = FALSE 
    ) %>%
      layout(
        height = 1200,  # Adjust overall height
        width = 900,   # Adjust overall width
        margin = list(t = 50, b = 50, l = 50, r = 50)  # Set margins
      )
    
  })
  ######################## Headmap plot ########################################
  score_input <- reactive({
  if (input$inputfile == "kg") {
    dbNSFP50_candi_unrel_1KG_hg38 %>% filter(chr == input$chromosome) %>%
      mutate(chr = replace_na(chr, 23))  # Reemplazar NA por "X"
  } else if (input$inputfile == "imputed") {
    dbNSFP50_IMP %>% filter(chr == input$chromosome) %>%
      mutate(chr = replace_na(chr, 23))  # Reemplazar NA por "X"
  } else {
    filter <- as.data.frame(dbNSFP50_candi_unrel_1KG_hg38$SNP)
    colnames(filter) <- "SNP"
    combi_data <- merge(dbNSFP50_IMP, filter, by = "SNP")
    combined_data <- combi_data %>% filter(chr == input$chromosome) %>%
      mutate(chr = replace_na(chr, 23))  # Reemplazar NA por "X"
  }
  })
  

  ################################################################## 
  
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(reshape2)
  
  output$heatmap_plot <- renderPlot({
    if (is.null(score_input()) || nrow(score_input()) == 0) {
      message("⚠️ score_input es NULL o vacío")
      return(NULL)
    }
    
    # Filtrar las columnas de interés
    sel_score <- score_input() %>% dplyr::select(
      SNP, chr, BP, ref, alt, genename,
      SIFT_score, SIFT_pred, SIFT_converted_rankscore,
      SIFT4G_score, SIFT4G_pred, SIFT4G_converted_rankscore,
      Polyphen2_HDIV_score, Polyphen2_HDIV_rankscore, Polyphen2_HDIV_pred,
      Polyphen2_HVAR_score, Polyphen2_HVAR_rankscore, Polyphen2_HVAR_pred,
      MutationTaster_score, MutationTaster_pred, 
      MutationAssessor_score, MutationAssessor_rankscore, MutationAssessor_pred,
      PROVEAN_score, PROVEAN_converted_rankscore, PROVEAN_pred, 
      REVEL_score, REVEL_rankscore,'GERP++_RS_rankscore',
      GERP_91_mammals, GERP_91_mammals_rankscore,
      phyloP17way_primate, phyloP17way_primate_rankscore,
      phastCons17way_primate, phastCons17way_primate_rankscore
    )
    
    df <- if (input$scoreorrank == "Score") {
      sel_score %>% dplyr::select(SNP, chr, BP,ends_with("_score"),GERP_91_mammals,phyloP17way_primate,phastCons17way_primate)
    } else {
      sel_score %>% dplyr::select(SNP, chr, BP,ends_with("_rankscore"))
    }
    
    # Reemplazar "_" por "\n" en los nombres de las columnas
    colnames(df) <- gsub("_", "\n", colnames(df))
    
    df <- df %>%
      mutate(across(-c(SNP, chr,BP), as.numeric)) %>%
      arrange(chr,BP) # Ordenar por chr
    
    score_matrix <- df %>% dplyr::select(-c(SNP, chr)) %>% as.matrix()
    rownames(score_matrix) <- df$SNP
    
    # Normalizar
    score_matrix_normalized <- apply(score_matrix, 2, function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    })
    
    
    # Crear heatmap
    heatmap_plot <- as.ggplot(function() {
      pheatmap::pheatmap(score_matrix_normalized, 
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               display_numbers = score_matrix, 
               number_format = "%.2f", 
               main = "Heatmap Scores", 
               fontsize_number = 10)
    })
    
    heatmap_plot
  })
  
  output$dynamic_heatmap_plot <- renderUI({
    n_rows <- nrow(score_input())
    
    # Altura base: 20 px por fila, mínimo 400 px
    plot_height <- max(400, n_rows * 20)
    
    plotOutput("heatmap_plot", height = paste0(plot_height, "px"))
  })
#################################################################################


  output$snp_table <- renderDT({
    df <- score_input() %>%
      filter(chr == input$chromosome) %>%
      dplyr::select(SNP, chr, BP, genename,MIM_disease) %>%
      arrange(chr, BP) %>%
      mutate(
        genename = sapply(strsplit(as.character(genename), ";"), function(x) {
          paste(unique(x), collapse = ";")  # Mantener solo elementos únicos
        }),
        SNP = paste0("<a href='https://www.ncbi.nlm.nih.gov/snp/", SNP, 
                     "' target='_blank'>", SNP, "</a>"),
        genename = ifelse(!is.na(genename), 
                          paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", genename, 
                                 "' target='_blank'>", genename, "</a>"), 
                          NA)
      )
    
    datatable(df, escape = FALSE, selection = "none")  # Deshabilita selección ya que los enlaces están embebidos
  })
  
  
  
#################################################################################

  ############################## GO Enrichment ################################
  ################################################################################
  
  ################################################################## 
  
  dbNSFP50_candi_unrel_1KG_hg38 <- dbNSFP50_candi_unrel_1KG_hg38 %>%
    mutate(chr = replace_na(chr, 23))  # Reemplazar NA por "X"
  
  dbNSFP50_IMP <- dbNSFP50_IMP %>%
    mutate(chr = replace_na(chr, 23))  # Reemplazar NA por "X"
  
  combined_data <- reactive({
    score_input1 <- dbNSFP50_candi_unrel_1KG_hg38
    score_input2 <- dbNSFP50_IMP
    filter <- as.data.frame(score_input2$SNP)
    colnames(filter) <- "SNP"
    combined_data <- merge(score_input1, filter, by = "SNP")
  })
  
  go_input <- reactive({
    score_input <- switch(
      input$inputfile,
      "kg" = dbNSFP50_candi_unrel_1KG_hg38,
      "imputed" = dbNSFP50_IMP,
      "common" = combined_data()
    )
    if (input$chromosome != 0) {
      score_input %>%
        filter(chr == input$chromosome)
    } else {
      score_input
    }
  })
  
  
  ###############################################################################   
  output$go_plot <- renderPlot({
    dataset <- go_input() 
    
    # Mapas de columnas por ontología
    ontology_map <- list(
      BP = dataset$GO_biological_process,
      CC = dataset$GO_cellular_component,
      MF = dataset$GO_molecular_function
    )
    
    all_results <- purrr::map_dfr(names(ontology_map), function(ont) {
      input_column <- ontology_map[[ont]]
      input_column <- na.omit(input_column)
      
      go_terms <- AnnotationDbi::select(
        GO.db,
        keys = keys(GO.db),
        columns = c("TERM", "ONTOLOGY"),
        keytype = "GOID"
      ) %>%
        filter(ONTOLOGY == ont)
      
      go_terms$TERM <- tolower(go_terms$TERM)
      all_go_terms <- unique(go_terms$TERM)
      
      input_list <- unique(input_column)
      input_terms <- unlist(str_split(input_list, ";")) %>%
        str_trim() %>%
        tolower()
      input_terms <- input_terms[input_terms != ""]
      
      sample_term_table <- table(input_terms)
      background_term_table <- table(
        factor(all_go_terms, levels = unique(c(names(sample_term_table), all_go_terms)))
      )
      
      if (length(sample_term_table) == 0 || length(background_term_table) == 0) return(NULL)
      
      results <- lapply(names(sample_term_table), function(term) {
        a <- sample_term_table[term]
        b <- sum(sample_term_table) - a
        c <- background_term_table[term]
        d <- sum(background_term_table) - c
        
        if (any(is.na(c(a, b, c, d))) || any(c(a, b, c, d) < 0)) return(NULL)
        
        mat <- matrix(c(a, b, c, d), nrow = 2)
        ft <- fisher.test(mat, alternative = "greater")
        
        data.frame(
          term = term,
          count_sample = a,
          count_bg = c,
          p_value = ft$p.value,
          ontology = ont
        )
      }) %>% bind_rows()
      
      if (nrow(results) == 0) return(NULL)
      
      results %>%
        arrange(p_value) %>%
        mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
        filter(p_adj < 0.05) %>%
        slice_head(n = 15) %>%
        mutate(term = str_trunc(term, 50))
    })
    
    if (nrow(all_results) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "Sin términos GO enriquecidos", size = 6) +
          theme_void()
      )
    }
    
    ggplot(all_results, aes(x = reorder(term, -p_adj), y = -log10(p_adj), fill = ontology)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_wrap(~ ontology, scales = "free") +
      labs(
        title = "Top 15 GO terms by ontology",
        x = "GO term",
        y = "-log10(FDR)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # título centrado, negrita, tamaño 18
      )
  })
  
  
  ############################################################################################################
  
  keeg_input <- reactive({
    print("⚙️ kegg_input se está ejecutando")
    score_input <- switch(
      input$inputfile,
      "kg" = dbNSFP50_candi_unrel_1KG_hg38,
      "imputed" = dbNSFP50_IMP,
      "common" = combined_data()
    )
    if (input$chromosome != 0) {
      score_input %>%
        filter(chr == input$chromosome)
    } else {
      score_input
    }
    
  })
  
  kegg_result_reactive <- reactive({
    # filter rsid tested at LFMM test
    df <- keeg_input() 
    
    genes_unicos <- df$genename %>%
      strsplit(";") %>%
      unlist() %>%
      unique() %>%
      na.omit()
    
    genes_entrez <- bitr(
      genes_unicos,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    enrichKEGG(
      gene = genes_entrez$ENTREZID,
      organism = "hsa",
      pvalueCutoff = input$kegg_pval
    )
  })
  
  output$kegg_plot <- renderPlot({
    
    kegg_result <- kegg_result_reactive()
    
    validate(need(nrow(kegg_result@result) > 0, "No enriched KEGG terms found."))
    
    bar <- barplot(kegg_result, showCategory = 15)
    cnet <- cnetplot(kegg_result, categorySize = "pvalue", foldChange = NULL)
    (bar + cnet) +
      plot_annotation(
        title = "KEGG Enrichment Analysis",
        theme = theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
        )
      )
  })
  
  output$kegg_text <- renderPrint({
    req(kegg_result_reactive())
    kegg_result_reactive()
  })
  ###############################################################################
  
  output$kegg_results_table <- renderDT({
    kegg_result <- kegg_result_reactive()
    
    validate(need(nrow(kegg_result@result) > 0, "No enriched KEGG terms found."))
    
    result_df <- kegg_result@result %>%
      mutate(
        across(where(is.numeric), ~ round(.x, 3)),
        # Enlace KEGG para la columna ID (pathway)
        ID = paste0(
          "<a href='https://www.kegg.jp/pathway/", ID, 
          "' target='_blank'>", ID, "</a>"
        ),
        # Enlaces múltiples para cada geneID, separados por ;
        geneID = sapply(strsplit(geneID, "/"), function(genes) {
          paste0(
            "<a href='https://www.kegg.jp/entry/hsa:", genes,
            "' target='_blank'>", genes, "</a>"
          ) %>% paste(collapse = "; ")
        })
      )
    
    datatable(
      result_df,
      escape = FALSE,  # necesario para mostrar HTML
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  
  ###############################################################################
  
  output$kegg_table <- renderDT({
    # Filtrar rsid testeados en LFMM
    df <- keeg_input() 
    
    # Filtrar filas completas y eliminar NAs
    df <- df %>%
      filter(
        !is.na(genename),
        !is.na(Pathway_KEGG_id),
        !is.na(Pathway_KEGG_full),
        !is.na(GO_biological_process),
        !is.na(GO_cellular_component),
        !is.na(GO_molecular_function)
      ) %>%
      filter(genename != "", Pathway_KEGG_id != "", Pathway_KEGG_full != "",GO_biological_process!= "",GO_cellular_component!= "",GO_molecular_function!= "")
    
    # Expandir simultáneamente genename y pathway usando pmap
    library(purrr)
    library(tidyr)
    
    # Expande simultáneamente genename y pathway usando pmap
    kegg_pairs <- df %>%
      mutate(row = pmap(list(genename, GO_biological_process, GO_cellular_component, GO_molecular_function, Pathway_KEGG_id, Pathway_KEGG_full), function(g, bp, cc, mf, id, full) {
        genes <- strsplit(g, ";")[[1]]
        bps   <- strsplit(bp, ";")[[1]]
        ccs   <- strsplit(cc, ";")[[1]]
        mfs   <- strsplit(mf, ";")[[1]]
        ids   <- strsplit(id, ";")[[1]]
        fulls <- strsplit(full, ";")[[1]]
        
        # Usa la longitud mínima de todos los vectores
        n <- min(length(genes), length(bps), length(ccs), length(mfs), length(ids), length(fulls))
        
        tibble::tibble(
          genename = genes[1:n], 
          GO_biological_process = bps[1:n],
          GO_cellular_component = ccs[1:n],
          GO_molecular_function = mfs[1:n],
          Pathway_KEGG_id = ids[1:n],
          Pathway_KEGG_full = fulls[1:n]
        )
      })) %>%
      dplyr::select(row) %>%
      unnest(row) %>%
      distinct()
    
    
    # Añadir los enlaces a genename y Pathway_KEGG_id
    kegg_pairs <- kegg_pairs %>%
      mutate(
        genename = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", genename, 
                          "' target='_blank'>", genename, "</a>"),
        Pathway_KEGG_id = paste0("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway+", Pathway_KEGG_id, 
                                 "' target='_blank'>", Pathway_KEGG_id, "</a>")
      )
    
    # Renderizar la tabla DT
    DT::datatable(kegg_pairs, escape = FALSE, rownames = FALSE, options = list(pageLength = 10))
  })
  
  ############################################################################################################
  
  output$kegg_wordcloud <- renderPlot({
    
    df <- keeg_input() 
    
    df <- df %>%
      filter(!is.na(Pathway_KEGG_full), Pathway_KEGG_full != "") %>%
      mutate(Pathway_KEGG_full = tolower(Pathway_KEGG_full)) %>%
      pull(Pathway_KEGG_full) %>%
      strsplit(";") %>%
      unlist() %>%
      trimws()
    
    # Contar frecuencia
    freq_table <- sort(table(df), decreasing = TRUE)
    
    # Wordcloud
    wordcloud(
      words = names(freq_table),
      freq = freq_table,
      scale = c(5, 0.8),
      min.freq = 2,
      max.words = 100,
      random.order = FALSE,
      rot.per = 0.15,
      colors = rainbow(10)
    )
  })
  
  output$kegg_barplot <- renderPlot({
    
    df <- keeg_input() 
    
    # Expandir los datos de pathways y genes
    genes_exp <- df %>%
      dplyr::select(genename) %>%
      filter(!is.na(genename)) %>%
      separate_rows(genename, sep = ";") %>%
      filter(genename != "") %>%
      mutate(row_id = row_number())
    
    pathways_exp <- df %>%
      dplyr::select(Pathway_KEGG_id, Pathway_KEGG_full) %>%
      filter(!is.na(Pathway_KEGG_id), !is.na(Pathway_KEGG_full)) %>%
      separate_rows(Pathway_KEGG_id, Pathway_KEGG_full, sep = ";") %>%
      filter(Pathway_KEGG_id != "", Pathway_KEGG_full != "") %>%
      mutate(row_id = row_number())
    
    # Combinar los datos de genes y pathways
    kegg_pairs <- inner_join(genes_exp, pathways_exp, by = "row_id") %>%
      dplyr::select(genename, Pathway_KEGG_id, Pathway_KEGG_full) %>%
      distinct() %>%
      filter(Pathway_KEGG_full != ".") # Excluir los pathways con "."
    
    # Agregar una columna para identificar el input
    kegg_pairs$input_source <- case_when(
      input$inputfile == "kg" ~ "KG",
      input$inputfile == "imputed" ~ "Imputed",
      TRUE ~ "Combined"
    )
    
    # Contar la frecuencia de cada pathway, agrupando por input_source
    df_freq <- kegg_pairs %>%
      count(Pathway_KEGG_full, input_source) %>%
      arrange(desc(n)) %>%
      slice_head(n = 20)  # Los 20 más frecuentes
    
    # Crear el gráfico de barras
    ggplot(df_freq, aes(x = reorder(Pathway_KEGG_full, n), y = n, fill = input_source)) +
      geom_col() +
      coord_flip() +
      labs(x = "Pathway", y = "Frecuencia", title = "Pathways KEGG más frecuentes por fuente de datos") +
      theme_minimal() +
      theme(legend.title = element_blank())  # Eliminar el título de la leyenda
  })
  
  ###############################################################################
  output$manhat_plot <- renderPlot({
  
  cr.PIR <- calc_candidate_regions(ihs_PIR,
                                   threshold = 4,
                                   pval = TRUE,
                                   window_size = 1E6,
                                   overlap = 1E5,
                                   min_n_extr_mrk = 2)
  cr.PIR
  
  cr.CEU <- calc_candidate_regions(ihs_CEU,
                                   threshold = 4,
                                   pval = TRUE,
                                   window_size = 1E6,
                                   overlap = 1E5,
                                   min_n_extr_mrk = 2)
  cr.CEU
  
  cr.IBS <- calc_candidate_regions(ihs_IBS,
                                   threshold = 4,
                                   pval = TRUE,
                                   window_size = 1E6,
                                   overlap = 1E5,
                                   min_n_extr_mrk = 2)
  cr.IBS
  
  ####################### function to plot manhattan by candidate regions
  manhattanplot <- function(df, cr = NULL, title = "") {
    df <- df %>%
      filter(!is.na(CHR), !is.na(POSITION), !is.na(LOGPVALUE)) %>%
      mutate(CHR = as.numeric(CHR),
             P = 10^(-LOGPVALUE)) %>%
      arrange(CHR, POSITION) %>%
      filter(P < 0.05)
    
    chr_lengths <- df %>%
      group_by(CHR) %>%
      dplyr::summarize(chr_len = max(POSITION), .groups = "drop") %>%
      mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len))
    
    df <- df %>%
      left_join(chr_lengths, by = "CHR") %>%
      mutate(BPcum = POSITION + tot)
    
    axis_df <- df %>%
      group_by(CHR) %>%
      dplyr::summarize(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop")
    
    p <- ggplot() +
      geom_point(data = df,
                 aes(x = BPcum, y = -log10(P), color = as.factor(CHR %% 2)),
                 alpha = 0.75, size = 1) +
      scale_color_manual(values = c("red", "green")) +
      scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
      geom_hline(yintercept = -log10(1e-4), linetype = "dashed", color = "blue") +
      labs(x = "Chromosome", y = "-log10(p-value)", title = title) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
      )
    
    if (!is.null(cr)) {
      cr <- cr %>%
        left_join(df %>% group_by(CHR) %>% dplyr::summarize(chr_offset = min(BPcum)), by = "CHR") %>%
        mutate(BPcum = START + chr_offset)
      
      # Añadir las líneas coloreadas sin tocar el mapeo global
      p <- p + 
        geom_vline(data = cr, 
                   aes(xintercept = BPcum),
                   color = cr$color,  # directamente desde el vector
                   linetype = "dashed", linewidth = 0.5)
    }
    
    return(p)
  }
  
  
  ihs_PIR2 <- ihs_PIR %>% dplyr::select(CHR, POSITION, LOGPVALUE)
  candidate_regions <- data.frame(CHR = cr.PIR$CHR, START = cr.PIR$START)
  candidate_regions$color <- c("black","black","black","grey","red","black","red","grey","black","black","grey","black","red","grey","red","red","black","red","grey")
  
  mhtt.PIR <- manhattanplot(ihs_PIR2, cr = candidate_regions, title = "PIR")
  
  ihs_CEU2 <- ihs_CEU %>% dplyr::select(CHR, POSITION, LOGPVALUE)
  candidate_regions <- data.frame(CHR = cr.CEU$CHR, START = cr.CEU$START)
  candidate_regions$color <- c("black","black","black","black","black","blue","black","grey","grey","blue","black","blue","grey")
  
  mhtt.CEU <- manhattanplot(ihs_CEU2, cr = candidate_regions, title = "CEU")
  
  ihs_IBS2 <- ihs_IBS %>% dplyr::select(CHR, POSITION, LOGPVALUE)
  candidate_regions <- data.frame(CHR = cr.IBS$CHR, START = cr.IBS$START)
  candidate_regions$color <- c("black","green","black","black","grey","green","green","black","grey","black","grey","green","black","grey","green","green","green","black")
  mhtt.IBS <- manhattanplot(ihs_IBS2, cr = candidate_regions, title = "IBS")
  
  # Mostrar juntos
  library(patchwork)
  mhtt.CEU / mhtt.IBS / mhtt.PIR + plot_layout(ncol = 1)
  })
  
  ################################################################################
}

# Run the application 
shinyApp(ui = ui, server = server)
