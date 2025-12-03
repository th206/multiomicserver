rm(list = ls())
options(shiny.maxRequestSize = 1024 * 1024 * 1024)

library(shiny)
library(shinyjs)
library(shinydashboard)
library(readxl)
library(dplyr)
library(tidyr)
library(visNetwork)
library(DT)
library(shinyWidgets)
library(igraph)
library(rlang)
library(writexl)
library(digest)

# ====================================================
# Helper Functions for Network Operations
# ====================================================

neighbors <- function(nodes, network) {
  nodes <- tolower(nodes)
  ix <- unlist(lapply(nodes, function(x) {
    union(
      which(tolower(network$edges$from_id) == x),
      which(tolower(network$edges$to_id) == x)
    )
  })) %>% unique()
  unique(c(network$edges$from_id[ix], network$edges$to_id[ix]))
}

stable_pick <- function(key, choices) {
  h <- digest(key, algo = "xxhash32")
  i <- strtoi(substr(h, 1, 7), base = 16L)
  choices[(i %% length(choices)) + 1L]
}

maxneighbors <- function(nodes, network, limit = 0) {
  if (limit == 0) limit <- 1E99
  nnodes     <- length(nodes)
  nnodeslast <- 0
  nodeslast  <- nodes
  
  while (nnodes > nnodeslast && nnodes < limit) {
    nodeslast   <- nodes
    nnodeslast  <- nnodes
    nodes       <- neighbors(nodes, network)
    nnodes      <- length(nodes)
    cat("maxneighbors iteration, nodes:", nnodes, "\n")
  }
  
  selected_nodes <- unique(nodeslast)
  
  if (length(selected_nodes) > limit) {
    cat("maxneighbors trimming nodes from", length(selected_nodes), "to", limit, "\n")
    
    node_scores <- data.frame(id = selected_nodes) %>%
      mutate(degree = sapply(id, function(x) {
        sum(network$edges$from_id == x | network$edges$to_id == x, na.rm = TRUE)
      })) %>%
      arrange(desc(degree))
    
    selected_nodes <- head(node_scores$id, limit)
  }
  
  cat("maxneighbors final nodes:", length(selected_nodes), "\n")
  selected_nodes
}

maxneighbors_noSTAT <- function(nodes, network, limit = 0) {
  if (limit == 0) limit <- 1E99
  nodesin    <- nodes
  nnodes     <- length(nodes)
  nnodeslast <- 0
  nodeslast  <- nodes
  
  while ((nnodes > nnodeslast) & (nnodes <= limit)) {
    nodeslast    <- nodes
    nnodeslast   <- nnodes
    nogrow_nodes <- nodes[grep("STAT: |GWAS: ", nodes)]
    grow_nodes   <- nodes[grep("STAT: |GWAS: ", nodes, invert = TRUE)]
    nodes        <- c(nogrow_nodes, neighbors(grow_nodes, network))
    nnodes       <- length(nodes)
  }
  
  unique(c(nodesin, nodeslast))
}

nodes2network <- function(nodes, network) {
  nodes <- tolower(nodes)
  cat("Input nodes to nodes2network:", paste(head(nodes, 5), collapse = ", "), "\n")
  
  ix_from <- unique(unlist(lapply(nodes, function(x)
    which(tolower(network$edges$from_id) == x))))
  ix_to   <- unique(unlist(lapply(nodes, function(x)
    which(tolower(network$edges$to_id) == x))))
  ix      <- unique(c(ix_from, ix_to))
  
  iy <- unique(unlist(lapply(nodes, function(x)
    which(tolower(network$nodes$id) == x |
            tolower(network$nodes$TRAITID) == x))))
  
  cat("Nodes matched in network$nodes:", length(iy), "\n")
  
  subnet <- list(
    edges = network$edges[ix, ],
    nodes = network$nodes[iy, ]
  )
  
  subnet$edges <- subnet$edges[
    subnet$edges$from_id %in% subnet$nodes$id &
      subnet$edges$to_id %in% subnet$nodes$id,
  ]
  
  cat("Subnet nodes:", nrow(subnet$nodes), "Edges:", nrow(subnet$edges), "\n")
  subnet
}

normalize_headers <- function(df) {
  nm <- names(df)
  keymap <- c(
    "traitid1"  = "TRAITID1", "traitid2"  = "TRAITID2",
    "pvalue"    = "PVALUE",   "beta"      = "BETA",    "cor"  = "COR",
    "id"        = "id",
    "traitid"   = "TRAITID",  "shortname" = "SHORTNAME", "plat" = "PLAT",
    "color"     = "COLOR",    "shape"     = "SHAPE"
  )
  for (i in seq_along(nm)) {
    k <- tolower(nm[i]); if (k %in% names(keymap)) nm[i] <- keymap[[k]]
  }
  names(df) <- nm
  df
}

processExcelSheet <- function(inFile, sheet, mand_cols, anno_cols, custom_cols) {
  data <- readxl::read_excel(path = inFile, sheet = sheet, .name_repair = "unique_quiet")
  data <- normalize_headers(data)
  
  has_traitid1 <- "TRAITID1" %in% names(data)
  has_traitid2 <- "TRAITID2" %in% names(data)
  validNetwork <- has_traitid1 && has_traitid2
  
  if (validNetwork) {
    has_p <- "PVALUE" %in% names(data)
    has_w <- any(c("BETA", "COR") %in% names(data))
    
    if (!has_p && has_w)  data$PVALUE <- 0.05
    if (has_p && !has_w)  data$BETA   <- 1
    if (!has_p && !has_w) { data$PVALUE <- 1; data$BETA <- 1 }
  }
  
  validAnno   <- all(toupper(anno_cols)   %in% names(data))
  validCustom <- all(toupper(custom_cols) %in% names(data))
  
  list(
    sheet        = sheet,
    data         = data,
    validNetwork = validNetwork,
    validAnno    = validAnno,
    validCustom  = validCustom
  )
}

processExcelFile <- function(inFile) {
  sheets <- readxl::excel_sheets(inFile)
  
  mand_cols   <- c("TRAITID1", "TRAITID2", "PVALUE", "BETA", "COR")
  anno_cols   <- c("TRAITID", "SHORTNAME", "PLAT")
  custom_cols <- c("PLAT", "COLOR", "SHAPE")
  
  datalist    <- list()
  sheets_ok   <- character()
  sheets_flag <- character()
  anno        <- NULL
  custom      <- NULL
  
  for (sheet in sheets) {
    res <- processExcelSheet(inFile, sheet, mand_cols, anno_cols, custom_cols)
    
    if (res$validNetwork) {
      data <- res$data
      
      wcol   <- intersect(c("BETA", "COR"), names(data))
      weight <- if (length(wcol)) data[[wcol[1]]] else rep(1, nrow(data))
      pval   <- data[["PVALUE"]]; if (is.null(pval)) pval <- rep(0.05, nrow(data))
      
      if (!"id" %in% names(data)) data$id <- paste0(sheet, "_", seq_len(nrow(data)))
      
      df <- data.frame(
        to     = trimws(data[["TRAITID1"]]),
        from   = trimws(data[["TRAITID2"]]),
        pvalue = as.numeric(pval),
        weight = as.numeric(weight),
        id     = data[["id"]],
        type   = sheet,
        stringsAsFactors = FALSE
      )
      
      datalist[[sheet]] <- df
      sheets_ok <- c(sheets_ok, sheet)
      
    } else if (res$validAnno) {
      
      anno <- res$data |>
        dplyr::mutate(TRAITID = trimws(TRAITID)) |>
        dplyr::select(dplyr::all_of(anno_cols))
      
    } else if (res$validCustom) {
      
      custom <- res$data |>
        dplyr::select(dplyr::all_of(custom_cols))
      
    } else {
      
      need <- c("TRAITID1", "TRAITID2")
      missing_cols <- setdiff(need, names(res$data))
      if (length(missing_cols) == 0)
        missing_cols <- "TRAITID1/TRAITID2 not detected (check spelling/case)."
      
      sheets_flag <- c(
        sheets_flag,
        sprintf(
          "[%s] Missing required columns: %s",
          sheet, paste(missing_cols, collapse = ", ")
        )
      )
    }
  }
  
  if (length(datalist) == 0) {
    stop("No valid network sheets found. Each network sheet must contain TRAITID1 and TRAITID2 (case-insensitive).")
  }
  
  if (is.null(anno)) {
    all_traits <- unique(unlist(lapply(
      datalist,
      function(df) c(as.character(df$to), as.character(df$from))
    )))
    anno <- data.frame(
      TRAITID   = all_traits,
      SHORTNAME = all_traits,
      PLAT      = rep("Unknown", length(all_traits)),
      stringsAsFactors = FALSE
    )
  }
  
  list(
    datalist    = datalist,
    sheets_ok   = sheets_ok,
    sheets_flag = sheets_flag,
    anno        = anno,
    custom      = custom
  )
}

constructNetwork <- function(datalist, anno, custom) {
  all_edges <- do.call(rbind, datalist)
  cat(
    "Initial edges:", nrow(all_edges),
    "Sample from/to:", head(all_edges$from, 3), "->", head(all_edges$to, 3), "\n"
  )
  
  all_edges$sign   <- sign(all_edges$weight)
  all_edges$weight <- abs(all_edges$weight)
  all_edges$pvalue <- as.numeric(as.character(all_edges$pvalue))
  
  if (any(is.na(all_edges$pvalue))) {
    warning("Some p-values could not be converted to numeric; setting them to 1E-100.")
    all_edges$pvalue[is.na(all_edges$pvalue)] <- 1E-100
  }
  
  all_edges$pvalue[all_edges$pvalue == 0] <- 1E-100
  
  edge_thickness <- -log10(all_edges$pvalue)
  
  if (min(edge_thickness) == max(edge_thickness)) {
    all_edges$width <- 1
  } else {
    seq_range <- round(seq(
      min(edge_thickness) - 1,
      max(edge_thickness) - 1,
      length.out = 10
    ))
    if (length(unique(seq_range)) < length(seq_range)) {
      seq_range <- unique(seq_range)
    }
    all_edges$width <- as.numeric(cut(
      edge_thickness,
      breaks = seq_range,
      labels = FALSE
    ))
  }
  
  all_nodes <- data.frame(
    TRAITID = unique(sort(c(all_edges$from, all_edges$to, anno$TRAITID))),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(anno)) {
    anno <- anno %>% distinct(TRAITID, .keep_all = TRUE)
    all_nodes <- left_join(all_nodes, anno, by = "TRAITID")
  } else {
    stop("Missing annotation. Please include a sheet with headers: TRAITID, SHORTNAME, PLAT")
  }
  
  all_nodes <- all_nodes %>%
    mutate(
      SHORTNAME = ifelse(is.na(SHORTNAME), TRAITID, SHORTNAME),
      PLAT      = ifelse(is.na(PLAT), "Unknown", PLAT)
    )
  
  all_edges$from_id <- all_edges$from
  all_edges$to_id   <- all_edges$to
  
  all_edges$from <- all_nodes$SHORTNAME[
    match(tolower(all_edges$from), tolower(all_nodes$TRAITID))
  ]
  all_edges$to <- all_nodes$SHORTNAME[
    match(tolower(all_edges$to), tolower(all_nodes$TRAITID))
  ]
  
  all_edges <- all_edges[!is.na(all_edges$from) & !is.na(all_edges$to), ]
  
  cat(
    "Edges after mapping:", nrow(all_edges),
    "Vitamin_D edges:",
    sum(
      tolower(all_edges$from_id) == "vitamin_d" |
        tolower(all_edges$to_id) == "vitamin_d"
    ),
    "\n"
  )
  
  all_nodes$id <- all_nodes$SHORTNAME
  names(all_nodes)[names(all_nodes) == "PLAT"] <- "plat"
  
  if (!is.null(custom)) {
    all_platforms     <- union(unique(all_nodes$plat), custom$PLAT)
    missing_platforms <- setdiff(all_platforms, custom$PLAT)
    
    missing_df <- data.frame(
      PLAT  = missing_platforms,
      COLOR = rep("#999999", length(missing_platforms)),
      SHAPE = rep("triangle", length(missing_platforms)),
      stringsAsFactors = FALSE
    )
    
    custom <- rbind(custom, missing_df)
    
    plat_colors <- setNames(custom$COLOR, custom$PLAT)
    plat_shapes <- setNames(custom$SHAPE, custom$PLAT)
    
    all_nodes$color <- plat_colors[all_nodes$plat]
    all_nodes$shape <- plat_shapes[all_nodes$plat]
    
  } else {
    plat_list <- unique(all_nodes$plat)
    
    default_colors <- setNames(rep("#999999", length(plat_list)), plat_list)
    default_shapes <- setNames(rep("star", length(plat_list)),   plat_list)
    
    default_colors[c(
      "DNA", "SOMA", "ALAMAR", "BRAIN", "BM", "LD", "CPG", "CLIN",
      "CM", "IgA", "RNA", "HD4", "IgG", "miRNA", "OLINK", "PGP",
      "PM", "SM", "UM", "STAT", "GWAS"
    )] <- c(
      "#23bbee", "#a62281", "#a62281", "#f2921f", "#ffc815", "#ffc815",
      "#145da9", "#a0b6a8", "#57ba47", "#e41d30", "#5c2d83", "#57ba47",
      "#e41d30", "#5c2d83", "#a62281", "#e41d30", "#57ba47", "#57ba47",
      "#57ba47", "#EEEEEE", "#EEEEEE"
    )
    
    default_shapes[c("UM", "CM", "SM", "DNA", "RNA", "CPG", "STAT", "GWAS")] <- c(
      "square", "square", "triangle", "diamond",
      "diamond", "diamond", "circle", "dot"
    )
    
    all_nodes$color <- default_colors[all_nodes$plat]
    all_nodes$shape <- default_shapes[all_nodes$plat]
  }
  
  all_nodes$group <- all_nodes$plat
  
  all_edges$title <- paste0(
    all_edges$type, ": ", all_edges$id,
    ", p=", all_edges$pvalue,
    ", beta=", ifelse(all_edges$sign > 0, "", "-"), all_edges$weight
  )
  
  all_edges$color <- ifelse(all_edges$sign > 0, "blue", "red")
  
  list(edges = all_edges, nodes = all_nodes)
}

# ====================================================
# Server Function
# ====================================================

server <- function(input, output, session) {
  # Define persistent storage paths
  data_file <- "data.rds"
  custom_settings_file <- "custom_settings.rds"
  traitLabelTrigger <- reactiveVal(0)
  labelTrigger <- reactiveVal(0)
  rv <- reactiveValues(
    customSelections = list(),
    useCases = list(),
    dataLoaded = FALSE,
    selectedTrait = NULL,
    useCaseActive = FALSE,
    triggerNetworkAfterViewLoad = FALSE,
    pendingFocus = NULL
  )
  networkTrigger <- reactiveVal(0)
  processClicked <- reactiveVal(FALSE)
  myValues <- reactiveValues(
    datalist = NULL,
    sheets_ok = NULL,
    sheets_flag = NULL,
    anno = NULL,
    custom = NULL,
    fullnet = NULL,
    platlist = NULL,
    listItems = NULL,
    allSheets = NULL
  )
  observeEvent(session, {
    updateTabItems(session, "tabs", "Tab6")
  }, once = TRUE)
  
  observe({
    if (isTRUE(processClicked())) {
      # Hide Data Formats (Tab4)
      shinyjs::hide(selector = "a[data-value='Tab4']")
      # Hide Import Data (Tab7)
      shinyjs::hide(selector = "a[data-value='Tab7']")
    } else {
      # Show them when Process has NOT been clicked (or after clear)
      shinyjs::show(selector = "a[data-value='Tab4']")
      shinyjs::show(selector = "a[data-value='Tab7']")
    }
  })
  
  # After a view is loaded, trigger the network in a clean reactive context
  observe({
    if (isTRUE(rv$triggerNetworkAfterViewLoad)) {
      # reset flag first to avoid loops
      rv$triggerNetworkAfterViewLoad <- FALSE
      
      # Now safe to bump the trigger and switch tabs
      networkTrigger(networkTrigger() + 1)
      updateTabItems(session, "tabs", "Tab1")
    }
  })
  
  
  ############## New ############
  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  # Dynamically show/disable import controls on Tab 7, based on rv$dataLoaded
  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  output$importControls <- renderUI({
    # If data has NOT yet been loaded/processed, show everything as normal:
    if (!processClicked()) {
      tagList(
        tags$head(tags$style(HTML("
        .page-wrap { max-width: 1200px; margin: 0 auto; }
        /* Make inputs/selects full width inside the box */
        .box .form-control,
        .box .selectize-control,
        .box .selectize-input,
        .box .input-group { width: 100% !important; }
        /* DataTable fills width */
        #dynamicTable { width: 100% !important; }
        /* Tighter look on small screens */
        @media (max-width: 768px){
          .page-wrap { padding: 0 12px; }
        }
        .box {
         box-shadow: 0 2px 6px rgba(0,0,0,0.08);
        border-radius: 10px;
        }
      "))),
      
      div(class = "page-wrap",
          
          fluidRow(
            column(width = 2,
                   selectInput("skin_color", "Change Skin Color", 
                               c("red", "purple", "green", "blue", "yellow"))
            ),
            column(width = 2,
                   textInput("titleInput", "Enter Server Title:", "Multi-omics Server")
            )
          ),
          fluidRow(
            
            # LEFT: Upload + summaries + preview
            column(
              width = 7,
              box(
                width = 12, solidHeader = TRUE, status = "primary",
                title = "Upload Data to Build a Network",
                
                # Guidance (kept compact on the left)
                div(style = "background:#f7f7fb;border-left:4px solid #5c6bc0;padding:10px 12px;margin-bottom:10px;border-radius:6px;",
                    tags$b("What file do I need?"), br(),
                    "Each row represents a relationship between two items.",
                    tags$ul(
                      tags$li(tags$b("Two columns, TRAITID1, and TRAITID2"), " are the linked nodes."),
                      tags$li("Optional: ", tags$code("PVALUE"), ", ", tags$code("BETA"), ", ", tags$code("COR")),
                      tags$li("Accepted: ", tags$code(".xlsx") )
                    )
                ),
                
                fileInput('file_in', 'Choose Excel file',
                          accept = c('.xlsx', '.xls'), multiple = FALSE),
                
                uiOutput("SHEETlist"),        
                uiOutput("summaryStats"),     
                textOutput("summary1"),       
                textOutput("summary3"),
                uiOutput("summary2"),
                
                hr(),
                DTOutput("dynamicTable"),
                hr(),
                
                fluidRow(
                  column(6, uiOutput("confirm_btn")),
                  column(6, uiOutput("downloadButtonUI"))
                ),
                
                #actionButton("clearData", "Clear Saved Data", class = "clear-data-btn")
              )
            ),
            # RIGHT: example + custom home page
            column(
              width = 5,
              box(
                width = 12, solidHeader = TRUE, status = "info",
                title = "Data Format Guide",
                
                tags$p("Your file should contain at least two columns defining pairwise links."),
                tags$ol(
                  tags$li(tags$b("TRAITID1 and TRAITID2 = nodes"), " (e.g., ", tags$code("GeneA, GeneB"), ")."),
                  tags$li("Optional columns (used for filtering/visuals): ",
                          tags$code("PVALUE"), ", ", tags$code("BETA"), ", ", tags$code("COR"), "."),
                  tags$li("Optional Sheet name: ", tags$code("ANNO"))
                ),
                tags$hr(),
                tags$p(tags$em("Examples: gene–gene (PPI), metabolite–enzyme, miRNA–target.")),
                
                # Example file download
                downloadButton("download_example", "Download Example", class = "btn btn-default",
                               style = "margin-top:6px;")
              ),
              
              box(
                width = 12, solidHeader = TRUE, status = "warning",
                title = "Custom Home Page (HTML/Markdown)",
                fileInput("home_upload", "Upload HTML/MD", accept = c(".html", ".md")),
                tags$p(style="color:#666; font-size:13px; margin:0;",
                       "If provided, it will be shown as the Home tab.")
              )
            )
          ),
          fluidRow(
            column(width = 2,
                   uiOutput("mySelection")
            )
          )
      )
      )
    } else {
      # If data HAS already been loaded/processed, replace these controls
      # with disabled versions (or simply hide them entirely). For example:
      tagList(
        fluidRow(
          column(width = 12,
                 # Show a message instead of the inputs:
                 tags$div(
                   style = "padding: 10px; background-color: #f0f0f0; border-radius: 5px;",
                   "Data has already been loaded. To upload a new file, please click 'Clear Saved Data' first."
                 )
          )
        )
      )
    }
  })
  
  ####################################

  output$download_example <- downloadHandler(
    filename = function() {
      "example_data.xlsx"   # name seen by the user
    },
    content = function(file) {
      # path to the file inside your app
      example_path <- "www/string_interactions_short.xlsx"
      
      # copy to the temp file `file`
      file.copy(from = example_path, to = file, overwrite = TRUE)
    }
  )
  
  
  # Log package versions
  observe({
    cat("Shiny version:", as.character(packageVersion("shiny")), "\n")
    cat("ShinyWidgets version:", as.character(packageVersion("shinyWidgets")), "\n")
  })
  
  loadData <- function(saved_data, saved_settings) {
    required_node_cols <- c("TRAITID", "id", "plat", "color", "shape")
    required_edge_cols <- c("from_id", "to_id", "from", "to")
    if (!is.null(saved_data$datalist) && 
        !is.null(saved_data$anno) && 
        !is.null(saved_data$fullnet) &&
        all(required_node_cols %in% names(saved_data$fullnet$nodes)) &&
        all(required_edge_cols %in% names(saved_data$fullnet$edges))) {
      isolate({
        myValues$datalist <- saved_data$datalist
        myValues$sheets_ok <- saved_data$sheets_ok
        myValues$allSheets <- saved_data$sheets_ok
        myValues$sheets_flag <- saved_data$sheets_flag
        myValues$anno <- saved_data$anno
        myValues$custom <- saved_data$custom
        myValues$fullnet <- saved_data$fullnet
        myValues$platlist <- c("ALL", sort(unique(saved_data$fullnet$nodes$plat)))
        rv$customSelections <- saved_settings$customSelections
        rv$dataLoaded <- TRUE
      })
      
      isolate({
        for (platform in names(rv$customSelections)) {
          if (!is.null(rv$customSelections[[platform]]$color)) {
            myValues$fullnet$nodes$color[myValues$fullnet$nodes$plat == platform] <- 
              rv$customSelections[[platform]]$color
          }
          if (!is.null(rv$customSelections[[platform]]$shape)) {
            myValues$fullnet$nodes$shape[myValues$fullnet$nodes$plat == platform] <- 
              rv$customSelections[[platform]]$shape
          }
        }
      })
      
      isolate({
        updateSelectInput(session, "plat", choices = myValues$platlist, selected = myValues$platlist[1])
        if (!is.null(myValues$fullnet)) {
          traitlist <- sort(myValues$fullnet$nodes$id)
          cat("loadData: Trait list length:", length(traitlist), "\n")
          cat("loadData: Sample traits:", paste(head(traitlist, 5), collapse = ", "), "\n")
          cat("loadData: Last few traits:", paste(tail(traitlist, 5), collapse = ", "), "\n")
          tryCatch({
            updateSelectInput(
              session,
              "trait",
              choices = traitlist,
              selected = traitlist[1]
            )
            rv$selectedTrait <- traitlist[1]
          }, error = function(e) {
            cat("Error updating trait dropdown in loadData:", e$message, "\n")
            showNotification(paste("Error updating trait dropdown:", e$message), type = "error")
          })
        }
      })
      
      cat("loadData: Platforms in fullnet$nodes:", paste(unique(myValues$fullnet$nodes$plat), collapse = ", "), "\n")
      cat("loadData: Total nodes:", nrow(myValues$fullnet$nodes), "\n")
      
      saveRDS(list(
        datalist = myValues$datalist,
        sheets_ok = myValues$sheets_ok,
        sheets_flag = myValues$sheets_flag,
        anno = myValues$anno,
        custom = myValues$custom,
        fullnet = myValues$fullnet
      ), data_file)
      saveRDS(list(customSelections = rv$customSelections), custom_settings_file)
      cat("Saved data to", data_file, "and settings to", custom_settings_file, "\n")
      
      updateTabItems(session, "tabs", "Tab2")
      cat("Loaded saved data and settings successfully\n")
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  observe({
    req(myValues$fullnet, input$plat)
    if (input$plat == "ALL") {
      traitlist <- sort(myValues$fullnet$nodes$id)
    } else if (input$plat == "Unknown") {
      ix <- which(myValues$fullnet$nodes$plat == "Unknown" | is.na(myValues$fullnet$nodes$plat))
      traitlist <- sort(myValues$fullnet$nodes$id[ix])
    } else {
      ix <- which(myValues$fullnet$nodes$plat == input$plat)
      traitlist <- sort(myValues$fullnet$nodes$id[ix])
    }
    if (length(traitlist) == 0) {
      traitlist <- c("No traits available")
    }
    cat("Updating trait list for platform", input$plat, ":", length(traitlist), "traits\n")
    tryCatch({
      updateSelectInput(
        session,
        inputId = "trait",
        choices = traitlist,
        selected = traitlist[1]
      )
    }, error = function(e) {
      cat("Error updating trait dropdown:", e$message, "\n")
      showNotification(paste("Error updating trait dropdown:", e$message), type = "error")
    })
  })
  

  # Clear data
  observeEvent(input$clearData, {
    if (file.exists(data_file)) file.remove(data_file)
    if (file.exists(custom_settings_file)) file.remove(custom_settings_file)
    rv$dataLoaded <- FALSE
    myValues$datalist <- NULL
    myValues$fullnet <- NULL
    myValues$sheets_ok <- NULL
    myValues$allSheets <- NULL
    myValues$sheets_flag <- NULL
    myValues$anno <- NULL
    myValues$custom <- NULL
    myValues$platlist <- NULL
    myValues$listItems <- NULL
    rv$customSelections <- list()
    processClicked(FALSE)
    showNotification("Saved data cleared.", type = "message")
    updateTabItems(session, "tabs", "Tab7")
  })
  
  # Process Excel file
  processedData <- reactive({
    req(input$file_in)
    if (file.exists(data_file)) file.remove(data_file)
    if (file.exists(custom_settings_file)) file.remove(custom_settings_file)
    
    tryCatch({
      result <- processExcelFile(input$file_in$datapath)
      if(length(result$datalist) == 0){
        sendSweetAlert(
          session = session,
          title = "Missing Mandatory Sheet",
          text = "The uploaded file does not contain any sheet with the mandatory headers: TRAITID1, TRAITID2, PVALUE, and one of BETA or COR.",
          type = "error"
        )
        return(NULL)
      }
      if(is.null(result$anno)) {
        sendSweetAlert(
          session = session,
          title = "Missing Annotation",
          text = "The uploaded file is missing the annotation sheet with headers: TRAITID, SHORTNAME, PLAT.",
          type = "error"
        )
        return(NULL)
      }
      myValues$datalist <- result$datalist
      myValues$sheets_ok <- result$sheets_ok
      myValues$allSheets <- result$sheets_ok
      myValues$sheets_flag <- result$sheets_flag
      myValues$anno <- result$anno
      myValues$custom <- result$custom
      rv$dataLoaded <- TRUE
      cat("Uploaded datalist sheets:", paste(names(result$datalist), collapse = ", "), "\n")
      result
    }, error = function(e) {
      sendSweetAlert(
        session = session,
        title = "Error processing file",
        text = paste("Error:", e$message),
        type = "error"
      )
      return(NULL)
    })
  })
  
  observeEvent(processedData(), {
    dataResult <- processedData()
    if (is.null(dataResult)) return(NULL)
    fullnet <- tryCatch({
      constructNetwork(dataResult$datalist, dataResult$anno, dataResult$custom)
    }, error = function(e) {
      sendSweetAlert(
        session = session,
        title = "Error constructing network",
        text = paste("Error:", e$message),
        type = "error"
      )
      return(NULL)
    })
    
    if (is.null(fullnet)) return(NULL)
    myValues$fullnet <- fullnet
    plat_list <- c("ALL", sort(unique(fullnet$nodes$plat)))
    myValues$platlist <- plat_list
    cat("Platforms in fullnet$nodes:", paste(unique(fullnet$nodes$plat), collapse = ", "), "\n")
    cat("Total nodes:", nrow(fullnet$nodes), "\n")
    vitamin_d_rows <- fullnet$nodes %>% dplyr::filter(id == "Vitamin_D" | TRAITID == "Vitamin_D")
    cat("Vitamin_D rows in fullnet$nodes:\n")
    print(vitamin_d_rows)
    isolate({
      updateSelectInput(session, "plat", choices = plat_list, selected = plat_list[1])
      traitlist <- sort(fullnet$nodes$id)
      cat("Initial trait list length:", length(traitlist), "\n")
      cat("Initial sample traits:", paste(head(traitlist, 5), collapse = ", "), "\n")
      cat("Initial last few traits:", paste(tail(traitlist, 5), collapse = ", "), "\n")
      cat("Initial is Vitamin_D in traitlist?", "Vitamin_D" %in% traitlist, "\n")
      tryCatch({
        updateSelectInput(
          session,
          "trait",
          choices = traitlist,
          selected = traitlist[1]
        )
        rv$selectedTrait <- traitlist[1]
      }, error = function(e) {
        cat("Error updating trait dropdown in processedData:", e$message, "\n")
        showNotification(paste("Error updating trait dropdown:", e$message), type = "error")
      })
    })

    
    output$summaryStats <- renderUI({
      # --- Safely read data ---
      nodes <- tryCatch(myValues$fullnet$nodes, error = function(e) NULL)
      edges <- tryCatch(myValues$fullnet$edges, error = function(e) NULL)
      
      # --- If missing, show gentle message ---
      if (is.null(nodes) || is.null(edges) || !"TRAITID" %in% names(nodes)) {
        return(div(style = "color:#6c757d; margin-top:6px;", "Waiting for data…"))
      }
      
      # --- Compute numbers ---
      trait_n <- format(length(unique(nodes$TRAITID)), big.mark = ",")
      assoc_n <- format(nrow(edges), big.mark = ",")
      
      # --- Styled dual metric display ---
      tags$div(
        style = "display:flex; gap:10px; margin-top:8px;",
        tags$div(
          style = "flex:1; background:#e3f2fd; border-left:4px solid #1565C0;
               padding:8px 12px; border-radius:8px; text-align:center;",
          HTML(sprintf("<b>Unique TRAIT IDs</b><br><span style='font-size:20px; color:#0d47a1;'>%s</span>", trait_n))
        ),
        tags$div(
          style = "flex:1; background:#e3f2fd; border-left:4px solid #42A5F5;
               padding:8px 12px; border-radius:8px; text-align:center;",
          HTML(sprintf("<b>Total Associations</b><br><span style='font-size:20px; color:#0d47a1;'>%s</span>", assoc_n))
        )
      )
    })
    
    
    observe({
      # only run once fullnet and its nodes exist
      req(myValues$fullnet)
      req(myValues$fullnet$nodes)
      
      rv$platform_counts <- myValues$fullnet$nodes %>%
        dplyr::group_by(plat) %>%
        dplyr::summarise(
          trait_count = dplyr::n_distinct(id),   # or TRAITID if you prefer
          .groups = "drop"
        ) %>%
        dplyr::rename(PLAT = plat)
    })
  })
  
  # Dynamic trait dropdown
  output$SelectTrait <- renderUI({
    req(myValues$fullnet)
    cat("Rendering SelectTrait, input$plat:", input$plat, "\n")
    trait_choices <- myValues$fullnet$nodes$id
    selectInput("trait", "Select Trait:", choices = c("", trait_choices), selected = NULL)
  })
  
  
  observeEvent(input$trait, {
    req(myValues$fullnet)
    cat("observeEvent input$trait triggered, input$trait:", input$trait, ", rv$selectedTrait:", rv$selectedTrait, ", input$current_node_selection:", input$current_node_selection, "\n")
    if (is.null(input$trait) || input$trait == "" || (!is.null(rv$selectedTrait) && input$trait == rv$selectedTrait)) {
      cat("Skipping input$trait update: empty, unchanged, or invalid, input$trait:", input$trait, ", rv$selectedTrait:", rv$selectedTrait, "\n")
      return()
    }
    cat("Processing input$trait:", input$trait, "\n")
    selected_traitid <- myValues$fullnet$nodes$TRAITID[
      tolower(myValues$fullnet$nodes$id) == tolower(input$trait)
    ]
    if (length(selected_traitid) == 0) {
      cat("No TRAITID found for input$trait:", input$trait, "\n")
      showNotification("Selected trait not found in network.", type = "error")
      return()
    }
    rv$selectedTrait <- input$trait
    rv$pendingFocus <- selected_traitid
    storage$focus <- NULL # Clear current focus
    storage$subnet <- NULL # Clear network
    storage$nodes <- NULL
    storage$edges <- NULL
    rv$useCaseActive <- FALSE # Reset use case
    # Clear node selection to prioritize dropdown
    session$sendCustomMessage(type = "clearNodeSelection", message = list())
    session$sendCustomMessage(type = "forceClearNodeSelection", message = list())
    # Add a slight delay to ensure client-side update
    session$sendCustomMessage(type = "forceClearNodeSelectionDelayed", message = list(delay = 100))
    cat("Set rv$selectedTrait to:", input$trait, ", rv$pendingFocus to TRAITID:", selected_traitid, 
        ", cleared storage$focus, storage$subnet, storage$nodes, storage$edges, rv$useCaseActive, and input$current_node_selection\n")
    # Trigger UI update
    networkTrigger(networkTrigger() + 1)
  })
  
  
  output$selectionStatus <- renderText({
    if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      paste("Selected Node:", input$current_node_selection)
    } else if (!is.null(rv$selectedTrait)) {
      paste("Selected Trait:", rv$selectedTrait)
    } else {
      "No selection"
    }
  })
  
  
  # Dynamic button label for Update Network
  output$updateNetworkLabel <- renderText({
    cat("Rendering updateNetworkLabel, input$current_node_selection:", input$current_node_selection, 
        ", input$trait:", input$trait, ", rv$selectedTrait:", rv$selectedTrait, 
        ", rv$pendingFocus:", rv$pendingFocus, ", storage$focus:", storage$focus, "\n")
    # Validate myValues$fullnet$nodes
    nodes_valid <- is.data.frame(myValues$fullnet$nodes) && nrow(myValues$fullnet$nodes) > 0
    cat("Is myValues$fullnet$nodes valid:", nodes_valid, "\n")
    selected_condition <- if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      "node_selection"
    } else if (!is.null(rv$pendingFocus)) {
      "pending_focus"
    } else if (!is.null(storage$focus)) {
      "current_focus"
    } else if (!is.null(input$plat) && input$plat != "") {
      "platform_selected"
    } else {
      "default"
    }
    cat("Selected condition:", selected_condition, "\n")
    if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      paste("Focus on", input$current_node_selection)
    } else if (!is.null(rv$pendingFocus)) {
      if (!nodes_valid) {
        cat("Error: myValues$fullnet$nodes is invalid, cannot map rv$pendingFocus:", rv$pendingFocus, "\n")
        return("Invalid network data")
      }
      selected_id <- myValues$fullnet$nodes$id[
        tolower(myValues$fullnet$nodes$TRAITID) == tolower(rv$pendingFocus)
      ]
      if (length(selected_id) == 0) {
        selected_id <- rv$selectedTrait %||% "Unknown"
        cat("No id found for rv$pendingFocus:", rv$pendingFocus, ", using fallback:", selected_id, "\n")
      }
      paste("Select Focus for", input$plat %||% "Unknown", ":", selected_id)
    } else if (!is.null(storage$focus)) {
      if (!nodes_valid) {
        cat("Error: myValues$fullnet$nodes is invalid, cannot map storage$focus:", storage$focus, "\n")
        return("Invalid network data")
      }
      selected_id <- myValues$fullnet$nodes$id[
        tolower(myValues$fullnet$nodes$TRAITID) == tolower(storage$focus)
      ]
      if (length(selected_id) == 0) {
        selected_id <- "Unknown"
        cat("No id found for storage$focus:", storage$focus, ", using fallback:", selected_id, "\n")
      }
      paste("Display Network for", input$plat %||% "Unknown", ":", selected_id)
    } else if (!is.null(input$plat) && input$plat != "") {
      paste("Select Trait for", input$plat)
    } else {
      "Select Platform and Trait"
    }
  })
  
  # Handle Update Network button
  observeEvent(input$updateNetworkBtn, {
    cat("updateNetworkBtn clicked, input$current_node_selection:", input$current_node_selection, 
        ", rv$pendingFocus:", rv$pendingFocus, ", input$trait:", input$trait, 
        ", rv$selectedTrait:", rv$selectedTrait, "\n")
    if (!is.null(rv$pendingFocus) && !is.null(input$trait) && input$trait == rv$selectedTrait) {
      storage$focus <- rv$pendingFocus
      session$sendCustomMessage(type = "clearNodeSelection", message = list())
      cat("Set storage$focus to rv$pendingFocus:", rv$pendingFocus, ", cleared input$current_node_selection\n")
    } else if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      storage$focus <- input$current_node_selection
      cat("Set storage$focus to input$current_node_selection:", input$current_node_selection, "\n")
    } else {
      cat("No valid focus, keeping storage$focus:", storage$focus, "\n")
    }
    networkTrigger(networkTrigger() + 1)
    cat("Triggered network update, storage$focus:", storage$focus, "\n")
  })
  

  observeEvent(input$confirmFocus, {
    removeModal()
    req(myValues$fullnet)
    if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0 && !is.null(rv$pendingFocus)) {
      selected_traitid <- rv$pendingFocus
      if (length(selected_traitid) == 0) {
        showNotification("Selected node not found in network.", type = "error")
        return()
      }
      storage$focus <- selected_traitid
      rv$useCaseActive <- FALSE
      showNotification(paste("Focusing network on node:", input$current_node_selection))
      cat("Focus node TRAITID:", selected_traitid, "\n")
    } else if (!is.null(rv$selectedTrait) && !is.null(input$trait)) {
      selected_traitid <- myValues$fullnet$nodes$TRAITID[
        tolower(myValues$fullnet$nodes$id) == tolower(input$trait)
      ]
      if (length(selected_traitid) == 0) {
        showNotification("Selected trait not found in network.", type = "error")
        return()
      }
      storage$focus <- selected_traitid
      rv$useCaseActive <- FALSE
      showNotification(paste("Displaying network for trait:", input$trait))
      cat("Display network for TRAITID:", selected_traitid, "\n")
    } else {
      showNotification("Please select a trait or node to focus on.", type = "error")
      return()
    }
    networkTrigger(networkTrigger() + 1)
  })
  

  observeEvent(input$saveUseCaseBtn, {
    if (!is.null(storage$subnet)) {
      showModal(modalDialog(
        title = "Save View",
        textInput(
          "useCaseName",
          "Enter a name for this view:",
          value = paste("View", length(rv$useCases) + 1)
        ),
        textAreaInput(
          "useCaseDescription",
          "Description (optional):",
          placeholder = "Add a short note about this saved view...",
          height = "80px",
          resize = "vertical"
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirmSaveUseCase", "Save")
        )
      ))
    } else {
      showModal(modalDialog(
        title = "Error",
        "No network is selected to save.",
        easyClose = TRUE
      ))
    }
  })
  

  observeEvent(input$confirmSaveUseCase, {
    req(input$useCaseName)
    
    rv$useCases[[input$useCaseName]] <- list(
      subnet      = storage$subnet,
      description = input$useCaseDescription
    )
    
    removeModal()
    showNotification(paste("Saved view:", input$useCaseName))
  })
  

  output$useCaseList <- renderUI({
    if (length(rv$useCases) > 0) {
      selectInput(
        "selectedUseCase",
        "Select a saved view:",
        choices = names(rv$useCases)
      )
    } else {
      h4("No saved views.")
    }
  })
  
  observeEvent(input$loadUseCaseBtn, {
    # 1) Prefer table selection, fallback to dropdown
    idx <- input$savedViewsTable_rows_selected
    if (length(idx) == 1 && length(rv$useCases) > 0) {
      view_names   <- names(rv$useCases)
      selectedName <- view_names[idx]
    } else {
      req(input$selectedUseCase)
      selectedName <- input$selectedUseCase
    }
    
    saved_obj <- rv$useCases[[selectedName]]
    if (is.null(saved_obj)) return(NULL)
    
    # old vs new format
    if (!is.null(saved_obj$nodes) && !is.null(saved_obj$edges)) {
      saved_net  <- saved_obj
      saved_desc <- NULL
    } else {
      saved_net  <- saved_obj$subnet
      saved_desc <- saved_obj$description
    }
    
    if (is.null(saved_net)) return(NULL)
    
    # restore subnet state
    storage$subnet <- saved_net
    
    if (nrow(saved_net$edges) == 0) {
      edges_df <- data.frame(
        from  = character(0),
        to    = character(0),
        title = character(0),
        color = character(0),
        width = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      edge_cols <- tolower(names(saved_net$edges))
      
      # guess source/target columns (same logic as in renderVisNetwork)
      from_col <- if ("from_id" %in% edge_cols) "from_id" else "from"
      to_col   <- if ("to_id"   %in% edge_cols) "to_id"   else "to"
      
      sign_vec <- if ("sign" %in% edge_cols) saved_net$edges$sign else rep(1, nrow(saved_net$edges))
      
      edges_df <- data.frame(
        from  = as.character(saved_net$edges[[from_col]]),
        to    = as.character(saved_net$edges[[to_col]]),
        title = paste(
          as.character(saved_net$edges[[from_col]]),
          "->",
          as.character(saved_net$edges[[to_col]])
        ),
        color = ifelse(sign_vec > 0, "blue", "red"),
        width = rep(1, nrow(saved_net$edges)),
        stringsAsFactors = FALSE
      )
      
      if (all(c("pvalue", "weight") %in% edge_cols)) {
        edges_df$title <- sprintf(
          "GGM_OLINK: %s->%s, p=%.2e, beta=%.3f",
          as.character(saved_net$edges[[from_col]]),
          as.character(saved_net$edges[[to_col]]),
          saved_net$edges$pvalue,
          saved_net$edges$weight
        )
      }
      
      edges_df <- dplyr::distinct(edges_df, from, to, .keep_all = TRUE)
    }
    
    nodes_raw <- saved_net$nodes
    n <- nrow(nodes_raw)
    
    if (!"id" %in% names(nodes_raw)) {
      stop("Saved view nodes table has no 'id' column.")
    }
    
    id_vec <- as.character(nodes_raw$id)
    
    if ("label" %in% names(nodes_raw)) {
      label_vec <- as.character(nodes_raw$label)
    } else {
      label_vec <- id_vec
    }
    
    if ("title" %in% names(nodes_raw)) {
      title_vec <- as.character(nodes_raw$title)
    } else {
      title_vec <- label_vec
    }
    
    color_vec <- if ("color" %in% names(nodes_raw)) {
      as.character(nodes_raw$color)
    } else {
      rep("#97c2fc", n)   
    }
    
    shape_vec <- if ("shape" %in% names(nodes_raw)) {
      as.character(nodes_raw$shape)
    } else {
      rep("dot", n)
    }
    
    size_vec <- if ("size" %in% names(nodes_raw)) {
      as.numeric(nodes_raw$size)
    } else {
      rep(25, n)
    }
    
    nodes_df <- data.frame(
      id    = id_vec,
      label = label_vec,
      title = title_vec,
      color = color_vec,
      shape = shape_vec,
      size  = size_vec,
      stringsAsFactors = FALSE
    )
    
    storage$edges <- edges_df
    storage$nodes <- nodes_df
    
    # you *can* keep focus info for other logic, but it won't drive the Saved View path
    if (nrow(saved_net$nodes) > 0 && "TRAITID" %in% names(saved_net$nodes)) {
      storage$focus <- saved_net$nodes$TRAITID[1]
    }
    
    rv$useCaseActive <- TRUE
    rv$pendingFocus  <- NULL
    
    rv$triggerNetworkAfterViewLoad <- TRUE
    
    msg <- paste("Loaded view:", selectedName)
    if (!is.null(saved_desc) && nzchar(saved_desc)) {
      msg <- paste0(msg, " — ", saved_desc)
    }
    showNotification(msg)
  })
  
  
  
  # Download current subnet (kept as-is; still called "use case" internally)
  output$downloadUseCase <- downloadHandler(
    filename = function() {
      paste("use_case_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      req(storage$subnet)
      saveRDS(storage$subnet, file)
    }
  )
  

  observeEvent(input$uploadUseCase, {
    req(input$uploadUseCase)
    use_case <- readRDS(input$uploadUseCase$datapath)
    
    # If user uploads plain subnet
    if (!is.null(use_case$nodes) && !is.null(use_case$edges)) {
      saved_net <- use_case
    } else if (!is.null(use_case$subnet)) {
      # If they upload an object wrapped like list(subnet = ..., description = ...)
      saved_net <- use_case$subnet
    } else {
      showNotification("Uploaded file is not a valid use-case/view object.", type = "error")
      return(NULL)
    }
    
    storage$subnet <- saved_net
    
    edges_df <- data.frame(
      from  = saved_net$edges$from,
      to    = saved_net$edges$to,
      title = saved_net$edges$title,
      color = saved_net$edges$color,
      width = saved_net$edges$width,
      stringsAsFactors = FALSE
    )
    
    nodes_df <- data.frame(
      id    = saved_net$nodes$id,
      label = saved_net$nodes$id,
      title = saved_net$nodes$id,
      color = saved_net$nodes$color,
      shape = saved_net$nodes$shape,
      stringsAsFactors = FALSE
    )
    nodes_df$size <- 25
    
    storage$edges <- edges_df
    storage$nodes <- nodes_df
    
    if (nrow(saved_net$nodes) > 0) {
      focus_traitid   <- saved_net$nodes$TRAITID[1]
      storage$focus   <- focus_traitid
      rv$pendingFocus <- NULL
      cat("Set storage$focus to:", focus_traitid, "\n")
      
      updateSelectInput(
        session,
        inputId  = "trait",
        selected = saved_net$nodes$id[1]
      )
    }
    
    rv$useCaseActive <- TRUE
    cat("Set rv$useCaseActive to TRUE, cleared rv$pendingFocus\n")
    
    networkTrigger(networkTrigger() + 1)
    updateTabItems(session, "tabs", "Tab1")
    showNotification("Saved view loaded successfully!")
  })
  
  # Table of saved views with descriptions
  output$savedViewsTable <- DT::renderDT({
    if (length(rv$useCases) == 0) return(NULL)
    
    df <- data.frame(
      View = names(rv$useCases),
      Description = sapply(rv$useCases, function(x) {
        # Backward compatible: if old format, no description
        if (is.list(x) && !is.null(x$description)) {
          x$description
        } else {
          ""
        }
      }),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      df,
      rownames = FALSE,
      selection = "single",
      options = list(
        pageLength = 5,
        lengthChange = FALSE,
        dom = "tip"   # table + info + pagination
      )
    )
  })
  
  # When user clicks a row in the table, sync the selectInput
  observeEvent(input$savedViewsTable_rows_selected, {
    idx <- input$savedViewsTable_rows_selected
    if (length(idx) == 1 && length(rv$useCases) > 0) {
      view_names <- names(rv$useCases)
      selected_name <- view_names[idx]
      updateSelectInput(
        session,
        "selectedUseCase",
        selected = selected_name
      )
    }
  })
  

  
  # Custom Home Page
  persistentHomePath <- "www/custom_home.html"
  observeEvent(input$homePageFile, {
    if (!is.null(input$homePageFile)) {
      ext <- tools::file_ext(input$homePageFile$name)
      if (ext == "html") {
        file.copy(input$homePageFile$datapath, persistentHomePath, overwrite = TRUE)
      } else if (ext == "md") {
        markdown_text <- readLines(input$homePageFile$datapath, warn = FALSE)
        html_content <- markdown::markdownToHTML(
          text = paste(markdown_text, collapse = "\n"),
          fragment.only = TRUE)
        writeLines(html_content, persistentHomePath)
      }
    }
  })
  
  output$customHomePage <- renderUI({
    invalidateLater(1000, session)
    if (file.exists(persistentHomePath)) {
      includeHTML(persistentHomePath)
    } else {
      includeHTML("www/default_home.html")
    }
  })
  

  stable_pick <- function(key, choices) {
    h <- digest(key, algo = "xxhash32")
    i <- strtoi(substr(h, 1, 7), base = 16L)
    choices[(i %% length(choices)) + 1L]
  }
  
  
  
  observe({
    req(myValues$fullnet)
    req(myValues$fullnet$nodes)
    
    list_items <- unique(myValues$fullnet$nodes$plat)
    
    # Same allowed sets as in the table
    shape_opts_all  <- c("circle", "diamond", "star", "square", "triangle")
    shape_defaults  <- setdiff(shape_opts_all, "circle")
    color_opts_all  <- c("Red", "Green", "Blue", "Purple",
                         "Grey", "Yellow", "Golden", "Pink", "Black")
    
    for (platform in list_items) {
      
      # SHAPE
      if (is.null(rv$customSelections[[platform]]$shape)) {
        rv$customSelections[[platform]]$shape <- stable_pick(platform, shape_defaults)
      }
      
      # COLOR
      if (is.null(rv$customSelections[[platform]]$color)) {
        rv$customSelections[[platform]]$color <- stable_pick(platform, color_opts_all)
      }
    }
  })
  
  
  output$dynamicTable <- renderDT({
    req(myValues$fullnet)
    
    list_items <- unique(myValues$fullnet$nodes$plat)
    myValues$listItems <- list_items
    
    shape_opts_all  <- c("circle", "diamond", "star", "square", "triangle")
    shape_defaults  <- setdiff(shape_opts_all, "circle")
    color_opts_all  <- c("Red", "Green", "Blue", "Purple",
                         "Grey", "Yellow", "Golden", "Pink", "Black")
    
    platform_label_map <- NULL
    if (!is.null(rv$platform_counts)) {
      platform_label_map <- setNames(
        paste0(rv$platform_counts$PLAT, " (", rv$platform_counts$trait_count, ")"),
        rv$platform_counts$PLAT
      )
    }
    
    data <- data.frame(
      #Platform = list_items,
      Platform = if (!is.null(platform_label_map)) {
        platform_label_map[list_items]
      } else {
        list_items
      },
      Shape = sapply(seq_along(list_items), function(i) {
        platform <- list_items[i]
        current <- if (!is.null(rv$customSelections[[platform]]$shape)) {
          rv$customSelections[[platform]]$shape
        } else {
          stable_pick(platform, shape_defaults)
        }
        options_html <- paste(
          sapply(shape_opts_all, function(opt) {
            if (opt == current) {
              sprintf('<option value="%s" selected>%s</option>', opt, opt)
            } else {
              sprintf('<option value="%s">%s</option>', opt, opt)
            }
          }),
          collapse = ""
        )
        sprintf('<select id="shape_%d">%s</select>', i, options_html)
      }),
      Color = sapply(seq_along(list_items), function(i) {
        platform <- list_items[i]
        current <- if (!is.null(rv$customSelections[[platform]]$color)) {
          rv$customSelections[[platform]]$color
        } else {
          stable_pick(platform, color_opts_all)
        }
        options_html <- paste(
          sapply(color_opts_all, function(opt) {
            if (opt == current) {
              sprintf('<option value="%s" selected>%s</option>', opt, opt)
            } else {
              sprintf('<option value="%s">%s</option>', opt, opt)
            }
          }),
          collapse = ""
        )
        sprintf('<select id="color_%d">%s</select>', i, options_html)
      }),
      stringsAsFactors = FALSE
    )
    
    datatable(
      data,
      class = "table table-striped table-hover",
      escape = FALSE,
      rownames = FALSE,
      options = list(
        dom = 't',
        paging = FALSE,
        info = FALSE,
        ordering = FALSE
      ),
      colnames = c("Platform (count)", "Shape", "Color")
    )
  })
  
  
  ###########
  observe({
    req(myValues$listItems)
    for(i in seq_along(myValues$listItems)) {
      platform <- myValues$listItems[i]
      shape_input <- input[[paste0("shape_", i)]]
      color_input <- input[[paste0("color_", i)]]
      if (!is.null(shape_input)) {
        if (is.null(rv$customSelections[[platform]])) {
          rv$customSelections[[platform]] <- list()
        }
        rv$customSelections[[platform]]$shape <- shape_input
      }
      if (!is.null(color_input)) {
        if (is.null(rv$customSelections[[platform]])) {
          rv$customSelections[[platform]] <- list()
        }
        rv$customSelections[[platform]]$color <- color_input
      }
    }
  })
  
  output$SHEETlist <- renderUI({
    if (!is.null(myValues$allSheets)) {
      tagList(
        #h4("Sheets that met the criteria"),
        checkboxGroupInput(
          inputId = "selectedSheetsInput",
          label = "Select sheets for further processing:",
          choices = myValues$allSheets,
          selected = myValues$sheets_ok
        )
      )
    }
  })
  
  storage <- reactiveValues(
    focus = NULL,
    subnet = NULL,
    edges = NULL,
    nodes = NULL
  )
  
  observeEvent(input$submitInfo, {
    req(myValues$datalist, myValues$fullnet)
    activeSheet <- input$sheetTabs
    if (!is.null(activeSheet)) {
      selected <- input[[paste0(activeSheet, "table_rows_selected")]]
      if (!is.null(selected) && length(selected) == 1) {
        sheetData <- myValues$datalist[[activeSheet]]
        focus_ids <- sheetData$from[selected]
        valid_ids <- focus_ids[tolower(focus_ids) %in% tolower(myValues$fullnet$nodes$TRAITID)]
        if (length(valid_ids) > 0) {
          storage$focus <- valid_ids
          selected_id <- myValues$fullnet$nodes$id[
            tolower(myValues$fullnet$nodes$TRAITID) == tolower(valid_ids)
          ]
          rv$selectedTrait <- selected_id
          rv$pendingFocus <- valid_ids
          session$sendCustomMessage(type = "clearNodeSelection", message = list())
          updateSelectInput(session, "trait", selected = selected_id)
          rv$useCaseActive <- FALSE
          updateTabItems(session, "tabs", "Tab1")
          showNotification(paste("Displaying network for table selection:", selected_id))
          cat("Table row selected, set storage$focus to:", valid_ids, ", rv$selectedTrait to:", selected_id, "\n")
          networkTrigger(networkTrigger() + 1)
        } else {
          showNotification("Selected trait not found in network.", type = "error")
          cat("No valid TRAITID found for table selection\n")
        }
      } else {
        showNotification("Please select exactly one row in the table.", type = "warning")
        cat("Invalid table row selection\n")
      }
    }
  })
  
  
  
  # Network rendering (degree neighbors + static layout for big graphs)
  output$network <- renderVisNetwork({
    networkTrigger()
    req(myValues$fullnet)
    

    focus_size     <- 25
    normal_size    <- 25
    focus_border   <- 3
    normal_border  <- 1
    incident_width <- 2
    normal_width   <- 1
    
    warn_threshold <- 400L   # only for warning, no hard cap
    
    ## --- helpers ----------------------------------------------------------------
    empty_edges <- function() {
      data.frame(
        from  = character(0),
        to    = character(0),
        title = character(0),
        color = character(0),
        width = numeric(0),
        stringsAsFactors = FALSE
      )
    }
    
    find_focus_id <- function(act_trait) {
      ndf <- myValues$fullnet$nodes
      hit <- ndf$id[tolower(ndf$TRAITID) == tolower(act_trait)]
      if (length(hit) == 1) return(hit)
      hit2 <- ndf$id[tolower(ndf$id) == tolower(act_trait)]
      if (length(hit2) == 1) return(hit2)
      character(0)
    }
    
    val_or_empty <- function(v, n_rows) {
      if (missing(v) || is.null(v)) return(character(n_rows))
      x <- as.character(v); x[is.na(x)] <- ""; trimws(x)
    }
    
    # render just the focus node (isolated)
    render_focus_only <- function(focus_id) {
      nodes_df <- myValues$fullnet$nodes %>%
        dplyr::mutate(id = as.character(id)) %>%
        dplyr::filter(id == focus_id)
      
      if (!"id" %in% names(nodes_df)) stop("Expected 'id' column in nodes")
      
      n <- nrow(nodes_df)
      label <- val_or_empty(nodes_df$label, n)
      if (all(label == "")) {
        label <- val_or_empty(nodes_df$SYMBOL, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$gene_symbol, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$TRAITID, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$name, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$Node, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$node, n)
      }
      label[label == ""] <- as.character(nodes_df$id)
      nodes_df$label <- label
      nodes_df$title <- nodes_df$label
      
      nodes_df$shape <- sapply(nodes_df$plat, function(p) rv$customSelections[[p]]$shape)
      nodes_df$color.background <- sapply(nodes_df$plat, function(p) rv$customSelections[[p]]$color)
      
      
      nodes_df$size        <- focus_size
      nodes_df$borderWidth <- focus_border
      
      edges_df <- empty_edges()
      
      storage$nodes  <- nodes_df
      storage$edges  <- edges_df
      storage$subnet <- list(nodes = data.frame(id = focus_id), edges = data.frame())
      
      visNetwork(nodes_df, edges_df, height = "800px", width = "100%") %>%
        visNodes(
          shadow  = list(enabled = TRUE, size = 10),
          scaling = list(min = 10, max = 30),
          font    = list(size = 16)
        ) %>%
        visLayout(randomSeed = 4711) %>%
        visPhysics(stabilization = FALSE) %>%
        visEdges(smooth = TRUE) %>%
        visOptions(
          #nodesIdSelection = list(enabled = TRUE),
          nodesIdSelection = FALSE,
          highlightNearest = FALSE
        ) %>%
        visEvents(
          selectNode   = "function(p){Shiny.onInputChange('current_node_selection', p.nodes);}",
          deselectNode = "function(p){Shiny.onInputChange('current_node_selection', []);}",
          click        = "function(p){Shiny.onInputChange('click', p.nodes[0]);}"
        ) %>%
        visInteraction(navigationButtons = FALSE)
    }
    
    ## --- inputs/state -----------------------------------------------------------
    # input$maxnodes is now "degree"
    deg_input <- suppressWarnings(as.integer(input$maxnodes))
    if (is.na(deg_input) || deg_input < 1) deg_input <- 1
    
    nodes_df <- NULL
    edges_df <- NULL
    
    ## --- use-case / Saved View branch -----------------------------------------
    if (rv$useCaseActive && !is.null(storage$subnet)) {
      
      nodes_df <- storage$nodes
      edges_df <- storage$edges
      
    } else {
      

      act_trait <- if (!is.null(storage$focus) && length(storage$focus) > 0) storage$focus else NULL
      #validate(need(!is.null(act_trait), "No active trait selected. Please select a trait and click 'Focus on'."))
      
      # If no active trait, clear the network
      if (is.null(act_trait) || !nzchar(act_trait)) {
        showNotification(
          "No active trait selected. Please select a trait and click 'Focus on'.",
          type = "message",
          duration = 5
        )
        
        # Make sure subnet state is cleared too
        storage$subnet <- NULL
        storage$nodes  <- NULL
        storage$edges  <- NULL
        
        return(NULL)
      }
      
      focus_id <- find_focus_id(act_trait)
      validate(need(length(focus_id) == 1, paste("Focus node not found for", act_trait)))
      focus_id <- as.character(focus_id)
      

      nodes_all <- myValues$fullnet$nodes %>%
        dplyr::mutate(id = as.character(id))
      
      edges_all <- myValues$fullnet$edges
      edge_cols <- tolower(names(edges_all))
      from_col  <- if ("from_id" %in% edge_cols) "from_id" else stop("No 'from_id' in edges")
      to_col    <- if ("to_id"   %in% edge_cols) "to_id"   else stop("No 'to_id' in edges")
      
      edges_all <- edges_all %>%
        dplyr::mutate(
          from_chr = as.character(.data[[from_col]]),
          to_chr   = as.character(.data[[to_col]])
        )
      
      # --- igraph distances ------------------------------------------------------
      res <- tryCatch({
        g <- igraph::graph_from_data_frame(
          d = data.frame(
            from = edges_all$from_chr,
            to   = edges_all$to_chr,
            stringsAsFactors = FALSE
          ),
          directed = FALSE,
          vertices = nodes_all %>%
            dplyr::rename(name = id)
        )
        
        if (!focus_id %in% igraph::V(g)$name) {
          stop(paste("Focus node", focus_id, "not found in graph vertices"))
        }
        
        dists <- igraph::distances(g, v = focus_id)[1, ]
        dists <- dists[is.finite(dists)]
        
        list(dists = dists)
      }, error = function(e) {
        cat("igraph error for", act_trait, ":", e$message, "\n")
        showNotification(
          paste("Error computing graph distances for", act_trait, ":", e$message),
          type = "error"
        )
        return(NULL)
      })
      
      if (is.null(res)) return(NULL)
      
      dists <- res$dists
      
      # If focus has no neighbors
      if (!length(dists) || all(dists == 0)) {
        return(render_focus_only(focus_id))
      }
      
      max_deg <- max(dists)
      
      # clamp selected degree
      if (deg_input > max_deg) deg_input <- max_deg
      
      # cumulative node counts within degree <= d
      total_by_deg <- sapply(1:max_deg, function(d) sum(dists <= d))
      
      # Update dropdown labels: "d (nNodes)"
      choices_vec <- as.character(1:max_deg)
      names(choices_vec) <- paste0(1:max_deg, " (", total_by_deg, ")")
      
      updateSelectInput(
        session,
        "maxnodes",
        label   = "Neighbors (degree)",
        choices = choices_vec,
        selected = as.character(deg_input)
      )
      
      chosen_deg     <- deg_input
      chosen_nodes_n <- total_by_deg[chosen_deg]
      
      warn_threshold <- 400L
      if (!is.na(chosen_nodes_n) && chosen_nodes_n > warn_threshold) {
        showNotification(
          paste0(
            "Selected degree ", chosen_deg, " includes ~", chosen_nodes_n,
            " nodes. Large networks may slow down or freeze your browser."
          ),
          type = "warning",
          duration = 8
        )
      }
      
      # nodes within chosen degree (includes focus at distance 0)
      keep_nodes <- names(dists)[dists <= chosen_deg]
      
      subnet_nodes <- nodes_all %>%
        dplyr::filter(id %in% keep_nodes)
      
      subnet_edges <- edges_all %>%
        dplyr::filter(from_chr %in% keep_nodes & to_chr %in% keep_nodes)
      
      subnet <- list(nodes = subnet_nodes, edges = subnet_edges)
      storage$subnet <- subnet
      
      ## --- edges (0-row safe) ---------------------------------------------------
      edge_cols_sub <- tolower(names(subnet$edges))
      from_col_sub  <- if ("from_id" %in% edge_cols_sub) "from_id" else from_col
      to_col_sub    <- if ("to_id"   %in% edge_cols_sub) "to_id"   else to_col
      
      if (nrow(subnet$edges) == 0) {
        edges_df <- empty_edges()
      } else {
        sign_vec <- if ("sign" %in% edge_cols_sub) subnet$edges$sign else rep(1, nrow(subnet$edges))
        edges_df <- data.frame(
          from  = as.character(subnet$edges[[from_col_sub]]),
          to    = as.character(subnet$edges[[to_col_sub]]),
          title = paste(
            as.character(subnet$edges[[from_col_sub]]),
            "->",
            as.character(subnet$edges[[to_col_sub]])
          ),
          color = ifelse(sign_vec > 0, "blue", "red"),
          width = rep(normal_width, nrow(subnet$edges)),
          stringsAsFactors = FALSE
        )
        if (all(c("pvalue", "weight") %in% edge_cols_sub)) {
          edges_df$title <- sprintf(
            "GGM_OLINK: %s->%s, p=%.2e, beta=%.3f",
            as.character(subnet$edges[[from_col_sub]]),
            as.character(subnet$edges[[to_col_sub]]),
            subnet$edges$pvalue,
            subnet$edges$weight
          )
        }
        edges_df <- dplyr::distinct(edges_df, from, to, .keep_all = TRUE)
      }
      
      ## --- nodes (preserve platform styling, force labels) ----------------------
      nodes_df <- subnet$nodes
      if (!"id" %in% names(nodes_df)) stop("Expected 'id' column in nodes")
      nodes_df$id <- as.character(nodes_df$id)
      
      n <- nrow(nodes_df)
      label <- val_or_empty(nodes_df$label, n)
      if (all(label == "")) {
        label <- val_or_empty(nodes_df$SYMBOL, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$gene_symbol, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$TRAITID, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$name, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$Node, n)
        if (all(label == "")) label <- val_or_empty(nodes_df$node, n)
      }
      label[label == ""] <- as.character(nodes_df$id)
      nodes_df$label <- label
      nodes_df$title <- nodes_df$label
      
      nodes_df$size        <- ifelse(tolower(nodes_df$id) == tolower(focus_id), focus_size, normal_size)
      nodes_df$borderWidth <- ifelse(tolower(nodes_df$id) == tolower(focus_id), focus_border, normal_border)
      
      if (nrow(edges_df) > 0 && length(focus_id) == 1) {
        edges_df$width[
          tolower(edges_df$from) == tolower(focus_id) |
            tolower(edges_df$to)   == tolower(focus_id)
        ] <- incident_width
      }
      
      valid_nodes <- nodes_df$id
      edges_df <- edges_df %>%
        dplyr::filter(from %in% valid_nodes & to %in% valid_nodes)
      
      storage$nodes <- nodes_df
      storage$edges <- edges_df
    }
    
    if (is.null(nodes_df) || nrow(nodes_df) == 0) {
      showNotification("Cannot render network: no valid nodes.", type = "error")
      return(NULL)
    }
    
    n_nodes <- nrow(nodes_df)
    
    # For big graphs, precompute static layout with igraph
    if (n_nodes > 150 && nrow(edges_df) > 0) {
      g_sub <- igraph::graph_from_data_frame(
        d = data.frame(
          from = as.character(edges_df$from),
          to   = as.character(edges_df$to),
          stringsAsFactors = FALSE
        ),
        directed = FALSE,
        vertices = nodes_df %>% dplyr::mutate(id = as.character(id)) %>%
          dplyr::rename(name = id)
      )
      
      lay <- igraph::layout_with_fr(g_sub)
      nodes_df$x <- lay[, 1] * 200
      nodes_df$y <- lay[, 2] * 200
    }
    
    # base network (common)
    net <- visNetwork(nodes_df, edges_df, height = "800px", width = "100%") %>%
      visNodes(
        shadow  = list(enabled = TRUE, size = 10),
        scaling = list(min = 10, max = 30),
        font    = list(size = 16)
      ) %>%
      visEdges(smooth = TRUE) %>%
      visOptions(
        #nodesIdSelection = list(enabled = TRUE),
        nodesIdSelection = FALSE,
        highlightNearest = FALSE
      ) %>%
      visInteraction(navigationButtons = FALSE)
    
    if (n_nodes > 150) {
      # LARGE graph: use static layout, no physics
      net <- net %>%
        visPhysics(enabled = FALSE) %>%
        visEvents(
          selectNode   = "function(p){Shiny.onInputChange('current_node_selection', p.nodes);}",
          deselectNode = "function(p){Shiny.onInputChange('current_node_selection', []);}",
          click        = "function(p){Shiny.onInputChange('click', p.nodes[0]);}"
        )
    } else {
      # SMALL/MEDIUM graph: let vis.js settle, then freeze physics
      net <- net %>%
        visLayout(randomSeed = 4711) %>%
        visPhysics(stabilization = list(enabled = TRUE, iterations = 300)) %>%
        visEvents(
          stabilized   = "function () { this.setOptions({ physics: { enabled: false } }); }",
          selectNode   = "function(p){Shiny.onInputChange('current_node_selection', p.nodes);}",
          deselectNode = "function(p){Shiny.onInputChange('current_node_selection', []);}",
          click        = "function(p){Shiny.onInputChange('click', p.nodes[0]);}"
        )
    }
    
    net
  })
  

  
  # Network table
  output$tbl.networktable <- DT::renderDataTable({
    req(storage$edges, storage$nodes, storage$focus)
    selected <- input$current_node_selection
    cat("Rendering tbl.networktable, selected node:", selected, ", storage$focus:", storage$focus, "\n")
    print(head(storage$edges, n = min(nrow(storage$edges), 5)))
    
    if (is.null(selected) || length(selected) == 0) {
      cat("No node selected, returning empty table\n")
      showNotification("No node selected.", type = "message")
      return(data.frame(Message = "No node selected"))
    }
    
    # Validate myValues$fullnet$nodes
    if (!is.data.frame(myValues$fullnet$nodes) || nrow(myValues$fullnet$nodes) == 0) {
      cat("Error: myValues$fullnet$nodes is not a valid data frame or is empty\n")
      showNotification("Network data is invalid or empty.", type = "error")
      return(data.frame(Message = "Invalid network data"))
    }
    
    selected_traitid <- myValues$fullnet$nodes$TRAITID[tolower(myValues$fullnet$nodes$id) == tolower(selected)]
    if (length(selected_traitid) == 0) {
      cat("No TRAITID found for selected node:", selected, "\n")
      showNotification("Selected node not found in network.", type = "error")
      return(data.frame(Message = "Node not found"))
    }
    
    central_traitid <- if (length(storage$focus) > 1) {
      cat("Multiple central TRAITIDs:", storage$focus, "using first\n")
      storage$focus[1]
    } else {
      storage$focus
    }
    
    node_mapping <- myValues$fullnet$nodes %>% dplyr::select(TRAITID, id)
    relevant_mapping <- node_mapping %>% dplyr::filter(tolower(TRAITID) %in% tolower(c(selected_traitid, central_traitid)))
    cat("Relevant node mapping:\n")
    print(relevant_mapping)
    
    invalid_from <- setdiff(storage$edges$from, storage$nodes$id)
    invalid_to <- setdiff(storage$edges$to, storage$nodes$id)
    cat("Invalid from nodes:", invalid_from, ", Invalid to nodes:", invalid_to, "\n")
    
    is_central <- tolower(selected_traitid) == tolower(central_traitid)
    selected_ids <- node_mapping$id[tolower(node_mapping$TRAITID) == tolower(selected_traitid)]
    central_ids <- node_mapping$id[tolower(node_mapping$TRAITID) == tolower(central_traitid)]
    
    central_fullnet_edges <- myValues$fullnet$edges %>% 
      dplyr::filter(tolower(from_id) == tolower(central_traitid) | tolower(to_id) == tolower(central_traitid))
    cat("Central fullnet edges:\n")
    print(head(central_fullnet_edges, n = min(nrow(central_fullnet_edges), 5)))
    
    filtered_edges <- storage$edges %>% 
      dplyr::filter(
        tolower(from) %in% tolower(selected_ids) |
          tolower(to) %in% tolower(selected_ids)
      ) %>% 
      dplyr::distinct(from, to, .keep_all = TRUE)
    
    if (nrow(filtered_edges) > 0) {
      cat("Filtered edges for selected node:\n")
      print(head(filtered_edges, n = min(nrow(filtered_edges), 5)))
    } else {
      cat("No edges found for selected node:", selected, "\n")
      filtered_edges <- data.frame(Message = "No edges found for selected node")
    }
    
    filtered_edges
  }, options = list(lengthChange = FALSE))
  
  
  
  output$EXPORTtable <- DT::renderDataTable({
    req(storage$subnet)
    zwi <- storage$subnet$edges
    out <- data.frame(
      ID = zwi$id,
      TYPE = zwi$type,
      TRAITID1 = zwi$from_id,
      TRAIT1 = zwi$from,
      BETA = zwi$sign * zwi$weight,
      PVALUE = zwi$pvalue,
      TRAITID2 = zwi$to_id,
      TRAIT2 = zwi$to,
      stringsAsFactors = FALSE
    )
    out
  }, filter = 'top', rownames = FALSE, selection = 'single', server = FALSE,
  extensions = 'Buttons',
  options = list(pageLength = 10, scrollX = TRUE, dom = 'Blfrtip',
                 buttons = list("copy", list(extend = "collection",
                                             buttons = c("csv", "excel", "pdf"),
                                             text = "Download"))))
  
  output$dynamicTitle <- renderText({
    if (is.null(input$titleInput)) "Multi-omics Server" else input$titleInput
  })
  
  output$confirm_btn <- renderUI({
    if (!is.null(myValues$allSheets)) {
      actionButton("confirmBtn", "Process", class = "submit-data-btn")
    }
  })
  
  processedCustomData <- eventReactive(input$confirmBtn, {
    req(input$selectedSheetsInput)
    list_items <- myValues$listItems
    if (is.null(list_items)) {
      list_items <- unique(myValues$fullnet$nodes$plat)
    }
    rows <- lapply(seq_along(list_items), function(i) {
      platform <- list_items[i]
      shape_val <- if (!is.null(rv$customSelections[[platform]]$shape)) 
        rv$customSelections[[platform]]$shape else "circle"
      color_val <- if (!is.null(rv$customSelections[[platform]]$color)) 
        rv$customSelections[[platform]]$color else "Red"
      data.frame(
        Platform = platform,
        Shape = shape_val,
        Color = color_val,
        stringsAsFactors = FALSE
      )
    })
    result <- do.call(rbind, rows)
    attr(result, "Sheets_selections") <- input$selectedSheetsInput
    result
  })
  
  observeEvent(input$confirmBtn, {
    withProgress(message = "Generating network...", value = 0, {
      incProgress(0.3, detail = "Processing network appearance settings...")
      Sys.sleep(1)
      
      customSettings <- processedCustomData()
      if (!is.null(customSettings)) {
        isolate({
          for (i in 1:nrow(customSettings)) {
            platform <- customSettings$Platform[i]
            newColor <- customSettings$Color[i]
            newShape <- customSettings$Shape[i]
            myValues$fullnet$nodes$color[myValues$fullnet$nodes$plat == platform] <- newColor
            myValues$fullnet$nodes$shape[myValues$fullnet$nodes$plat == platform] <- newShape
          }
        })
      }
      
      Sys.sleep(1)
      incProgress(0.6, detail = "Updating selected sheets...")
      
      myValues$sheets_ok <- input$selectedSheetsInput
      
      saveRDS(list(
        datalist   = myValues$datalist,
        sheets_ok  = myValues$sheets_ok,
        sheets_flag = myValues$sheets_flag,
        anno       = myValues$anno,
        custom     = myValues$custom,
        fullnet    = myValues$fullnet
      ), data_file)
      
      saveRDS(list(customSelections = rv$customSelections), custom_settings_file)
      cat("Saved data to", data_file, "and settings to", custom_settings_file, "\n")
      
      incProgress(1, detail = "Switching to Tables tab...")
      Sys.sleep(1)
    })
    rv$dataLoaded <- TRUE
    
    updateTabItems(session, "tabs", "Tab2")
    processClicked(TRUE)
  })

  
  output$tabs <- renderUI({
    req(myValues$datalist)
    sheets_ok <- myValues$sheets_ok
    if (is.null(sheets_ok) || length(sheets_ok) == 0) {
      return(h4("No sheets selected for processing."))
    }
    nTabs <- length(sheets_ok)
    myTabs <- lapply(seq_len(nTabs), function(i) {
      sheetName <- sheets_ok[i]
      tabPanel(title = sheetName,
               DT::dataTableOutput(paste0(sheetName, "table"))
      )
    })
    tabBox(id = "sheetTabs", width = 12, !!!myTabs)
  })
  
  observe({
    req(myValues$datalist)
    sheets_ok <- myValues$sheets_ok
    for(sheet in sheets_ok) {
      local({
        sh <- sheet
        output[[paste0(sh, "table")]] <- DT::renderDataTable({
          myValues$datalist[[sh]]
        }, options = list(pageLength = 10, scrollX = TRUE), 
        rownames = FALSE, selection = 'single')
      })
    }
  })
  
  output$SAMPLEtable <- renderUI({
    req(myValues$sheets_flag)
    tagList(
      h4("Sheets That Did Not Qualify"),
      DT::dataTableOutput("failedTable")
    )
  })
  
  output$failedTable <- DT::renderDataTable({
    data.frame(Sheet_Flag = myValues$sheets_flag, stringsAsFactors = FALSE)
  }, options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  
  output$processedTable <- DT::renderDataTable({
    req(processedCustomData())
    processedCustomData()
  })
  
  output$dynamic_legend <- renderUI({
    req(myValues$fullnet)
    legend_data <- unique(myValues$fullnet$nodes[, c("plat", "color", "shape")])
    tagList(
      tags$h4("Legend"),
      lapply(1:nrow(legend_data), function(i) {
        platform <- legend_data$plat[i]
        color <- legend_data$color[i]
        shape <- tolower(legend_data$shape[i])
        shapeDiv <- switch(shape,
                           "diamond" = tags$div(
                             style = sprintf("width: 20px; height: 20px; background-color: %s; transform: rotate(45deg); margin-right: 5px;", color)
                           ),
                           "circle" = tags$div(
                             style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 50%%; margin-right: 5px;", color)
                           ),
                           "star" = tags$span(
                             "★",
                             style = sprintf("color: %s; font-size: 20px; margin-right: 5px;", color)
                           ),
                           "triangle" = tags$div(
                             style = sprintf("width: 0; height: 0; border-left: 10px solid transparent; border-right: 10px solid transparent; border-bottom: 20px solid %s; margin-right: 5px;", color)
                           ),
                           tags$div(
                             style = sprintf("width: 20px; height: 20px; background-color: %s; margin-right: 5px;", color)
                           )
        )
        tags$div(
          style = "display: flex; align-items: center; margin-bottom: 5px;",
          shapeDiv,
          tags$span(platform, style = "font-size: 14px;")
        )
      })
    )
  })
  
  observe({
    session$sendCustomMessage("change_skin", input$skin_color)
  })
  
  output$networkTableTitle <- renderUI({
    cat("Rendering networkTableTitle, storage$focus:", storage$focus, ", input$current_node_selection:", input$current_node_selection, "\n")
    if (!is.null(storage$focus)) {
      selected_id <- myValues$fullnet$nodes$id[
        tolower(myValues$fullnet$nodes$TRAITID) == tolower(storage$focus)
      ]
      h4(paste("Associations for", selected_id))
    } else {
      h4("")
    }
  })
  
  
  # Network instructions
  output$networkInstructions <- renderUI({
    req(storage$subnet$nodes)
    cat("Rendering network instructions, subnet nodes:", nrow(storage$subnet$nodes), "\n")
    div(
      style = "margin-top: 10px; padding: 10px; font-size: 14px; color: #333; background-color: #f8f8f8; border-radius: 5px;",
      "To explore the network: (1) Select a trait from the dropdown or click a node to view its associations in the table below. (2) Click 'Focus on' to update the network to show the selected trait or node's neighbors. (3) Drag nodes to reposition them for better visualization without changing the network."
    )
  })
  
  
  observeEvent(input$resetNetwork, {
    cat("Reset Network clicked\n")
    storage$focus <- NULL
    storage$subnet <- NULL
    storage$nodes <- NULL
    storage$edges <- NULL
    rv$selectedTrait <- NULL
    rv$pendingFocus <- NULL
    session$sendCustomMessage(type = "clearNodeSelection", message = list())
    updateSelectInput(session, "trait", selected = character(0))
    updateSelectInput(session, "plat", selected = "Unknown")
    rv$useCaseActive <- FALSE
    networkTrigger(networkTrigger() + 1)
    showNotification("Network reset.", type = "message")
    cat("Cleared storage$focus, rv$selectedTrait, rv$pendingFocus, input$current_node_selection, input$trait, input$plat\n")
  })
  
  observeEvent(input$current_node_selection, {
    cat("input$current_node_selection changed to:", input$current_node_selection, "\n")
    if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      # Validate selected node
      if (!input$current_node_selection %in% myValues$fullnet$nodes$id) {
        cat("Error: Selected node", input$current_node_selection, "not found in myValues$fullnet$nodes$id\n")
        showNotification("Selected node not found in network.", type = "error")
        session$sendCustomMessage(type = "clearNodeSelection", message = list())
        return()
      }
      cat("Node selected, updating UI but not network, selected node:", input$current_node_selection, "\n")
      rv$pendingFocus <- NULL # Clear pending focus to prioritize node selection
      rv$selectedTrait <- NULL # Clear selected trait for consistency
    } else {
      cat("No node selected, clearing selection\n")
    }
  })
  
  # ─────────────────────────────────────────────────────────────────────
  # Download handler
  # ─────────────────────────────────────────────────────────────────────
  output$downloadExcel <- downloadHandler(
    filename = function() {
      paste0("all_data_with_customizations_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(myValues$sheets_ok, myValues$datalist, myValues$anno)
      
      # 1) Build a named list of data.frames for each selected sheet
      sheets_to_write <- list()
      for (sh in myValues$sheets_ok) {
        sheets_to_write[[sh]] <- myValues$datalist[[sh]]
      }
      
      # 2) Add the annotation table
      sheets_to_write[["Annotation"]] <- myValues$anno
      
      # 3) Add “Failed Sheets” if there were any
      if (!is.null(myValues$sheets_flag) && length(myValues$sheets_flag) > 0) {
        sheets_to_write[["Failed_Sheets"]] <-
          data.frame(Sheet_Flag = myValues$sheets_flag, stringsAsFactors = FALSE)
      }
      
      # 4) Build a “Customizations” sheet from processedCustomData()
      cust_df <- processedCustomData()
      if (!is.null(cust_df)) {
        sheets_sel <- attr(cust_df, "Sheets_selections")
        cust_df2 <- cust_df
        cust_df2$Selected_Sheets <- paste(sheets_sel, collapse = ", ")
        sheets_to_write[["Customizations"]] <- cust_df2
      }
      
      # 5) write them all into a single .xlsx
      writexl::write_xlsx(sheets_to_write, path = file)
    }
  )
}
# ====================================================
# UI Definition
# ====================================================

js <- "Shiny.addCustomMessageHandler('change_skin', function(skin) {
  document.body.className = skin;
});"

max_nodes_list <- c(1, 2, 4, 6, 8, 10)

ui <- dashboardPage(
  skin = "red",
  
  dashboardHeader(
    title = textOutput("dynamicTitle")
  ),
  
  dashboardSidebar(
    width = 150,
    sidebarMenu(
      id       = "tabs",
      selected = "Tab6",
      menuItem("Home",         tabName = "Tab6", icon = icon("th")),
      menuItem("Data Formats", tabName = "Tab4", icon = icon("th")),
      menuItem("Import Data",  tabName = "Tab7", icon = icon("th")),
      menuItem("Tables",       tabName = "Tab2", icon = icon("th")),
      menuItem("Network",      tabName = "Tab1", icon = icon("dashboard")),
      menuItem("Export",       tabName = "Tab3", icon = icon("th")),
      menuItem("Saved-views",  tabName = "Tab5", icon = icon("th"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    
    # -------------------------------------------------
    # JS handlers & custom styles
    # -------------------------------------------------
    tags$head(
      tags$script(HTML("
        Shiny.addCustomMessageHandler('change_skin', function(skin) {
          var body = document.body;
          body.className = body.className.replace(/\\bskin-\\S+/g, '');
          body.className += ' skin-' + skin;
        });

        function clearSelectionHandler(type, delay) {
          return function(message) {
            setTimeout(function() {
              var network = document.getElementById('network').chart;
              if (network) {
                console.log('Executing ' + type + ', current nodes:', network.getSelectedNodes());
                network.unselectAll();
                Shiny.onInputChange('current_node_selection', []);
                console.log(type + ' completed:', network.getSelectedNodes());
              } else {
                console.log(type + ' failed: network not found');
                Shiny.onInputChange('current_node_selection', []);
              }
            }, delay || 0);
          }
        }

        Shiny.addCustomMessageHandler('clearNodeSelection', clearSelectionHandler('clearNodeSelection', 0));
        Shiny.addCustomMessageHandler('forceClearNodeSelection', clearSelectionHandler('forceClearNodeSelection', 0));
        Shiny.addCustomMessageHandler('forceClearNodeSelectionDelayed', function(message) {
          clearSelectionHandler('forceClearNodeSelectionDelayed', message.delay)(message);
        });
        Shiny.addCustomMessageHandler('forceServerClearNodeSelection', function(message) {
          Shiny.onInputChange('current_node_selection', []);
        });
      ")),
      
      tags$style(HTML("
        .custom-usecase-btn {
          background-color: #FF5733;
          color: white;
          border: none;
          border-radius: 5px;
          padding: 10px 20px;
          font-size: 16px;
        }
        .custom-usecase-btn:hover { background-color: #E04E2F; }

        .focus-network-btn {
          background-color: #33A1FF;
          color: white;
          border: none;
          border-radius: 5px;
          padding: 10px 20px;
          font-size: 16px;
        }
        .focus-network-btn:hover { background-color: #2A8BD6; }

        .download-usecase-btn {
          background-color: #33A1FF;
          color: white;
          border: none;
          border-radius: 5px;
          padding: 10px 20px;
          font-size: 16px;
        }
        .download-usecase-btn:hover { background-color: #2A8BD6; }

        .clear-data-btn {
          background-color: #0000FF;
          color: white;
          border: none;
          border-radius: 5px;
          padding: 10px 20px;
          margin-top: 10px;
          font-size: 16px;
        }
        .clear-data-btn:hover { background-color: #E02F2F; }

        .reset-network-btn {
          background-color: #0000FF;
          color: white;
          border: none;
          border-radius: 5px;
          padding: 10px 20px;
          font-size: 16px;
        }
        .reset-network-btn:hover { background-color: #E02F2F; }

        .metric-card {
          background: linear-gradient(135deg, #eaf3ff 0%, #f8fbff 100%);
          border-left: 4px solid #1565C0;
          padding: 12px 16px;
          border-radius: 8px;
          margin: 8px 0;
          box-shadow: 0 1px 3px rgba(0,0,0,0.08);
        }
        .metric-card .shiny-text-output {
          font-weight: 700;
          font-size: 18px;
          color: #0d47a1;
        }

        .btn-toolbar-top .btn {
          font-weight: 600;
          border-radius: 6px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
          transition: all 0.2s ease-in-out;
        }
        .btn-toolbar-top .btn:hover { transform: translateY(-2px); }
      "))
    ),
    
    # -------------------------------------------------
    # TABS
    # -------------------------------------------------
    tabItems(
      
      # ----------------------------
      # Tab 1 — NETWORK
      # ----------------------------
      tabItem(
        tabName = "Tab1",
        fluidPage(
          
          fluidRow(
            column(
              width = 3,
              selectInput("plat", "Select Platform", choices = NULL)
            ),
            column(
              width = 3,
              uiOutput("SelectTrait"),
              actionButton(
                "updateNetworkBtn",
                label = textOutput("updateNetworkLabel", inline = TRUE),
                class = "focus-network-btn"
              )
            ),
            column(
              width = 2,
              selectInput(
                "maxnodes",
                "Degree (neighbors)",
                choices  = max_nodes_list,
                selected = max_nodes_list[2]
              )
            )
          ),
          
          fluidRow(
            column(width = 9, visNetworkOutput("network")),
            column(width = 3, uiOutput("dynamic_legend"))
          ),
          
          fluidRow(
            column(
              width = 10,
              conditionalPanel(
                condition = "input.current_node_selection && input.current_node_selection.length > 0",
                uiOutput("networkTableTitle"),
                DT::dataTableOutput("tbl.networktable")
              )
            )
          ),
          
          fluidRow(
            column(
              width = 12,
              actionButton("saveUseCaseBtn", "Save as View", class = "custom-usecase-btn"),
              downloadButton("downloadUseCase", "Download View", class = "download-usecase-btn")
            )
          ),
          
          fluidRow(
            column(
              width = 12,
              box(
                width       = NULL,
                title       = "Network Instructions",
                status      = "info",
                solidHeader = TRUE,
                uiOutput("networkInstructions"),
                style       = "font-size: 14px; color: #333;"
              )
            )
          )
        )
      ),
      
      # ----------------------------
      # Tab 2 — TABLES
      # ----------------------------
      tabItem(
        tabName = "Tab2",
        fluidPage(
          fluidRow(
            h3("Select a row from the table and click Focus network"),
            actionButton("submitInfo", "Focus network on selection", class = "focus-network-btn"),
            uiOutput("tabs")
          )
        )
      ),
      
      # ----------------------------
      # Tab 3 — EXPORT
      # ----------------------------
      tabItem(
        tabName = "Tab3",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              h3("Export"),
              p("Download or view the network data of the active view"),
              DT::dataTableOutput("EXPORTtable"),
              hr(),
              downloadButton("downloadExcel", "Download All Data as Excel", class = "btn btn-primary")
            )
          )
        )
      ),
      
      # ----------------------------
      # Tab 4 — DATA FORMATS
      # ----------------------------
      tabItem(
        tabName = "Tab4",
        fluidPage(
          fluidRow(
            column(width = 12, h3("HowTo"), includeHTML("howto.html"))
          )
        )
      ),
      
      # ----------------------------
      # Tab 5 — SAVED VIEWS
      # ----------------------------
      tabItem(
        tabName = "Tab5",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              h3("Saved Views"),
              uiOutput("useCaseList"),
              actionButton("loadUseCaseBtn", "Load Selected View", class = "custom-usecase-btn"),
              fileInput("uploadUseCase", "Load Saved View (RDS)", accept = ".rds"),
              br(),
              h4("Saved Views Overview"),
              DT::DTOutput("savedViewsTable")
            )
          )
        )
      ),
      
      # ----------------------------
      # Tab 6 — HOME PAGE
      # ----------------------------
      tabItem(
        tabName = "Tab6",
        fluidPage(
          uiOutput("customHomePage")
        )
      ),
      
      # ----------------------------
      # Tab 7 — IMPORT DATA
      # ----------------------------
      tabItem(
        tabName = "Tab7",
        fluidPage(
          uiOutput("importControls"),
          
          fluidRow(
            column(width = 5, uiOutput("SAMPLEtable"))
          ),
          
          fluidRow(
            column(width = 4)
          ),
          
          fluidRow(
            column(width = 12)
          )
        )
      )
    )
  )
)

# ====================================================
# Run the Application
# ====================================================

shinyApp(ui = ui, server = server)