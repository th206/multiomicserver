library(shiny)
library(shinydashboard)
library(readxl)
library(dplyr)
library(tidyr)
library(visNetwork)
library(DT)
library(shinyWidgets)
library(igraph)

# ====================================================
# Helper Functions for Network Operations
# ====================================================

neighbors <- function(nodes, network) {
  nodes <- tolower(nodes)
  ix <- unlist(lapply(nodes, function(x) {
    union(which(tolower(network$edges$from_id) == x), which(tolower(network$edges$to_id) == x))
  })) %>% unique()
  unique(c(network$edges$from_id[ix], network$edges$to_id[ix]))
}

maxneighbors <- function(nodes, network, limit = 0) {
  if (limit == 0) limit <- 1E99
  nnodes <- length(nodes)
  nnodeslast <- 0
  nodeslast <- nodes
  while ((nnodes > nnodeslast) & (nnodes <= limit)) {
    nodeslast <- nodes
    nodes <- neighbors(nodes, network)
    nnodeslast <- nnodes
    nnodes <- length(nodes)
  }
  nodeslast
}

maxneighbors_noSTAT <- function(nodes, network, limit = 0) {
  if (limit == 0) limit <- 1E99
  nodesin <- nodes
  nnodes <- length(nodes)
  nnodeslast <- 0
  nodeslast <- nodes
  while ((nnodes > nnodeslast) & (nnodes <= limit)) {
    nodeslast <- nodes
    nnodeslast <- nnodes
    nogrow_nodes <- nodes[grep("STAT: |GWAS: ", nodes)]
    grow_nodes <- nodes[grep("STAT: |GWAS: ", nodes, invert = TRUE)]
    nodes <- c(nogrow_nodes, neighbors(grow_nodes, network))
    nnodes <- length(nodes)
  }
  unique(c(nodesin, nodeslast))
}

nodes2network <- function(nodes, network) {
  nodes <- tolower(nodes)
  cat("Input nodes to nodes2network:", paste(head(nodes, 5), collapse = ", "), "\n")
  ix_from <- unique(unlist(lapply(nodes, function(x) which(tolower(network$edges$from_id) == x))))
  ix_to <- unique(unlist(lapply(nodes, function(x) which(tolower(network$edges$to_id) == x))))
  ix <- unique(c(ix_from, ix_to))
  iy <- unique(unlist(lapply(nodes, function(x) which(tolower(network$nodes$id) == x | tolower(network$nodes$TRAITID) == x))))
  cat("Nodes matched in network$nodes:", length(iy), "\n")
  subnet <- list(edges = network$edges[ix, ], nodes = network$nodes[iy, ])
  subnet$edges <- subnet$edges[subnet$edges$from_id %in% subnet$nodes$id & subnet$edges$to_id %in% subnet$nodes$id, ]
  cat("Subnet nodes:", nrow(subnet$nodes), "Edges:", nrow(subnet$edges), "\n")
  subnet
}

processExcelSheet <- function(inFile, sheet, mand_cols, anno_cols, custom_cols) {
  data <- read_excel(path = inFile, sheet = sheet, .name_repair = "unique_quiet")
  names_lower <- tolower(names(data))
  has_traitid1 <- "traitid1" %in% names_lower
  has_traitid2 <- "traitid2" %in% names_lower
  validNetwork <- has_traitid1 && has_traitid2
  if(validNetwork) {
    if(!("pvalue" %in% names_lower)) {
      data$PVALUE <- 0.05
    }
    if(!(any("beta" %in% names_lower) || any("cor" %in% names_lower))) {
      data$BETA <- 1
    }
  }
  validAnno <- all(tolower(anno_cols) %in% names_lower)
  validCustom <- all(tolower(custom_cols) %in% names_lower)
  list(
    sheet = sheet,
    data = data,
    validNetwork = validNetwork,
    validAnno = validAnno,
    validCustom = validCustom
  )
}

processExcelFile <- function(inFile) {
  sheets <- excel_sheets(inFile)
  mand_cols <- c("TRAITID1", "TRAITID2", "PVALUE", "BETA", "COR")
  anno_cols <- c("TRAITID", "SHORTNAME", "PLAT")
  custom_cols <- c("PLAT", "COLOR", "SHAPE")
  datalist <- list()
  sheets_ok <- c()
  sheets_flag <- c()
  anno <- NULL
  custom <- NULL
  for (sheet in sheets) {
    res <- processExcelSheet(inFile, sheet, mand_cols, anno_cols, custom_cols)
    if (res$validNetwork) {
      data <- res$data
      weight_idx <- which(tolower(names(data)) %in% c("beta", "cor"))[1]
      if (!("id" %in% tolower(names(data)))) {
        data$id <- paste0(sheet, "_", seq_len(nrow(data)))
      } else {
        id_idx <- which(tolower(names(data)) == "id")
        names(data)[id_idx] <- "id"
      }
      df <- data.frame(
        to = trimws(data[["TRAITID1"]]),
        from = trimws(data[["TRAITID2"]]),
        pvalue = data[["PVALUE"]],
        weight = data[[weight_idx]],
        id = data$id,
        type = sheet,
        stringsAsFactors = FALSE
      )
      datalist[[sheet]] <- df
      sheets_ok <- c(sheets_ok, sheet)
    } else if (res$validAnno) {
      anno <- res$data %>% mutate(TRAITID = trimws(TRAITID)) %>% select(all_of(anno_cols))
    } else if (res$validCustom) {
      custom <- res$data %>% select(all_of(custom_cols))
    } else {
      missing_cols <- setdiff(mand_cols, names(res$data))
      sheets_flag <- c(sheets_flag, paste(sheet, ": missing columns -", paste(missing_cols, collapse = ", ")))
    }
  }
  if (length(datalist) == 0) {
    stop("Invalid format or missing data. Please ensure at least one sheet has the mandatory headers: TRAITID1, TRAITID2, PVALUE, and one of BETA or COR.")
  }
  if (is.null(anno)) {
    all_traits <- unique(unlist(lapply(datalist, function(df) {
      c(as.character(df$to), as.character(df$from))
    })))
    anno <- data.frame(
      TRAITID = all_traits,
      SHORTNAME = all_traits,
      PLAT = rep("Unknown", length(all_traits)),
      stringsAsFactors = FALSE
    )
  }
  list(
    datalist = datalist,
    sheets_ok = sheets_ok,
    sheets_flag = sheets_flag,
    anno = anno,
    custom = custom
  )
}

constructNetwork <- function(datalist, anno, custom) {
  all_edges <- do.call(rbind, datalist)
  cat("Initial edges:", nrow(all_edges), "Sample from/to:", head(all_edges$from, 3), "->", head(all_edges$to, 3), "\n")
  all_edges$sign <- sign(all_edges$weight)
  all_edges$weight <- abs(all_edges$weight)
  all_edges$pvalue <- as.numeric(as.character(all_edges$pvalue))
  if(any(is.na(all_edges$pvalue))) {
    warning("Some p-values could not be converted to numeric; setting them to 1E-100.")
    all_edges$pvalue[is.na(all_edges$pvalue)] <- 1E-100
  }
  all_edges$pvalue[all_edges$pvalue == 0] <- 1E-100
  edge_thickness <- -log10(all_edges$pvalue)
  if (min(edge_thickness) == max(edge_thickness)) {
    all_edges$width <- 1
  } else {
    seq_range <- round(seq(min(edge_thickness) - 1, max(edge_thickness) - 1, length.out = 10))
    if(length(unique(seq_range)) < length(seq_range)) {
      seq_range <- unique(seq_range)
    }
    all_edges$width <- as.numeric(cut(edge_thickness, breaks = seq_range, labels = FALSE))
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
      PLAT = ifelse(is.na(PLAT), "Unknown", PLAT)
    )
  all_edges$from_id <- all_edges$from
  all_edges$to_id <- all_edges$to
  all_edges$from <- all_nodes$SHORTNAME[match(tolower(all_edges$from), tolower(all_nodes$TRAITID))]
  all_edges$to <- all_nodes$SHORTNAME[match(tolower(all_edges$to), tolower(all_nodes$TRAITID))]
  all_edges <- all_edges[!is.na(all_edges$from) & !is.na(all_edges$to), ]
  cat("Edges after mapping:", nrow(all_edges), "Vitamin_D edges:", 
      sum(tolower(all_edges$from_id) == "vitamin_d" | tolower(all_edges$to_id) == "vitamin_d"), "\n")
  all_nodes$id <- all_nodes$SHORTNAME
  names(all_nodes)[names(all_nodes) == "PLAT"] <- "plat"
  if (!is.null(custom)) {
    all_platforms <- union(unique(all_nodes$plat), custom$PLAT)
    missing_platforms <- setdiff(all_platforms, custom$PLAT)
    missing_df <- data.frame(
      PLAT = missing_platforms,
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
    default_shapes <- setNames(rep("star", length(plat_list)), plat_list)
    default_colors[c("DNA", "SOMA", "ALAMAR", "BRAIN", "BM", "LD", "CPG", "CLIN",
                     "CM", "IgA", "RNA", "HD4", "IgG", "miRNA", "OLINK", "PGP",
                     "PM", "SM", "UM", "STAT", "GWAS")] <-
      c("#23bbee", "#a62281", "#a62281", "#f2921f", "#ffc815", "#ffc815",
        "#145da9", "#a0b6a8", "#57ba47", "#e41d30", "#5c2d83", "#57ba47",
        "#e41d30", "#5c2d83", "#a62281", "#e41d30", "#57ba47", "#57ba47",
        "#57ba47", "#EEEEEE", "#EEEEEE")
    default_shapes[c("UM", "CM", "SM", "DNA", "RNA", "CPG", "STAT", "GWAS")] <-
      c("square", "square", "triangle", "diamond", "diamond", "diamond", "circle", "dot")
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
  
  rv <- reactiveValues(
    customSelections = list(),
    useCases = list(),
    dataLoaded = FALSE,
    selectedTrait = NULL
  )
  networkTrigger <- reactiveVal(0)
  
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
          cat("loadData: Is Vitamin_D in traitlist?", "Vitamin_D" %in% traitlist, "\n")
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
      vitamin_d_rows <- myValues$fullnet$nodes %>% dplyr::filter(id == "Vitamin_D" | TRAITID == "Vitamin_D")
      cat("loadData: Vitamin_D rows in fullnet$nodes:\n")
      print(vitamin_d_rows)
      
      output$summary1 <- renderText({
        paste0("Total number of unique TRAIT IDs: ", length(unique(myValues$fullnet$nodes$TRAITID)))
      })
      output$summary3 <- renderText({
        paste0("Total number of associations: ", nrow(myValues$fullnet$edges))
      })
      
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
    } else {
      ix <- which(myValues$fullnet$nodes$plat == input$plat)
      traitlist <- sort(myValues$fullnet$nodes$id[ix])
    }
    cat("Updating trait list for platform", input$plat, ":", length(traitlist), "traits\n")
    cat("Sample traits:", paste(head(traitlist, 5), collapse = ", "), "\n")
    cat("Last few traits:", paste(tail(traitlist, 5), collapse = ", "), "\n")
    cat("Is Vitamin_D in traitlist?", "Vitamin_D" %in% traitlist, "\n")
    tryCatch({
      updateSelectInput(
        session,
        inputId = "trait",
        choices = traitlist,
        selected = traitlist[1]
      )
    }, error = function(e) {
      cat("Error updating trait dropdown in observe:", e$message, "\n")
      showNotification(paste("Error updating trait dropdown:", e$message), type = "error")
    })
  })
  
  # Load saved data prompt
  observe({
    if (file.exists(data_file) && file.exists(custom_settings_file) && !rv$dataLoaded) {
      showModal(modalDialog(
        title = "Load Previous Data",
        "Saved data found. Would you like to continue with the previous data?",
        footer = tagList(
          actionButton("loadSavedDataYes", "Yes"),
          actionButton("loadSavedDataNo", "No")
        )
      ))
    } else if (!rv$dataLoaded) {
      updateTabItems(session, "tabs", "Tab6")
    }
  })
  
  observeEvent(input$loadSavedDataYes, {
    tryCatch({
      saved_data <- readRDS(data_file)
      saved_settings <- readRDS(custom_settings_file)
      success <- loadData(saved_data, saved_settings)
      if (!success) {
        showNotification("Saved data is invalid or missing required columns. Please upload a new file.", type = "error")
        updateTabItems(session, "tabs", "Tab7")
      }
      removeModal()
    }, error = function(e) {
      showNotification(paste("Error loading saved data:", e$message), type = "error")
      updateTabItems(session, "tabs", "Tab7")
      removeModal()
    })
  })
  
  observeEvent(input$loadSavedDataNo, {
    if (file.exists(data_file)) file.remove(data_file)
    if (file.exists(custom_settings_file)) file.remove(custom_settings_file)
    rv$dataLoaded <- FALSE
    updateTabItems(session, "tabs", "Tab6")
    removeModal()
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
    output$summary1 <- renderText({
      paste0("Total number of unique TRAIT IDs: ", length(unique(fullnet$nodes$TRAITID)))
    })
    output$summary3 <- renderText({
      paste0("Total number of associations: ", nrow(fullnet$edges))
    })
  })
  
  # Dynamic trait dropdown
  output$SelectTrait <- renderUI({
    cat("Attempting to render SelectTrait UI\n")
    if (is.null(myValues$fullnet)) {
      cat("myValues$fullnet is NULL, skipping SelectTrait rendering\n")
      showNotification("No network data available. Please upload a valid Excel file.", type = "error")
      return(div("No traits available. Please upload data."))
    }
    traitlist <- sort(myValues$fullnet$nodes$id)
    cat("SelectTrait: Trait list length:", length(traitlist), "\n")
    cat("SelectTrait: Sample traits:", paste(head(traitlist, 5), collapse = ", "), "\n")
    cat("SelectTrait: Is Vitamin_D in traitlist?", "Vitamin_D" %in% traitlist, "\n")
    tryCatch({
      selectInput(
        inputId = "trait",
        label = "Select Trait",
        choices = traitlist,
        selected = traitlist[1]
      )
    }, error = function(e) {
      cat("Error rendering SelectTrait:", e$message, "\n")
      showNotification(paste("Error rendering trait dropdown:", e$message), type = "error")
      return(div("Error rendering trait dropdown."))
    })
  })
  
  # Debug trait selection
  observeEvent(input$trait, {
    req(input$trait, myValues$fullnet)
    cat("input$trait changed to:", input$trait, "\n")
    # Find the TRAITID for the selected trait
    selected_traitid <- myValues$fullnet$nodes$TRAITID[
      tolower(myValues$fullnet$nodes$id) == tolower(input$trait)
    ]
    if (length(selected_traitid) == 0) {
      cat("No TRAITID found for input$trait:", input$trait, "\n")
      showNotification("Selected trait not found in network.", type = "error")
      return()
    }
    cat("Setting storage$focus to TRAITID:", selected_traitid, "\n")
    storage$focus <- selected_traitid
    rv$selectedTrait <- input$trait
  })
  
  
  # Dynamic button label for Update Network
  output$updateNetworkLabel <- renderText({
    if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      paste("Focus on", input$current_node_selection)
    } else if (is.null(input$plat) || is.null(rv$selectedTrait)) {
      "Select Platform and Trait"
    } else {
      paste("Display Network for", input$plat, ":", rv$selectedTrait)
    }
  })
  
  # Handle Update Network button
  observeEvent(input$updateNetworkBtn, {
    req(myValues$fullnet)
    if (!is.null(input$current_node_selection) && length(input$current_node_selection) > 0) {
      selected_traitid <- myValues$fullnet$nodes$TRAITID[
        tolower(myValues$fullnet$nodes$id) == tolower(input$current_node_selection)
      ]
      if (length(selected_traitid) == 0) {
        showNotification("Selected node not found in network.", type = "error")
        return()
      }
      storage$focus <- selected_traitid
      showNotification(paste("Focusing network on node:", input$current_node_selection))
      cat("Focus node TRAITID:", selected_traitid, "\n")
    } else {
      req(input$trait)
      selected_traitid <- myValues$fullnet$nodes$TRAITID[
        tolower(myValues$fullnet$nodes$id) == tolower(input$trait)
      ]
      if (length(selected_traitid) == 0) {
        showNotification("Selected trait not found in network.", type = "error")
        return()
      }
      storage$focus <- selected_traitid
      showNotification(paste("Displaying network for trait:", input$trait))
      cat("Display network for TRAITID:", selected_traitid, "\n")
    }
    networkTrigger(networkTrigger() + 1)
  })
  
  # Use Case Handlers
  observeEvent(input$saveUseCaseBtn, {
    if (!is.null(storage$subnet)) {
      showModal(modalDialog(
        title = "Save Use Case",
        textInput("useCaseName", "Enter a name for this use case:",
                  value = paste("Use Case", length(rv$useCases) + 1)),
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
    rv$useCases[[input$useCaseName]] <- storage$subnet
    removeModal()
    showNotification(paste("Saved use case:", input$useCaseName))
  })
  
  output$useCaseList <- renderUI({
    if(length(rv$useCases) > 0) {
      selectInput("selectedUseCase", "Select a saved use case:",
                  choices = names(rv$useCases))
    } else {
      h4("No saved use cases.")
    }
  })
  
  observeEvent(input$loadUseCaseBtn, {
    req(input$selectedUseCase)
    saved_net <- rv$useCases[[input$selectedUseCase]]
    if (!is.null(saved_net)) {
      storage$subnet <- saved_net
      edges_df <- data.frame(
        from = saved_net$edges$from,
        to = saved_net$edges$to,
        title = saved_net$edges$title,
        color = saved_net$edges$color,
        width = saved_net$edges$width,
        stringsAsFactors = FALSE
      )
      nodes_df <- data.frame(
        id = saved_net$nodes$id,
        label = saved_net$nodes$id,
        title = saved_net$nodes$id,
        color = saved_net$nodes$color,
        shape = saved_net$nodes$shape,
        stringsAsFactors = FALSE
      )
      nodes_df$size <- 25
      storage$edges <- edges_df
      storage$nodes <- nodes_df
      output$network <- renderVisNetwork({
        visNetwork(nodes_df, edges_df, height = "800px", width = "100%") %>%
          visNodes(shadow = list(enabled = TRUE, size = 10),
                   scaling = list(min = 10, max = 30)) %>%
          visLayout(randomSeed = 4711) %>%
          visOptions(nodesIdSelection = list(enabled = TRUE, style = 'width: 1px; height: 1px;')) %>%
          visPhysics(stabilization = FALSE) %>%
          visEdges(smooth = TRUE) %>%
          visOptions(collapse = TRUE, highlightNearest = TRUE) %>%
          visEvents(select = "function(nodes) {
            Shiny.onInputChange('current_node_selection', nodes.nodes);
          }") %>%
          visEvents(click = "function(nodes) {
            Shiny.onInputChange('click', nodes.nodes[0]);
          }") %>% 
          visInteraction(navigationButtons = TRUE)
      })
      updateTabItems(session, "tabs", "Tab1")
      showNotification(paste("Loaded use case:", input$selectedUseCase))
    }
  })
  
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
    storage$subnet <- use_case
    edges_df <- data.frame(
      from = use_case$edges$from,
      to = use_case$edges$to,
      title = use_case$edges$title,
      color = use_case$edges$color,
      width = use_case$edges$width,
      stringsAsFactors = FALSE
    )
    nodes_df <- data.frame(
      id = use_case$nodes$id,
      label = use_case$nodes$id,
      title = use_case$nodes$id,
      color = use_case$nodes$color,
      shape = use_case$nodes$shape,
      stringsAsFactors = FALSE
    )
    storage$edges <- edges_df
    storage$nodes <- nodes_df
    output$network <- renderVisNetwork({
      visNetwork(nodes_df, edges_df, height = "800px", width = "100%") %>%
        visNodes(shadow = list(enabled = TRUE, size = 10),
                 scaling = list(min = 10, max = 30)) %>%
        visLayout(randomSeed = 4711) %>%
        visOptions(nodesIdSelection = list(enabled = TRUE, style = 'width: 1px; height: 1px;')) %>%
        visPhysics(stabilization = FALSE) %>%
        visEdges(smooth = TRUE) %>%
        visOptions(collapse = TRUE, highlightNearest = TRUE) %>%
        visEvents(select = "function(nodes) {
          Shiny.onInputChange('current_node_selection', nodes.nodes);
        }") %>%
        visEvents(click = "function(nodes) {
          Shiny.onInputChange('click', nodes.nodes[0]);
        }") %>% 
        visInteraction(navigationButtons = TRUE)
    })
    updateTabItems(session, "tabs", "Tab1")
    showNotification("Use case loaded successfully!")
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
      tags$img(src = "about.jpg", width = "100%")
    }
  })
  
  # Dynamic Table for Customization
  output$dynamicTable <- renderDT({
    req(myValues$fullnet)
    list_items <- unique(myValues$fullnet$nodes$plat)
    myValues$listItems <- list_items
    data <- data.frame(
      Item = list_items,
      Shape = sapply(seq_along(list_items), function(i) {
        platform <- list_items[i]
        current <- if (!is.null(rv$customSelections[[platform]]$shape)) 
          rv$customSelections[[platform]]$shape else "circle"
        opts <- c("circle", "diamond", "star", "square", "triangle")
        options_html <- paste(
          sapply(opts, function(opt) {
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
        current <- if (!is.null(rv$customSelections[[platform]]$color)) 
          rv$customSelections[[platform]]$color else "Red"
        opts <- c("Red", "Green", "Blue", "Purple", "Grey", "Yellow", "Golden", "Pink", "Black")
        options_html <- paste(
          sapply(opts, function(opt) {
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
      escape = FALSE,
      rownames = FALSE,
      options = list(
        dom = 't',
        paging = FALSE,
        drawCallback = JS("function(settings) {
          Shiny.unbindAll(this.api().table().node());
          Shiny.bindAll(this.api().table().node());
        }")
      )
    )
  }, server = FALSE)
  
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
        h4("Sheets that met the criteria"),
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
  
  # Table row selection
  observeEvent(input$submitInfo, {
    req(myValues$datalist)
    activeSheet <- input$sheetTabs
    if (!is.null(activeSheet)) {
      selected <- input[[paste0(activeSheet, "table_rows_selected")]]
      if (!is.null(selected) && length(selected) == 1) {
        sheetData <- myValues$datalist[[activeSheet]]
        focus_ids <- sheetData$from[selected]
        valid_ids <- focus_ids[tolower(focus_ids) %in% tolower(myValues$fullnet$nodes$TRAITID)]
        if (length(valid_ids) > 0) {
          storage$focus <- valid_ids
          networkTrigger(networkTrigger() + 1)
          updateTabItems(session, "tabs", "Tab1")
        } else {
          showNotification("Selected trait not found in network.", type = "error")
        }
      } else {
        showNotification("Please select exactly one row in the table.", type = "warning")
      }
    }
  })
  
  # Network rendering
  # Network rendering
  output$network <- renderVisNetwork({
    networkTrigger()
    req(myValues$fullnet)
    maxnodes <- as.numeric(input$maxnodes)
    act_trait <- if (!is.null(storage$focus) && length(storage$focus) > 0) {
      storage$focus
    } else {
      NULL
    }
    cat("Active trait for network:", act_trait, "\n")
    cat("myValues$fullnet$edges columns:", paste(names(myValues$fullnet$edges), collapse = ", "), "\n")
    cat("myValues$fullnet$nodes columns:", paste(names(myValues$fullnet$nodes), collapse = ", "), "\n")
    if (is.null(act_trait)) {
      cat("No active trait, cannot render network\n")
      showNotification("No active trait selected. Please select a trait or node.", type = "error")
      return(NULL)
    }
    neighbor_nodes <- neighbors(act_trait, myValues$fullnet)
    cat("Neighbor nodes:", length(neighbor_nodes), "IDs:", paste(head(neighbor_nodes, 5), collapse = ", "), "\n")
    tryCatch({
      subnet <- neighbor_nodes %>% 
        maxneighbors_noSTAT(myValues$fullnet, limit = maxnodes) %>% 
        nodes2network(myValues$fullnet)
      cat("Subnet nodes:", length(subnet$nodes$id), "Edges:", nrow(subnet$edges), "\n")
      cat("subnet$edges columns:", paste(names(subnet$edges), collapse = ", "), "\n")
      cat("subnet$nodes columns:", paste(names(subnet$nodes), collapse = ", "), "\n")
      cat("subnet$edges sample (first 5 rows):\n")
      print(head(subnet$edges, n = min(nrow(subnet$edges), 5)))
      cat("subnet$nodes sample (first 5 rows):\n")
      print(head(subnet$nodes, n = min(nrow(subnet$nodes), 5)))
      storage$subnet <- subnet
      edge_cols <- tolower(names(subnet$edges))
      node_cols <- tolower(names(subnet$nodes))
      from_col <- if ("from" %in% edge_cols) {
        names(subnet$edges)[which(edge_cols == "from")]
      } else if ("from_id" %in% edge_cols) {
        names(subnet$edges)[which(edge_cols == "from_id")]
      } else {
        cat("No 'from' or 'from_id' in subnet$edges columns\n")
        showNotification("Missing 'from' or 'from_id' in edge data.", type = "error")
        return(NULL)
      }
      to_col <- if ("to" %in% edge_cols) {
        names(subnet$edges)[which(edge_cols == "to")]
      } else if ("to_id" %in% edge_cols) {
        names(subnet$edges)[which(edge_cols == "to_id")]
      } else {
        cat("No 'to' or 'to_id' in subnet$edges columns\n")
        showNotification("Missing 'to' or 'to_id' in edge data.", type = "error")
        return(NULL)
      }
      id_col <- if ("id" %in% node_cols) {
        names(subnet$nodes)[which(node_cols == "id")]
      } else {
        cat("No 'id' in subnet$nodes columns\n")
        showNotification("Missing 'id' in node data.", type = "error")
        return(NULL)
      }
      edges_df <- data.frame(
        from = subnet$edges[[from_col]],
        to = subnet$edges[[to_col]],
        title = paste(subnet$edges[[from_col]], "->", subnet$edges[[to_col]]),
        color = "blue",
        width = 1,
        stringsAsFactors = FALSE
      )
      if ("pvalue" %in% edge_cols && "weight" %in% edge_cols) {
        edges_df$title <- sprintf(
          "GGM_OLINK: %s->%s, p=%.2e, beta=%.3f",
          subnet$edges[[from_col]], subnet$edges[[to_col]],
          subnet$edges$pvalue, subnet$edges$weight
        )
      }
      if ("color" %in% edge_cols) edges_df$color <- subnet$edges[[names(subnet$edges)[which(edge_cols == "color")]]]
      if ("width" %in% edge_cols) edges_df$width <- subnet$edges[[names(subnet$edges)[which(edge_cols == "width")]]]
      edges_df <- edges_df %>% distinct(from, to, .keep_all = TRUE)
      cat("edges_df created with", nrow(edges_df), "rows\n")
      cat("edges_df columns:", paste(names(edges_df), collapse = ", "), "\n")
      cat("edges_df sample:", paste(head(edges_df$from, 3), "->", head(edges_df$to, 3), collapse = ", "), "\n")
      
      # Validate node IDs in edges
      valid_nodes <- subnet$nodes[[id_col]]
      invalid_edges <- !edges_df$from %in% valid_nodes | !edges_df$to %in% valid_nodes
      if (any(invalid_edges)) {
        cat("Invalid edges detected, removing", sum(invalid_edges), "edges\n")
        edges_df <- edges_df[!invalid_edges, ]
      }
      
      nodes_df <- data.frame(
        id = subnet$nodes[[id_col]],
        label = subnet$nodes[[id_col]],
        title = subnet$nodes[[id_col]],
        color = "grey",
        shape = "circle",
        borderWidth = 1,
        stringsAsFactors = FALSE
      )
      if ("color" %in% node_cols) nodes_df$color <- subnet$nodes[[names(subnet$nodes)[which(node_cols == "color")]]]
      if ("shape" %in% node_cols) nodes_df$shape <- subnet$nodes[[names(subnet$nodes)[which(node_cols == "shape")]]]
      
      # Fix size calculation to avoid length mismatch
      matching_nodes <- myValues$fullnet$nodes$id[
        tolower(myValues$fullnet$nodes$TRAITID) %in% tolower(act_trait)
      ]
      cat("Matching nodes for act_trait:", paste(matching_nodes, collapse = ", "), "\n")
      nodes_df$size <- ifelse(
        tolower(nodes_df$id) %in% tolower(matching_nodes),
        35,
        25
      )
      # Ensure size vector is same length as nodes_df
      if (length(nodes_df$size) != nrow(nodes_df)) {
        cat("Size vector length mismatch, setting default size\n")
        nodes_df$size <- rep(25, nrow(nodes_df))
      }
      
      # Log column lengths for debugging
      cat("nodes_df column lengths: id=", length(nodes_df$id), 
          ", label=", length(nodes_df$label), 
          ", title=", length(nodes_df$title), 
          ", color=", length(nodes_df$color), 
          ", shape=", length(nodes_df$shape), 
          ", size=", length(nodes_df$size), "\n")
      
      cat("nodes_df created with", nrow(nodes_df), "rows\n")
      cat("nodes_df columns:", paste(names(nodes_df), collapse = ", "), "\n")
      storage$edges <- edges_df
      storage$nodes <- nodes_df
      visNetwork(nodes_df, edges_df, height = "800px", width = "100%") %>%
        visNodes(
          shadow = list(enabled = TRUE, size = 10),
          scaling = list(min = 10, max = 30),
          borderWidth = 1,
          borderWidthSelected = 4,
          color = list(
            border = "#000000",
            highlight = list(border = "#FFFF00"),
            hover = list(border = "#FFFF00")
          )
        ) %>%
        visLayout(randomSeed = 4711) %>%
        visOptions(
          nodesIdSelection = list(enabled = TRUE, style = 'width: 1px; height: 1px;'),
          collapse = TRUE,
          highlightNearest = TRUE
        ) %>%
        visPhysics(stabilization = FALSE) %>%
        visEdges(smooth = TRUE) %>%
        visIgraphLayout(layout = "layout_with_fr") %>%
        visEvents(
          selectNode = "function(params) { Shiny.onInputChange('current_node_selection', params.nodes); }",
          deselectNode = "function(params) { Shiny.onInputChange('current_node_selection', []); }",
          click = "function(params) { Shiny.onInputChange('click', params.nodes[0]); }"
        ) %>% 
        visInteraction(navigationButtons = TRUE)
    }, error = function(e) {
      cat("Error in network rendering:", e$message, "\n")
      showNotification(paste("Error rendering network:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Network table
  output$tbl.networktable <- DT::renderDataTable({
    req(storage$edges, storage$nodes, storage$focus)
    selected <- input$current_node_selection
    print(head(storage$edges, n = min(nrow(storage$edges), 5)))
    if (is.null(selected) || length(selected) == 0) {
      cat("No node selected, returning empty table\n")
      showNotification("No node selected.", type = "message")
      return(data.frame(Message = "No node selected"))
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
    node_mapping <- myValues$fullnet$nodes %>% select(TRAITID, id)
    relevant_mapping <- node_mapping %>% filter(tolower(TRAITID) %in% tolower(c(selected_traitid, central_traitid)))
    print(relevant_mapping)
    invalid_from <- setdiff(storage$edges$from, storage$nodes$id)
    invalid_to <- setdiff(storage$edges$to, storage$nodes$id)
    is_central <- tolower(selected_traitid) == tolower(central_traitid)
    selected_ids <- node_mapping$id[tolower(node_mapping$TRAITID) == tolower(selected_traitid)]
    central_ids <- node_mapping$id[tolower(node_mapping$TRAITID) == tolower(central_traitid)]
    central_fullnet_edges <- myValues$fullnet$edges %>% 
      filter(tolower(from_id) == tolower(central_traitid) | tolower(to_id) == tolower(central_traitid))
    print(central_fullnet_edges)
    filtered_edges <- storage$edges %>% 
      filter(
        tolower(from) %in% tolower(selected_ids) |
          tolower(to) %in% tolower(selected_ids)
      ) %>% 
      distinct(from, to, .keep_all = TRUE)
    if (nrow(filtered_edges) > 0) {
      print(head(filtered_edges, n = min(nrow(filtered_edges), 5)))
    } else {
      cat("No edges found\n")
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
          for(i in 1:nrow(customSettings)){
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
        datalist = myValues$datalist,
        sheets_ok = myValues$sheets_ok,
        sheets_flag = myValues$sheets_flag,
        anno = myValues$anno,
        custom = myValues$custom,
        fullnet = myValues$fullnet
      ), data_file)
      saveRDS(list(customSelections = rv$customSelections), custom_settings_file)
      cat("Saved data to", data_file, "and settings to", custom_settings_file, "\n")
      incProgress(1, detail = "Switching to Tables tab...")
      Sys.sleep(1)
      updateTabItems(session, "tabs", "Tab2")
    })
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
    req(input$current_node_selection)
    h4(paste("Associations for", input$current_node_selection))
  })
  
  # Network instructions
  output$networkInstructions <- renderUI({
    req(storage$subnet$nodes)
    cat("Rendering network instructions, subnet nodes:", nrow(storage$subnet$nodes), "\n")
    div(
      style = "margin-top: 10px; padding: 10px; font-size: 14px; color: #333; background-color: #f8f8f8; border-radius: 5px;",
      "To explore the network: (1) Unselect the current node by clearing the dropdown above and select a new trait from the dropdown, or (2) Click a node in the network to select and focus on it."
    )
  })
}

# ====================================================
# UI Definition
# ====================================================

js <- "Shiny.addCustomMessageHandler('change_skin', function(skin) {
  document.body.className = skin;
});"
max_nodes_list <- c(1,20,40,60,80,100,150,200)

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = textOutput("dynamicTitle")),
  dashboardSidebar(width = 150,
                   sidebarMenu(id = "tabs", selected = "Tab6",
                               menuItem("Home", tabName = "Tab6", icon = icon("th")),
                               menuItem("Data Formats", tabName = "Tab4", icon = icon("th")),
                               menuItem("Import Data", tabName = "Tab7", icon = icon("th")),
                               menuItem("Tables", tabName = "Tab2", icon = icon("th")),
                               menuItem("Your Network", tabName = "Tab1", icon = icon("dashboard")),
                               menuItem("Export", tabName = "Tab3", icon = icon("th")),
                               menuItem("Use-cases", tabName = "Tab5", icon = icon("th"))
                   )
  ),
  dashboardBody(
    tags$script(HTML("
      Shiny.addCustomMessageHandler('change_skin', function(skin) {
        var body = document.body;
        body.className = body.className.replace(/\\bskin-\\S+/g, '');
        body.className += ' skin-' + skin;
      });
    ")),
    tags$head(tags$script(js),
              tags$style(HTML("
                .custom-usecase-btn {
                  background-color: #FF5733;
                  color: white;
                  border: none;
                  border-radius: 5px;
                  padding: 10px 20px;
                  font-size: 16px;
                }
                .custom-usecase-btn:hover {
                  background-color: #E04E2F;
                }
                .focus-network-btn {
                  background-color: #33A1FF;
                  color: white;
                  border: none;
                  border-radius: 5px;
                  padding: 10px 20px;
                  font-size: 16px;
                }
                .focus-network-btn:hover {
                  background-color: #2A8BD6;
                }
                .download-usecase-btn {
                  background-color: #33A1FF;
                  color: white;
                  border: none;
                  border-radius: 5px;
                  padding: 10px 20px;
                  font-size: 16px;
                }
                .download-usecase-btn:hover {
                  background-color: #2A8BD6;
                }
                .clear-data-btn {
                  background-color: #0000FF;
                  color: white;
                  border: none;
                  border-radius: 5px;
                  padding: 10px 20px;
                  font-size: 16px;
                  margin-top: 10px;
                }
                .clear-data-btn:hover {
                  background-color: #E02F2F;
                }
                .submit-data-btn {
                  background-color: #000000;
                  color: white;
                  border: none;
                  border-radius: 5px;
                  padding: 10px 20px;
                  font-size: 16px;
                  margin-top: 10px;
                }
                .clear-data-btn:hover {
                  background-color: #E02F2F;
                }
              "))
    ),
    tabItems(
      tabItem(tabName = "Tab1",
              fluidPage(
                fluidRow(
                  column(width = 2,
                         selectInput(inputId = "plat", label = "Select Platform", choices = NULL)
                  ),
                  column(width = 6,
                         uiOutput("SelectTrait"),
                         actionButton("updateNetworkBtn", 
                                      label = textOutput("updateNetworkLabel", inline = TRUE),
                                      class = "focus-network-btn")
                  ),
                  column(width = 2,
                         selectInput(inputId = "maxnodes", label = "Max nodes",
                                     choices = max_nodes_list, selected = max_nodes_list[2])
                  )
                ),
                fluidRow(
                  column(width = 9,
                         visNetworkOutput("network")
                  ),
                  column(width = 3,
                         uiOutput("dynamic_legend")
                  )
                ),
                fluidRow(
                  column(width = 10,
                         conditionalPanel(
                           condition = "input.current_node_selection && input.current_node_selection.length > 0",
                           uiOutput("networkTableTitle"),
                           DT::dataTableOutput('tbl.networktable')
                         )
                  )
                ),
                fluidRow(
                  column(width = 12,
                         actionButton("saveUseCaseBtn", "Save as Use-Case", class = "custom-usecase-btn"),
                         downloadButton("downloadUseCase", "Download Use-Case", class = "download-usecase-btn")
                  )
                ),
                fluidRow(
                  column(width = 12,
                         box(
                           width = NULL,
                           title = "Network Instructions",
                           status = "info",
                           solidHeader = TRUE,
                           uiOutput("networkInstructions"),
                           style = "font-size: 14px; color: #333;"
                         )
                  )
                )
              )
      ),
      tabItem(tabName = "Tab2",
              fluidPage(
                fluidRow(
                  h3("Select a row from the table and click Focus network"),
                  actionButton(inputId = "submitInfo", label = "Focus network on selection", class = "focus-network-btn"),
                  uiOutput("tabs")
                )
              )
      ),
      tabItem(tabName = "Tab3",
              fluidPage(
                fluidRow(
                  column(width = 12,
                         h3("Export"),
                         p("Download or view the network data of the active view"),
                         DT::dataTableOutput("EXPORTtable")
                  )
                )
              )
      ),
      tabItem(tabName = "Tab4",
              fluidPage(
                fluidRow(
                  column(width = 12,
                         h3("HowTo"),
                         includeHTML("howto.html")
                  )
                )
              )
      ),
      tabItem(tabName = "Tab5",
              fluidPage(
                fluidRow(
                  column(width = 12,
                         h3("Saved Use Cases"),
                         #uiOutput("useCaseList"),
                         #actionButton("loadUseCaseBtn", "Load Selected Use Case", class = "custom-usecase-btn"),
                         #fileInput("uploadUseCase", "Load Use Case", accept = ".rds")
                  )
                )
              )
      ),
      tabItem(tabName = "Tab6",
              fluidPage(
                uiOutput("customHomePage")
              )
      ),
      tabItem(tabName = "Tab7",
              fluidPage(
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
                  column(width = 4,
                         box(width = 12,
                             fileInput('file_in', 'Choose Excel File with required sheets',
                                       accept = c(".xlsx"), multiple = FALSE),
                             textOutput("summary1"),
                             textOutput("summary3"),
                             uiOutput("summary2"),
                             hr(),
                             DTOutput("dynamicTable"),
                             hr(),
                             uiOutput("SHEETlist"),
                             uiOutput("confirm_btn"),
                             uiOutput("downloadButtonUI"),
                             actionButton("clearData", "Clear Saved Data", class = "clear-data-btn")
                         )
                  ),
                  column(width = 5,
                         uiOutput("SAMPLEtable")
                  )
                ),
                fluidRow(
                  column(width = 2,
                         uiOutput("mySelection")
                  )
                ),
                fluidRow(
                  column(4, 
                         fileInput("homePageFile", 
                                   "Upload Custom Home Page (HTML or Markdown)", 
                                   accept = c(".html", ".md"))
                  )
                ),
                fluidRow(
                  column(12, 
                         helpText("Upload a custom HTML or Markdown file. The file will be rendered automatically in the Home tab. If nothing is uploaded, a default image is shown.")
                  )
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