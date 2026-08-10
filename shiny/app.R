# Shiny app to run BPVS-FE or BPVS-RE predictions

library(shiny)

# Attempt to use package functions when installed, otherwise source from repo
if (requireNamespace("BPVS", quietly = TRUE)) {
  load_bpvs_model <- BPVS::load_bpvs_model
  predict_bpvs <- BPVS::predict_bpvs
} else {
  source(file.path("..", "R", "load_model_and_predict.R"), local = TRUE)
}

model_cache <- list(FE = NULL, RE = NULL)

ui <- fluidPage(
  titlePanel("BPVS: FE / RE model prediction"),
  sidebarLayout(
    sidebarPanel(
      selectInput("type", "Model type:", choices = c("FE", "RE"), selected = "FE"),
      fileInput("file", "Upload CSV with features (header, comma)", accept = c(".csv")),
      checkboxInput("show_head", "Show head of uploaded data", TRUE),
      actionButton("go", "Predict"),
      downloadButton("download", "Download results (CSV)")
    ),
    mainPanel(
      verbatimTextOutput("status"),
      tableOutput("sample"),
      tableOutput("preds"),
      verbatimTextOutput("selected")
    )
  )
)

server <- function(input, output, session) {
  data_in <- reactive({
    req(input$file)
    read.csv(input$file$datapath, stringsAsFactors = FALSE)
  })

  output$status <- renderText({
    # try load model for selected type if not cached
    type <- input$type
    if (is.null(model_cache[[type]])) {
      model_try <- tryCatch({ load_bpvs_model(type) }, error = function(e) e)
      if (inherits(model_try, "error")) {
        "Model file not found. Place inst/extdata/bpvs_fe_model.rds and/or bpvs_re_model.rds or install the package."
      } else {
        model_cache[[type]] <<- model_try
        paste("Model", type, "loaded.")
      }
    } else {
      paste("Model", type, "cached.")
    }
  })

  output$sample <- renderTable({
    req(input$file)
    if (isTRUE(input$show_head)) head(data_in(), 5) else NULL
  })

  results <- eventReactive(input$go, {
    req(input$file)
    df <- data_in()
    type <- input$type
    model_local <- model_cache[[type]]
    if (is.null(model_local)) {
      model_local <- tryCatch(load_bpvs_model(type), error = function(e) stop(e$message))
    }
    res <- predict_bpvs(newdata = df, type = type, model = model_local)
    # Normalize results into a data.frame for display/download
    preds <- res$predictions
    if (is.data.frame(preds)) outdf <- preds else outdf <- data.frame(prediction = preds)
    attr(outdf, "selected") <- res$selected
    outdf
  })

  output$preds <- renderTable({
    req(results())
    results()
  })

  output$selected <- renderText({
    req(results())
    sel <- attr(results(), "selected")
    if (is.null(sel)) "No selected covariates stored in model." else paste("Selected:", paste(sel, collapse = ", "))
  })

  output$download <- downloadHandler(
    filename = function() paste0("bpvs_results_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- results()
      # include selected covariates as a header comment
      sel <- attr(df, "selected")
      if (!is.null(sel)) writeLines(paste("# selected:", paste(sel, collapse = ", ")), con = file)
      write.table(df, file = file, sep = ",", row.names = FALSE, col.names = TRUE, append = !is.null(sel))
    }
  )
}

shinyApp(ui, server)
