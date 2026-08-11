# ======================================================================
# BPVS Shiny App
# Bayesian Panel Variable Selection -- point-and-click interface
# Pathairat Pastpipatkul and Htwe Ko (2025)
# ======================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(BPVS)
  library(dplyr)
  library(ggplot2)
})

ui <- fluidPage(
  titlePanel("Bayesian Panel Variable Selection (BPVS)"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Data"),
      radioButtons("datasrc", "Data source:",
                   choices = c("Use bundled 'developing' example" = "bundled",
                               "Upload CSV" = "upload")),
      conditionalPanel(
        condition = "input.datasrc == 'upload'",
        fileInput("file", "Choose CSV file",
                  accept = c(".csv", "text/csv"))
      ),
      selectInput("outcome", "Outcome variable (y)", choices = character(0)),
      selectInput("idvar", "Cross-sectional id variable", choices = character(0)),
      selectInput("timevar", "Time variable", choices = character(0)),
      h4("Model"),
      checkboxGroupInput("effects", "Specifications:",
                         choices = c("Fixed effects (unit dummies)" = "fe",
                                     "Random effects (random intercept)" = "re"),
                         selected = c("fe", "re")),
      numericInput("iter", "Iterations", value = 2000, min = 200, step = 100),
      numericInput("warmup", "Warm-up", value = 1000, min = 100, step = 100),
      numericInput("chains", "Chains", value = 2, min = 1, max = 4, step = 1),
      numericInput("seed", "Seed", value = 123, step = 1),
      actionButton("run", "Run BPVS", class = "btn-primary",
                   icon = icon("play")),
      br(), br(),
      verbatimTextOutput("status")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Rankings",
                 h4("BPVS-FE"),
                 tableOutput("rank_fe"),
                 h4("BPVS-RE"),
                 tableOutput("rank_re")),
        tabPanel("Inclusion Probabilities",
                 plotOutput("pip_fe", height = "380px"),
                 plotOutput("pip_re", height = "380px")),
        tabPanel("Coefficient Estimates",
                 plotOutput("coef_fe", height = "380px"),
                 plotOutput("coef_re", height = "380px")),
        tabPanel("FE vs RE & Hausman",
                 plotOutput("pip_compare", height = "380px"),
                 h4("Bayesian Hausman"),
                 tableOutput("bayes_hausman"),
                 h4("Frequentist Hausman (plm)"),
                 tableOutput("hausman_plm")),
        tabPanel("Fit",
                 verbatimTextOutput("fit_stats"))
      )
    )
  )
)

server <- function(input, output, session) {

  values <- reactiveValues(
    data = NULL,
    res  = NULL,
    err  = NULL
  )

  observeEvent(input$datasrc, {
    if (input$datasrc == "bundled") {
      values$data <- BPVS::developing
    }
  })

  observeEvent(input$file, {
    req(input$file)
    dat <- tryCatch(
      read.csv(input$file$datapath, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(dat)) {
      values$err <- "Could not read the uploaded file."
      values$data <- NULL
    } else {
      values$err <- NULL
      values$data <- dat
    }
  })

  observe({
    dat <- values$data
    if (is.null(dat)) return()
    num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
    chr_cols <- names(dat)[!vapply(dat, is.numeric, logical(1))]
    updateSelectInput(session, "outcome",
                      choices = setdiff(num_cols, chr_cols))
    updateSelectInput(session, "idvar", choices = names(dat),
                      selected = if ("id" %in% names(dat)) "id" else names(dat)[1])
    updateSelectInput(session, "timevar", choices = names(dat),
                      selected = if ("year" %in% names(dat)) "year" else NULL)
  })

  observeEvent(input$run, {
    req(values$data, input$outcome, input$idvar)
    values$err <- NULL
    values$res <- NULL

    dat <- values$data
    predictors <- setdiff(
      names(dat)[vapply(dat, is.numeric, logical(1))],
      c(input$outcome, input$idvar, input$timevar)
    )
    if (length(predictors) == 0) {
      values$err <- "No numeric predictors available."
      return()
    }

    dat_scaled <- dat
    for (v in predictors) {
      m <- mean(dat_scaled[[v]], na.rm = TRUE)
      s <- sd(dat_scaled[[v]], na.rm = TRUE)
      dat_scaled[[v]] <- (dat_scaled[[v]] - m) / s
    }

    fmla <- stats::as.formula(
      paste(input$outcome, "~", paste(predictors, collapse = " + "))
    )

    values$status <- "Fitting (this can take a while)..."
    tryCatch({
      values$res <- BPVS::bpvs_pipeline(
        formula     = fmla,
        data        = dat_scaled,
        effects     = input$effects,
        outcome_var = input$outcome,
        id_var      = input$idvar,
        time_var    = input$timevar,
        chains      = input$chains,
        iter        = input$iter,
        warmup      = input$warmup,
        cores       = 1,
        seed        = input$seed,
        refresh     = 0
      )
      values$status <- "Done."
    }, error = function(e) {
      values$err <- paste("Fitting failed:", conditionMessage(e))
      values$status <- NULL
    })
  })

  output$status <- renderPrint({
    if (!is.null(values$err)) values$err else values$status
  })

  output$rank_fe <- renderTable({
    req(values$res$rank_fe)
    values$res$rank_fe[, c("Rank", "Variable", "Estimate", "PIP", "Selection")]
  }, digits = 4)

  output$rank_re <- renderTable({
    req(values$res$rank_re)
    values$res$rank_re[, c("Rank", "Variable", "Estimate", "PIP", "Selection")]
  }, digits = 4)

  output$pip_fe <- renderPlot({
    req(values$res$rank_fe)
    BPVS::plot_inclusion(values$res$rank_fe, title = "BPVS-FE: PIPs")
  })
  output$pip_re <- renderPlot({
    req(values$res$rank_re)
    BPVS::plot_inclusion(values$res$rank_re, title = "BPVS-RE: PIPs")
  })

  output$coef_fe <- renderPlot({
    req(values$res$rank_fe)
    BPVS::plot_coefficients(values$res$rank_fe, title = "BPVS-FE: Coefficient Estimates")
  })
  output$coef_re <- renderPlot({
    req(values$res$rank_re)
    BPVS::plot_coefficients(values$res$rank_re, title = "BPVS-RE: Coefficient Estimates")
  })

  output$pip_compare <- renderPlot({
    req(values$res$rank_fe, values$res$rank_re)
    BPVS::plot_pip_compare(values$res$rank_fe, values$res$rank_re)
  })

  output$bayes_hausman <- renderTable({
    req(values$res$comparison)
    values$res$comparison$bayesian_hausman
  }, digits = 4)

  output$hausman_plm <- renderTable({
    req(values$res$comparison$hausman_plm)
    as.data.frame(values$res$comparison$hausman_plm)
  }, digits = 4)

  output$fit_stats <- renderPrint({
    req(values$res$comparison)
    cat("FE: WAIC =", values$res$comparison$stats_fe$WAIC,
        "(se", values$res$comparison$stats_fe$WAIC_SE, "),",
        "Bayes R2 =", values$res$comparison$stats_fe$Bayes_R2, "\n")
    cat("RE: WAIC =", values$res$comparison$stats_re$WAIC,
        "(se", values$res$comparison$stats_re$WAIC_SE, "),",
        "Bayes R2 =", values$res$comparison$stats_re$Bayes_R2, "\n")
  })
}

shinyApp(ui = ui, server = server)
