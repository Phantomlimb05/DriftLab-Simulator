#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

library(shiny)

# ---------------------------
# Packages
# ---------------------------
if (!requireNamespace("shiny", quietly = TRUE)) stop("Please install 'shiny'")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'")
if (!requireNamespace("DT", quietly = TRUE)) stop("Please install 'DT'")

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ---------------------------
# Core simulator
# ---------------------------
# One generation update for allele A frequency p under selection, mutation, migration
wf_step <- function(p, N, s = 0, h = 0.5, u = 0, v = 0, m = 0, p_mig = 0) {
  # Selection: fitnesses
  w_AA <- 1 + s
  w_Aa <- 1 + h * s
  w_aa <- 1
  # After selection (deterministic change in gamete pool)
  mean_w <- p^2 * w_AA + 2 * p * (1 - p) * w_Aa + (1 - p)^2 * w_aa
  p_sel  <- (p^2 * w_AA + p * (1 - p) * w_Aa) / mean_w
  # Mutation (A -> a with rate u; a -> A with rate v)
  p_mut  <- (1 - u) * p_sel + v * (1 - p_sel)
  # Migration (fraction m replaced by migrants at frequency p_mig)
  p_migr <- (1 - m) * p_mut + m * p_mig
  # Genetic drift (sampling of 2N gametes)
  k <- rbinom(1, size = 2 * N, prob = p_migr)
  p_next <- k / (2 * N)
  p_next
}

# Simulate a full trajectory across G generations
simulate_wf <- function(N, p0, G, s = 0, h = 0.5, u = 0, v = 0, m = 0, p_mig = 0) {
  p <- numeric(G + 1)
  p[1] <- p0
  if (p0 < 0 || p0 > 1) stop("p0 must be in [0,1]")
  if (N <= 0 || G <= 0) stop("N and G must be positive")
  for (t in 1:G) {
    p[t + 1] <- wf_step(p[t], N = N, s = s, h = h, u = u, v = v, m = m, p_mig = p_mig)
  }
  p
}

# Simulate R replicates; return long tibble
simulate_replicates <- function(R, N, p0, G, s, h, u, v, m, p_mig) {
  reps <- vector("list", R)
  for (i in seq_len(R)) {
    reps[[i]] <- tibble(
      replicate = i,
      generation = 0:G,
      p = simulate_wf(N, p0, G, s, h, u, v, m, p_mig)
    )
  }
  bind_rows(reps)
}

# Compute fixation stats at generation G
fixation_summary <- function(df_long) {
  last <- df_long %>% group_by(replicate) %>% slice_tail(n = 1) %>% ungroup()
  n <- nrow(last)
  fixA <- sum(last$p >= 1 - 1e-12)
  fixa <- sum(last$p <= 1e-12)
  tibble(
    replicates = n,
    fixed_A = fixA,
    fixed_a = fixa,
    seg_segregating = n - fixA - fixa,
    P_fix_A = fixA / n,
    P_fix_a = fixa / n
  )
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("body { max-width: 1200px; margin: 0 auto; } .small-note { font-size: 12px; color: #666; }"))),
  titlePanel("Wright–Fisher Model Simulator"),
  sidebarLayout(
    sidebarPanel(
      h4("Population & Time"),
      numericInput("N", "Diploid population size (N)", value = 100, min = 10, step = 10),
      sliderInput("G", "Generations", min = 10, max = 5000, value = 200, step = 10),
      sliderInput("p0", "Initial freq of allele A (p0)", min = 0, max = 1, value = 0.5, step = 0.01),
      h4("Forces"),
      sliderInput("s", "Selection coefficient (s) for AA vs aa", min = -0.5, max = 0.5, value = 0, step = 0.01),
      sliderInput("h", "Dominance (h): Aa fitness = 1 + h*s", min = 0, max = 1, value = 0.5, step = 0.05),
      sliderInput("u", "Mutation A→a (u)", min = 0, max = 0.1, value = 0, step = 0.0001),
      sliderInput("v", "Mutation a→A (v)", min = 0, max = 0.1, value = 0, step = 0.0001),
      sliderInput("m", "Migration fraction (m)", min = 0, max = 0.5, value = 0, step = 0.01),
      sliderInput("p_mig", "Migrant allele frequency (p_mig)", min = 0, max = 1, value = 0.5, step = 0.01),
      h4("Replicates"),
      sliderInput("R", "Number of replicates", min = 1, max = 500, value = 50, step = 1),
      numericInput("seed", "Random seed (optional)", value = 123, min = 0, step = 1),
      actionButton("run", "Run simulation", class = "btn-primary"),
      hr(),
      downloadButton("download_csv", "Download trajectories (CSV)"),
      div(class = "small-note", "All parameters apply to every replicate. Plot caps at 500 trajectories for readability.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Trajectories",
                 br(),
                 plotOutput("traj_plot", height = "420px")),
        tabPanel("Summary",
                 br(),
                 fluidRow(
                   column(6, tableOutput("fix_table")),
                   column(6, plotOutput("end_hist", height = "300px"))
                 ),
                 br(),
                 DTOutput("rep_table")
        ),
        tabPanel("Help",
                 br(),
                 HTML("<h4>Model details</h4>
                 <ul>
                   <li><b>Wright–Fisher</b>: Non-overlapping generations, random mating, genetic drift via binomial sampling of 2N gametes.</li>
                   <li><b>Selection</b>: Fitnesses AA=1+s, Aa=1+h·s, aa=1; selection applied before drift.</li>
                   <li><b>Mutation</b>: A→a at rate u; a→A at rate v; applied after selection.</li>
                   <li><b>Migration</b>: Fraction m of the gene pool replaced by migrants with frequency p<sub>mig</sub>.</li>
                 </ul>
                 <p>Order of forces per generation: selection → mutation → migration → drift.</p>
                 <h4>Tips</h4>
                 <ul>
                   <li>Set s=0,u=v=0,m=0 to visualise pure drift.</li>
                   <li>Use larger N to slow drift; smaller N speeds fixation/loss.</li>
                   <li>Positive s increases fixation probability of A (dependent on h and p0).</li>
                 </ul>
                 ")
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  # Run simulations when button is clicked
  sims <- eventReactive(input$run, {
    validate(
      need(input$N > 0 && input$N == as.integer(input$N), "N must be a positive integer."),
      need(input$G > 0, "Generations must be positive."),
      need(input$p0 >= 0 && input$p0 <= 1, "p0 must be in [0,1].")
    )
    if (!is.null(input$seed) && !is.na(input$seed)) set.seed(input$seed)
    df <- simulate_replicates(
      R = input$R,
      N = as.integer(input$N),
      p0 = input$p0,
      G = as.integer(input$G),
      s = input$s,
      h = input$h,
      u = input$u,
      v = input$v,
      m = input$m,
      p_mig = input$p_mig
    )
    df
  }, ignoreInit = TRUE)
  
  # Trajectory plot
  output$traj_plot <- renderPlot({
    df <- sims()
    req(df)
    # Limit the number of drawn trajectories for readability
    R_draw <- min(length(unique(df$replicate)), 500)
    df_draw <- df %>% filter(replicate %in% seq_len(R_draw))
    ggplot(df_draw, aes(x = generation, y = p, group = replicate)) +
      geom_line(alpha = 0.5) +
      geom_hline(yintercept = c(0,1), linetype = "dashed") +
      labs(x = "Generation", y = "Allele A frequency (p)",
           title = sprintf("%d trajectories (N=%s, p0=%.2f, s=%.2f, h=%.2f)",
                           R_draw, format(input$N, big.mark=","), input$p0, input$s, input$h)) +
      theme_minimal(base_size = 13)
  })
  
  # Fixation summary table
  output$fix_table <- renderTable({
    df <- sims(); req(df)
    fixation_summary(df)
  }, digits = 4)
  
  # Histogram of final frequencies
  output$end_hist <- renderPlot({
    df <- sims(); req(df)
    last <- df %>% group_by(replicate) %>% slice_tail(n=1)
    ggplot(last, aes(x = p)) +
      geom_histogram(bins = 30) +
      labs(x = "Final allele A frequency", y = "Count",
           title = "Distribution of final frequencies (last generation)") +
      theme_minimal(base_size = 12)
  })
  
  # Replicate table (final state)
  output$rep_table <- renderDT({
    df <- sims(); req(df)
    last <- df %>% group_by(replicate) %>% slice_tail(n=1) %>% ungroup() %>%
      transmute(replicate, final_p = p,
                state = case_when(p <= 1e-12 ~ "Lost (a fixed)",
                                  p >= 1 - 1e-12 ~ "Fixed (A fixed)",
                                  TRUE ~ "Segregating"))
    datatable(last, rownames = FALSE, options = list(pageLength = 10))
  })
  
  # CSV download
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("wf_trajectories_N", input$N, "_R", input$R, "_G", input$G, ".csv")
    },
    content = function(file) {
      df <- sims(); req(df)
      readr::write_csv(df, file)
    }
  )
}

shinyApp(ui, server)
