rm(list = ls())

if(!require(tidyverse)) install.packages(tidyverse, repos = "https://cloud.r-project.org", quiet = TRUE)
if(!require(sn)) install.packages(sn, repos = "https://cloud.r-project.org", quiet = TRUE)

load('mapping_objects.rdata')

age.df <- age.df |> 
  mutate(label = paste0(crc.id, ", ", swfsc.id)) |> 
  arrange(crc.id, swfsc.id)

ui <- fluidPage(
  titlePanel("Age Confidence Rating Mapping"),
  sidebarLayout(
    sidebarPanel(
      p("Change the sliders to choose a percentage of the range between the Uniform distribution (solid grey line) and full Skew Normal distribution (yellow line) for each confidence rating"),
      p("You can click on the slider and use left and right arrow keys to change it by increments of 0.01"),
      sliderInput("cr4", "Confidence rating 4", min = 0.01, max = 0.99, value = 0.75, step = 0.01, ticks = FALSE),
      sliderInput("cr3", "Confidence rating 3", min = 0.01, max = 0.99, value = 0.5, step = 0.01, ticks = FALSE),
      sliderInput("cr2", "Confidence rating 2", min = 0.01, max = 0.99, value = 0.25, step = 0.01, ticks = FALSE),
      sliderInput("cr1", "Confidence rating 1", min = 0.01, max = 0.99, value = 0.05, step = 0.01, ticks = FALSE),
      actionButton("reset", "Reset")
    ),
    mainPanel(
      wellPanel(
        column(8, selectInput(
          'cr.filter', 
          'Filter individuals by confidence rating',
          choices = c('All', as.character(sort(unique(age.df$age.confidence)))),
          selected = 'All',
          selectize = FALSE
        )),
        column(8, selectInput("id", "Individual", choices = list(1), selected = 1)),
        fluidRow(
          column(3, actionButton("prev", "Prev")), 
          column(3, actionButton("nxt", "Next"))
        )
      ),
      wellPanel(
        p("Solid line is the actual confidence rating of the selected individual"),
        plotOutput("distPlot", height = 600),
        sliderInput("xlim", "Age limits", min = 0, max = 80, value = c(0, 80), width = "100%")
      )
    )
  )
)

server <- function(input, output, session) {
  df <- reactive({
    if(input$cr.filter == 'All') {
      age.df 
    } else {
      filter(age.df, age.confidence == as.numeric(input$cr.filter))
    }
  })
  
  p.vec <- reactive(c(input$cr1, input$cr2, input$cr3, input$cr4, 1))
  
  x.lims <- reactive(input$xlim)
  
  labid <- reactive({
    i <- as.numeric(input$id)
    updateSliderInput(session, "xlim", value = c(df()$age.min[i], df()$age.max[i]))
    i
  })
  
  observeEvent(input$reset, {
    updateSliderInput(session, "cr4", value = 0.75)
    updateSliderInput(session, "cr3", value = 0.5)
    updateSliderInput(session, "cr2", value = 0.25)
    updateSliderInput(session, "cr1", value = 0.05)
  })
  
  observe({
    updateSelectInput(
      session = session, 
      inputId = "id", 
      choices = setNames(as.list(1:nrow(df())), df()$label)
    )
  })

  observeEvent(input$prev, {
    i <- as.numeric(input$id) - 1
    if(i > 1) updateSelectInput(session, inputId = "id", selected = as.character(i))
  })
  observeEvent(input$nxt, {
    i <- as.numeric(input$id) + 1
    if(i < nrow(df())) updateSelectInput(session, inputId = "id", selected = as.character(i))
  })
  
  output$distPlot <- renderPlot({
    i <- labid()
    wt <- p.vec()
    
    id.df <- df()[i, ]
    age.vec <- seq(id.df$age.min, id.df$age.max, length.out = 1000)
    dens.df <- sapply(wt, function(x) {
      cwsnDensity(
        age = age.vec, 
        age.range = id.df$age.range, 
        dp = c(id.df$sn.location, id.df$sn.scale, id.df$sn.shape), 
        wt = x
      )
    }) |> 
      cbind(age.vec) |> 
      as.data.frame() |>
      setNames(c(1:length(wt), "age")) |> 
      pivot_longer(-age, names_to = "Confidence", values_to = "density") |> 
      mutate(
        this.cr = Confidence == as.character(id.df$age.confidence),
        Confidence = as.factor(Confidence)
      ) |> 
      arrange(Confidence, age)
    
    ggplot(dens.df, aes(age, density)) +
      geom_hline(yintercept = 1 / id.df$age.range, color = "grey") +
      geom_vline(xintercept = id.df$age.best, linewidth = 1.5) +
      geom_vline(
        xintercept = unlist(id.df[, c("age.min", "age.max")]), 
        linetype = "dashed"
      ) +
      geom_line(
        aes(color = Confidence, linetype = ifelse(this.cr, "dashed", "solid")), 
        linewidth = 2
      ) +
      scale_color_manual(values = conf.colors) +
      guides(linetype = "none") +
      labs(x = "Age", y = "Density") +
      xlim(x.lims()) +
      theme_minimal() +
      theme(
        legend.position = "top",
        text = element_text(size = 14)
      )
  })
}

shinyApp(ui = ui, server = server)