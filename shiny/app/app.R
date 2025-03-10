# packages <- c("readxl","writexl","ggsci","MASS","dplyr","openxlsx","readr","magrittr","purrr","lubridate","rlang","scales","DT","tidyverse","shinydashboard","shiny")
# lapply(packages, require, character.only = TRUE) 

library(dplyr)
library(tidyverse)
library(shiny)
library(shinydashboard)
library(DT)
library(JuliaCall)

showtext::showtext_auto(enable = TRUE)
options(dplyr.summarise.inform=F)
julia_setup()
julia_command(paste0("include(\"", getwd(), "/functions.jl\")"))

# UI ----
# UI ----
ui <- dashboardPage(
  dashboardHeader(title = "Power calculation for test-negative design",
                  titleWidth = 500),
  dashboardSidebar(disable = TRUE), 
  dashboardBody(
    fluidRow(
        tags$style(HTML("
        .main-header .logo {
        font-weight: bold !important;
        background-color: #252525 !important;
        }
        .main-header .navbar {
        background-color: #4d4d4d !important;
        }
        .main-sidebar {
        background-color: #222 !important; 
        }
        .content-wrapper {
        padding-left: 30px !important;
        padding-right: 30px !important;
        background-color: #878787 !important; 
      　}
        .box-title {
        font-weight: bold; 
        }
        input[type='number']::-webkit-outer-spin-button,
        input[type='number']::-webkit-inner-spin-button {
          -webkit-appearance: none;
          margin: 0;
        }
      ")),
    fluidRow(
      box(
        title = "Assumption of proportion and VE of each vaccine status",
        width = 12,
        DTOutput("input_table") 
      )),
    fluidRow(
      box(
        title = "Number of simulations",
        width = 3,
        numericInput("n_sim", h4(""), value = NA, min = 1, max = 1000)
      )),
    fluidRow(
      box(
        title = "Number of recruitments",
        width = 3,
        numericInput("n_rec", h4(""), value = NA, min = 1, max = 50000)
      )),
    fluidRow(
      box(
        title = "Diagnostic test",
        width = 3,
        numericInput("sens", h4("Sensitivity"), value = NA, min = 0, max = 1),
        numericInput("spec", h4("Specificity"), value = NA, min = 0, max = 1),
      )),
    fluidRow(
      box(
        title = "Probability of TD among individuals presenting with TD-related symptoms",
        width = 7,
        numericInput("prob_case", h4(""), value = NA, min = 0, max = 1)
      )),
    fluidRow(
        column(12, div(
          style = "margin-bottom: 15px;",
          actionButton("run", "Run", 
                       style = "background-color: #252525; color: white;")
        ))
      ),
    fluidRow(
        box(
          title = "Estimated Power",
          width = 12,
          plotOutput("myplot")
        )
      )
    )
  )
)

# Server ----
server <- function(input, output, session) {

  data <- reactiveVal(data.frame(
    Category = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7"),
    Description = c("Control", "-6 months, Dose 3", "6-12 months, Dose 3", "12-24 months, Dose 3","-6 months, Dose 4", "6-12 months, Dose 4", "12 months-, Dose 4"),
    Proportion = c(0.328, 0.026, 0.026, 0.045, 0.063, 0.063, 0.449),
    VE = c(0, 0.5, 0.3, 0.1, 0.5, 0.3, 0.1)
  ))
  
  output$input_table <- renderDT({
    datatable(
      data(),
      editable = TRUE, 
      options = list(dom = 't') 
    )
  })
  
  input_ready <- reactiveVal(FALSE)
  observe({
    if (!is.null(input$n_sim) && !is.null(input$n_rec) &&
        !is.null(input$sens) && !is.null(input$spec) && !is.null(input$prob_case)) {
      input_ready(TRUE)
    }
  })
  
  observeEvent(input$input_table_cell_edit, {
    info <- input$input_table_cell_edit
    updated_data <- data()
    updated_data[info$row, info$col] <- as.numeric(info$value)  
    data(updated_data)  
  })
  
  observeEvent(input$run, {
    req(input_ready())
    output$myplot <- renderPlot({
      n_iter <- 100
      n_sim <- input$n_sim
      n_rec <- input$n_rec
      prob_case <- input$prob_case
      sens <- input$sens
      spec <- input$spec
      
      df <- data()
      
      prop_neg_g <- df %>% pull(Proportion)
      ve_g <- df %>% pull(VE)
      cat <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7")
      
      withProgress(message = "Calculating", value = 0, {
        
        res <- lapply(1:n_sim, function(x) {
          Sys.sleep(0.1)
          setProgress(value = x/n_sim, 
                      detail = sprintf("%d/%d", x, n_sim))
          julia_call("get_est", as.integer(n_rec), prob_case, prop_neg_g, ve_g, sens, spec, "true")
        })
        
        df_est <- do.call(rbind, res)
      })
      
      df_power <- df_est %>% 
        group_by(n, group, null) %>% 
        mutate(sdif_bs = case_when(sdif_bs == "1"|sdif_bs == "TRUE" ~ 1,
                                   sdif_bs == "0"|sdif_bs == "FALSE" ~ 0)) %>% 
        mutate(power_true = sum(l_true > 0)/n_sim,
               power_obs = sum(l_obs > 0)/n_sim,
               power_cor = sum(sdif_bs)/n_sim
        ) %>% 
        distinct(n, group, .keep_all = T) %>% 
        dplyr::select(n, group, 13:16) %>% 
        gather(4:6, key = type, value = power)
      
      g_power <- 
        ggplot()+
        geom_bar(data = df_power, aes(x = group, y = power, fill = type), 
                 stat = "identity",
                 position = "dodge",
                 colour = "black")+
        scale_fill_manual(values = c("#f46d43","#3288bd","#66c2a5"), labels = c('Corrected', 'Observed', 'True'))+
        labs(x = "Category", y = "Power") +
        geom_hline(yintercept = 0.8, linetype = 2)+
        theme_bw()+
        theme(panel.grid.major.y = element_blank(),
              axis.text = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 12, face = "bold"),
              legend.title = element_blank(),
              axis.title = element_text(size = 15, face = "bold")
        )
      
      g_power
    })
  })
}

# Run app ----
shinyApp(ui,server)





