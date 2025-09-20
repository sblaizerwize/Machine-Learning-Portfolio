# Load packages ----------------------------------------------------------------
library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(survival)

# Define UI --------------------------------------------------------------------
ui <- page_sidebar(
  
  title = "ANÁLISIS DE PACIENTES CON CÁNCER GASTROINTESTINAL | PEC4 | SAD",
  
  sidebar = sidebar(
    # Carga del archivo
    fileInput("file", "Carga Archivo CSV"),
    
    hr(), 
    p("Plot Supervivencia vs. Variables Clínicas"),
    
    # Select variable for y-axis
    selectInput(
      inputId = "y",
      label = tags$p("Eje Y:", style = "font-size: 80%; color: gray;"),
      choices = NULL
    ),
    
    # Select variable for x-axis
    selectInput(
      inputId = "x",
      label = tags$p("Eje X:", style = "font-size: 80%; color: gray;"),
      choices = NULL
    ),
    
    hr(), 
    p("Plot Distribución de Edades e IMC"),
    
    # Select variable for histogram
    selectInput(
      inputId = "z",
      label = tags$p("Categoría:", style = "font-size: 80%; color: gray;"),
      choices = NULL
    ),
    
    hr(), 
    p("Plot Supervivencia Kaplan Mier"),
    
    # Select variable for Kaplan Mier
    selectInput(
      inputId = "w",
      label = tags$p("Categoría:", style = "font-size: 80%; color: gray;"),
      choices = NULL
    ),
  ),
  
  layout_columns(
    card(
      card_header("Resumen Estadístico"),
      verbatimTextOutput(outputId = "summary")
    ),
    card(
      card_header("Supervivencia de Kaplan Meier"),
      plotOutput(outputId = "probability")
    ),
    card(
      card_header("Distribución de Edades e IMC"),
      plotOutput(outputId = "histogram")
    ),
    card(
      card_header("Supervivencia según variables clínicas"),
      plotOutput(outputId = "survival")
    ),
    col_widths = c(5, 7, 5, 7),  # c(8, 4, 5, 7)
    row_heights = c(1, 1)        # c(1, 1)
  )
)

# Define server ----------------------------------------------------------------
server <- function(input, output, session) {
  
  # Definir función calcular_imc()
  calcular_imc <- function(peso, altura) {
    peso / (altura / 100)^2
  }
  
  # Lectura reactiva del archivo CSV
  tcga_gastro <- reactive({
    req(input$file)
    read.csv(input$file$datapath, header=TRUE, sep = ";")
  })
  
  # Transformación reactiva de datos y asignación a df_gastro
  df_gastro <- reactive({
    # Llamar objeto reactivo tcga_gastro()
    req(tcga_gastro()) 
    # Asignar transformaciones
    tcga_gastro() %>%
      # Seleccionar columnas de interés
      select(gender, vital_status, days_to_death, days_to_last_followup, 
             age_at_initial_pathologic_diagnosis, weight, height, venous_invasion, 
             lymphatic_invasion, perineural_invasion_present, cancer_type, stage_event_pathologic_stage)  %>%
      # Definir os_time
      mutate(
        os_time = ifelse(
          !is.na(days_to_death),
          days_to_death,
          days_to_last_followup)) %>%
      # Redefinir os_time
      mutate(
        os_time = ifelse(is.infinite(os_time), NA, os_time)
      ) %>%
      # Definir vital_status
      mutate(
        vital_status = ifelse((is.na(vital_status) & !is.na(days_to_last_followup)), "Alive", vital_status)
      ) %>%
      # Definir os_event
      mutate(
        os_event = ifelse(vital_status == "Dead", 1, 0) 
      ) %>%
      # Transformar cancer_type
      mutate(
        cancer_type = gsub("TCGA-","", cancer_type)
      ) %>%
      # Transformar stage_event_pathologic_stage
      mutate(
        stage_event_pathologic_stage = gsub("Stage ","", stage_event_pathologic_stage)
      ) %>%
      # Convertir a factor las siguientes variables
      mutate(
        across(c(gender, venous_invasion, lymphatic_invasion,
                 stage_event_pathologic_stage, cancer_type), as.factor)
      ) %>%
      # Crear categorías de edades
      mutate(
        age_category = cut(age_at_initial_pathologic_diagnosis,
                           breaks = c(-Inf,59,69,79,Inf),
                           labels = c("Personas_60-","Personas_60","Personas_70","Personas_80+"),
                           include.lowest = FALSE)
      ) %>%
      # Calcular IMC
      mutate(
        imc = calcular_imc(weight, height)
      ) %>%
      # Crear categorías IMC
      mutate(
        imc_category = cut(imc,
                           breaks = c(-Inf, 18.5, 25, 30, Inf),
                           labels = c("Bajo_Peso", "Normal", "Sobrepeso", "Obesidad"))
      )
  })
  
  kaplan_mier_cat <- reactive({
    # Llamar objetos reactivos 
    req(df_gastro(), input$w) 
    
    df <- df_gastro()
     
    # Crear objeto de supervivencia
    survival_obj <- Surv(time = df$os_time, event=df$os_event)
    formula_km   <- as.formula(paste("survival_obj ~", input$w))
    
    # Ajustar el model Kaplan-Meier por categoría seleccionada
    survfit(formula_km, type="kaplan-meier", data=df)
  })
  
  
  # Actualización dinámica de selectInputs según columnas disponibles
  observeEvent(df_gastro(), {

    # Seleccionar variables para el gráfico Supervivencia según variables clínicas
    # Eje x
    varclinic_x <- c("cancer_type", "stage_event_pathologic_stage", "venous_invasion", 
                           "lymphatic_invasion", "perineural_invasion_present")
    # Eje y
    varclinic_y <- c("os_time", "age_at_initial_pathologic_diagnosis")
    
    # Seleccionar variables para el gráfico de histogramas
    var_hist <- c("age_category", "imc_category")
    
    # Seleccionar variables para el gráfico de supervivencia Kaplan Meier
    var_kaplan <- c("age_category", "gender")

    updateSelectInput(
      session,
      "x",
      choices = varclinic_x,
      selected = varclinic_x[1]
    )

    updateSelectInput(
      session,
      "y",
      choices = varclinic_y,
      selected = varclinic_y[1]
    )

    updateSelectInput(
      session,
      "z",
      choices = var_hist,
      selected = var_hist[1])
    
    updateSelectInput(
      session,
      "w",
      choices = var_kaplan,
      selected = var_kaplan[1])
  })
  
  # Resumen estadístico
  output$summary <- renderPrint({
    req(df_gastro())
    summary(df_gastro())
  })
  
  # Histograma
  output$histogram <- renderPlot({
    # Llamar objetos reactivos
    req(df_gastro(), input$z)
    
    # Filtrar NAs
    df <- df_gastro()
    df_filtrado <- df[!is.na(df[[input$z]]), ]
    
    # Graficar Histograma
    ggplot(data=df_filtrado, 
      aes_string(x = input$z, fill = "gender"))+  
      geom_bar(position = "dodge", alpha = 0.7)+
      labs(title = paste("Distribución de", input$z), x = input$z, y = "Frecuencia") +
      theme_minimal()+
      #scale_fill_brewer(palette="PRGn")+
      theme(plot.title = element_text(size=12, color='Darkblue', face='bold', hjust = 0.5),
            axis.text.x = element_blank())
  })
  
  # Gráfico de Kaplan Mier
  output$probability <- renderPlot({
    # Llamar objetos reactivos
    req(kaplan_mier_cat(), input$w)
    
    fit <- kaplan_mier_cat()
    
    # Graficar Supervivencia de Kaplan Mier
    plot(fit, ylab = "Probabilidad", xlab='Tiempo (días)',
         mark.time = TRUE, col = hue_pal()(length(fit$strata)),
         main = paste("Supervivencia según", input$w),
         col.main = "darkblue",   
         cex.main = 1.0,          
         font.main = 2)           
    
    legend("bottomleft", inset = c(0, 0),
           legend = names(fit$strata),
           fill = hue_pal()(length(fit$strata)),
           bty = "n",
           cex = 0.5)
  })
  
  # Supervivencia según variables clínicas
  output$survival <- renderPlot({
    # Llamar objetos reactivos
    req(df_gastro(), input$x, input$y)
    
    # Filtrar NAs
    df <- df_gastro()
    df_filtrado <- df[!is.na(df[[input$x]]), ] 
    
    # Graficar Boxplot
    ggplot(data=df_filtrado, 
      aes_string(x = input$x, y = input$y, fill = input$x))+  
      geom_boxplot(col = 'black')+
      labs(title = paste("Supervivencia en función de", input$x), x = input$x, y = input$y) +
      theme_minimal()+
      #scale_fill_brewer(palette="PRGn")+
      theme(plot.title = element_text(size=12, color='Darkblue', face='bold', hjust = 0.5),
            axis.text.x = element_blank())
  })
}

# Create the Shiny app object --------------------------------------------------
shinyApp(ui = ui, server = server)
