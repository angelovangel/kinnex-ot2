library(shiny)
library(shinydashboard)
library(shinyjs)
#library(shinyalert)
library(plater)
library(tibble)
library(stringr)
library(reactable)
library(dplyr)
library(rmarkdown)
library(curl)
library(plater)

# for use on rack and aluminium block 
wells24 <- lapply(1:6, function(x) {str_c(LETTERS[1:4], x)}) %>% unlist()
wells96 <- lapply(1:12, function(x) {str_c(LETTERS[1:8], x)}) %>% unlist()


tab1 <-  fluidPage(
  fluidRow(
    box(width = 3, status = "warning", solidHeader = FALSE, height = 420,
      collapsible = F,
      column(12,
        selectizeInput('protocol', 'Kinnex protocol', 
                       choices = c('PacBio 16S' = '16s', 'PacBio full-length RNA' = 'flrna', 'PacBio single-cell RNA' = 'scrna'), 
                       selected = 'PacBio 16S', multiple = FALSE)),
      column(6, 
        selectizeInput('plex', 'Kinnex plex', choices = c(8,12,16), selected = 12, multiple = FALSE)),
      column(6, 
        selectizeInput('nsamples', 'Number of samples', choices = c(1:6), selected = 1, multiple = FALSE)),
      column(12, 
        selectizeInput('ncycles', 'PCR cycles', choices = c(9:12), selected = 9, multiple = FALSE)),
      column(6, 
        numericInput('MMvol', 'Mastermix vol', value = 22.5, step = 0.1, min = 10, max = 25)),
      column(6, 
        numericInput('primervol', 'Primer premix vol', value = 2.5, step = 0.1, min = 1, max = 5)),
      column(12, 
        downloadButton('download_script', 'OT2 script', style = 'margin-top:20px'),
        actionButton('deck', 'Deck', style = 'margin-top:20px'),
        uiOutput('show_protocol', inline = T)
      )
    ),
    box(width = 9, status = "warning", solidHeader = FALSE, height = 420,
      title = htmlOutput('htmlout'), collapsible = F,
      fluidRow(
        column(width = 4,
               tags$p('24 tuberack eppendorf 1.5ml'),
               reactableOutput('rack')
               ),
        column(width = 8,
               tags$p('PCR plate on ODTC'),
               reactableOutput('reaction_plate')
               )
        )
      )
  ),
  fluidRow(
    box(width = 12, color = "NA", solidHeader = FALSE, title = 'Instructions', collapsible = T, collapsed = T,
        tags$div("For each sample, prepare Mastermix and place in positions marked s1, s2, etc. 
               Place empty 1.5ml tubes in positions p1, p2, etc. The pooled PCR products for 's1' will be placed in 'p1'."
        ),
        tags$div('Mastermix preparation:'),
        reactableOutput('mastermix', width = "25%")
    )
  )
)
  

tab2 <- fluidRow(
  box(width = 12, status = "warning", solidHeader = FALSE, title = "PacBio Kinnex PCR Opentrons script preview", collapsible = F,
      verbatimTextOutput('protocol_preview')
  )
)


ui <- dashboardPage(skin = 'yellow',
                    #useShinyalert(),
                    
                    header = dashboardHeader(title = "PacBio Kinnex PCR for Opentrons", titleWidth = 400),
                    sidebar = dashboardSidebar(disable = T),
                    
                    body = dashboardBody(
                      useShinyjs(),
                      tabsetPanel(
                        tabPanel(title = "Protocol setup", icon = icon("table"), tab1),
                        tabPanel(title = "Opentrons script preview", icon = icon('code'), tab2)
                        #tabPanel(title = 'Instructions', icon = icon('list'), tab3)
                      )
                    )
)


# server #
server = function(input, output, session) {
  
  ### read template
  protocol_url <- "https://raw.githubusercontent.com/angelovangel/opentrons/main/protocols/08-kinnex-pcr.py"
  
  if (curl::has_internet()) {
    con <- url(protocol_url)
    protocol_template <- readLines(con, warn = F)
    close(con)
  } else {
    protocol_template <- readLines('08-kinnex-pcr.py', warn = F)
  }
  
  
  ### REACTIVES
  
  rack_df <- reactive({
    data.frame(
      wells = wells24, 
      samples = 
        c(paste0('s', seq(1,input$nsamples)),
          rep('.', 16-as.numeric(input$nsamples)),
          paste0('p', seq(1,input$nsamples)),
          rep('.', 8-as.numeric(input$nsamples))
          ),
      #palette.colors(4, palette = 'Set3', alpha = 0.3)
      color = c(
        palette.colors(input$nsamples, palette = 'Set3', alpha = 0.4),
        rep(NA, 16-as.numeric(input$nsamples)),
        palette.colors(input$nsamples, palette = 'Set3', alpha = 0.4),
        rep(NA, 8-as.numeric(input$nsamples))
      )
    )
  })
  
  rack <- reactive({
    plater::view_plate(
      data = rack_df(),
      well_ids_column = 'wells', 
      columns_to_display = 'samples', 24
    )
  })
  
  
  plate_df <- reactive({
    data.frame(
      wells = wells96, 
      samples = make_s(rack_df()$samples[1:6], plex = as.numeric(input$plex)),
      color = rep('', 96))
  })
  
  reaction_plate <- reactive({
    plater::view_plate(
      data = plate_df(),
      well_ids_column = 'wells', 
      columns_to_display = 'samples', 96
    )
  })
  
  myprotocol <- reactive({
    
    MMwells <- paste0("['", str_flatten(wells24[1:input$nsamples], collapse = "','"), "']")
    poolwells <- paste0("['", str_flatten(wells24[17:(17+as.numeric(input$nsamples)-1)], collapse = "','"), "']") # A5 is wells24[17]
    
    str_replace(
      string = protocol_template, pattern = 'ncycles = .*', replacement = paste0('ncycles = ', input$ncycles)
      ) %>% 
      str_replace(
        pattern = 'primervol = .*', replacement = paste0('primervol = ', input$primervol)
        ) %>%
      str_replace(
        pattern = 'MMvol = .*', replacement = paste0('MMvol = ', input$MMvol)
      ) %>%
      str_replace(
        pattern = 'MMwells = .*', replacement = paste0('MMwells = ', MMwells)
      ) %>%
      str_replace(
        pattern = 'poolwells = .*', replacement = paste0('poolwells = ', poolwells)
      ) %>%
      str_replace(
        pattern = 'plex = .*', replacement = paste0('plex = ', input$plex)
      )
  })
  
  ### OBSERVERS
  observeEvent(input$deck, {
    showModal(
      modalDialog(title = 'Opentrons deck preview',
                  HTML('<img src="deck.png">'),
                  size = 'l', easyClose = T, 
      )
    )
  })
  
  observeEvent(input$protocol, {
    updateNumericInput(
      session = session, 
      'MMvol', 
      value = case_when(
       input$protocol == 'flrna' ~ 20,
       input$protocol == 'scrna' ~ 15,
       input$protocol == '16s' ~ 22.5,
      )
    )
  })
  
  
  ## OUTPUTS
  output$mastermix <- renderReactable({
    pcrvol <- input$MMvol + input$primervol
    df <- tibble(
      Component = c('Water', 'Kinnex PCR mix(2x)', 'Template'),
      Volume = c(pcrvol*5, pcrvol*8.5, pcrvol*1.8)
    )
    
    reactable(
        df,
        columns = list(
          Component = colDef(footer = 'Total'),
          Volume = colDef(format = colFormat(digits = 1), footer = function(values) sprintf("%.1f", sum(values)))
        ),
        defaultColDef = colDef(footerStyle = list(fontWeight = "bold"))
      )
  })
    
  output$rack <- renderReactable({
    reactable(
      rack()$samples,
      highlight = T, wrap = F, 
      bordered = T, compact = T, fullWidth = T, sortable = F,
      defaultColDef = 
        colDef(minWidth = 40,html = TRUE,
               headerStyle = list(background = "#f7f7f8", fontSize = '80%'), 
               style = function(value) {
                 color <- rack_df()$color[ match(value, rack_df()$sample) ]
                 list(background = color)
               }
        )
    )
  })
  
  output$reaction_plate <- renderReactable({
    reactable(
      reaction_plate()$samples,
      highlight = T, wrap = F, 
      bordered = T, compact = T, fullWidth = T, sortable = F,
      defaultColDef = 
        colDef(minWidth = 40,html = TRUE, headerStyle = list(background = "#f7f7f8", fontSize = '80%'),
               style = function(value) {
                 color <- rack_df()$color[ match(value, rack_df()$sample) ]
                 list(background = color)
               }
        )
    )
  })
  
  ## Decide which protocol to show
  
  output$show_protocol <- renderUI({
    actionButton(
      'showpdf', 'Protocol', 
      style = 'margin-top:20px', 
      onclick = paste0("window.open('", input$protocol, ".pdf')")
      )
  })
  
  output$htmlout <- renderUI({
   paste0("Preview for tuberack and PCR plate")
  })
  
  output$protocol_preview <- renderPrint({
    write(myprotocol(), file = "")
  })
  
  ### DOWNLOADS
  output$download_script <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), '-kinnex-pcr.py')
    },
    content = function(con) {
      # at download time, replace name so that it appears on the Opentrons app
      replacement <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), '-kinnex-pcr.py')
      write(myprotocol() %>%
              str_replace(pattern = "08-kinnex-pcr.py", 
                          replacement = replacement), 
            con)
    }
  )
}


shinyApp(ui, server)