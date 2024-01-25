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

source('global.R')

# for use on rack and aluminium block 
wells24 <- lapply(1:6, function(x) {str_c(LETTERS[1:4], x)}) %>% unlist()
wells96 <- lapply(1:12, function(x) {str_c(LETTERS[1:8], x)}) %>% unlist()


tab1 <-  fluidPage(
  fluidRow(
    box(width = 3, status = "warning", solidHeader = FALSE, height = 550,
      collapsible = F,
      column(12,
             selectizeInput('left_pipet', 'Left pipette (1-channel)', 
                            choices = c('P20 Gen2' = 'p20_single_gen2','P300 Gen2' = 'p300_single_gen2'), 
                            selected = 'p20_single_gen2')),
      column(12,
        selectizeInput('protocol', 'Kinnex protocol', 
                       choices = c('PacBio 16S' = '16s', 'PacBio full-length RNA' = 'flrna', 'PacBio single-cell RNA' = 'scrna'), 
                       selected = 'PacBio 16S', multiple = FALSE)),
      column(12,
        selectizeInput('pipet', 'Left pipette', 
                       choices = c('p20_single_gen2', 'p300_single_gen2'), 
                       selected = 'p20_single_gen2')),
      column(6, 
        selectizeInput('plex', 'Kinnex plex', choices = c(8,10,12,14,16), selected = 12, multiple = FALSE)),
      column(6, 
        selectizeInput('nsamples', 'Number of samples', choices = c(1:6), selected = 1, multiple = FALSE)),
      column(12, 
        selectizeInput('ncycles', 'PCR cycles', choices = c(9:12), selected = 9, multiple = FALSE)),
      column(6, 
        numericInput('MMvol', 'Mastermix vol', value = 22.5, step = 0.1, min = 10, max = 25)),
      column(6, 
        numericInput('primervol', 'Primer premix vol', value = 2.5, step = 0.1, min = 1, max = 5)),
      column(12, 
        downloadButton('download_script', 'OT2 script', style = 'margin-top:5px'),
        #actionButton('deck', 'Deck', style = 'margin-top:5px'),
        uiOutput('show_protocol', inline = T)
      )
    ),
    box(width = 9, status = "warning", solidHeader = FALSE, height = 550,
      title = htmlOutput('htmlout'), collapsible = F,
      fluidRow(
        column(width = 4,
               tags$p('24 tuberack eppendorf 1.5ml - samples and pools'),
               reactableOutput('rack'),
               tags$hr(),
               tags$p('24 aluminium block 0.5ml - Kinnex primers'),
               reactableOutput('block')
               ),
        column(width = 8,
               tags$p('PCR plate on ODTC'),
               reactableOutput('reaction_plate')
               )
        )
      )
  ),
  fluidRow(
    box(width = 12, color = "NA", solidHeader = FALSE, title = 'Instructions', collapsible = T, collapsed = F,
        tags$div("For each sample, prepare Mastermix and place in positions marked s1, s2, etc. 
               Place empty 1.5ml tubes in positions p1, p2, etc. Place Kinnex primer mixes in aluminuium block as shown. The pooled PCR products for 's1' will be placed in 'p1'."
        ),
        tags$div('Mastermix preparation (for 1 sample):'),
        reactableOutput('mastermix', width = "25%")
    )
  )
)
  

tab2 <- fluidRow(
  box(width = 12, status = "warning", solidHeader = FALSE, title = "PacBio Kinnex PCR Opentrons script preview", collapsible = F,
    verbatimTextOutput('protocol_preview')
  )
)

tab3 <- fluidRow(
  box(
    width = 12, status = "warning", solidHeader = FALSE, 
    title = tags$div(actionButton('simulate', 'Simulate run'), actionButton('clear', 'Clear output')), 
    collapsible = F,
    #actionButton('simulate', 'Simulate run'),
    verbatimTextOutput('stdout')   
  )
)

tab4 <- fluidRow(
  box(width =12, status = 'warning', solidHeader = FALSE, title = 'Deck view', collapsible = F,
      htmlOutput('deck')
      )
)


ui <- dashboardPage(skin = 'black',
                    #useShinyalert(),
                    
                    header = dashboardHeader(title = "PacBio Kinnex PCR for Opentrons", titleWidth = 400),
                    sidebar = dashboardSidebar(disable = T),
                    
                    body = dashboardBody(
                      useShinyjs(),
                      tabsetPanel(
                        tabPanel(title = "Protocol setup", icon = icon("vials"), tab1),
                        tabPanel(title = "Opentrons script", icon = icon('code'), tab2),
                        tabPanel(title = 'Opentrons simulate', icon = icon('code'), tab3),
                        tabPanel(title = 'Deck view', icon = icon('border-none'), tab4)
                        #tabPanel(title = 'Instructions', icon = icon('list'), tab3)
                      )
                    )
)


# server #
server = function(input, output, session) {
  
  # add opentrons_simulate path
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, Sys.getenv('OPENTRONS_PATH'), sep = ":"))
  
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
  mmix_react <- reactiveValues(
    water = '132-X ul',
    kinnexmix = '165 ul',
    template = 'X ul (35 ng)',
    total = '297 ul'
  )
  
  rack_df <- reactive({
    data.frame(
      wells = wells24, 
      samples = 
        c(paste0('<b>s', seq(1,input$nsamples)),
          rep('.', 16-as.numeric(input$nsamples)),
          paste0('<b>p', seq(1,input$nsamples)),
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
      columns_to_display = 'samples', plate_size = 24
    )
  })
  
  block_df <- reactive({
    kinnexprimers <- 
      paste0(
        '<b>',
        c(
          LETTERS[ 1:input$plex - 1 ],
          paste0(LETTERS[as.numeric(input$plex)], 'Q')
        )
      )
    data.frame(
      wells = wells24,
      primers = 
        c(kinnexprimers,
          rep('.', 24 - length(kinnexprimers))
        )
    )
  })
  
  block <- reactive({
    plater::view_plate(
      data = block_df(),
      well_ids_column = 'wells',
      columns_to_display = 'primers', plate_size = 24
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
      columns_to_display = 'samples', plate_size = 96
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
      ) %>%
      str_replace(
        pattern = 'left_pipette = .*', replacement = paste0("left_pipette = ", "'", input$left_pipet, "'")
      )
  })
  
  ### OBSERVERS
  observeEvent(input$clear, {
    shinyjs::html(id = "stdout", "")
  })
  
  observeEvent(input$simulate, {
    # clear stdout
    shinyjs::html(id = "stdout", "")
    # check if opentrons_simulate is in path
    if (system2('which', args = 'opentrons_simulate') == 1) {
      shinyjs::html(id = 'stdout', "opentrons_simulate executable not found. Set the OPENTRONS_PATH variable to the opentrons path.")
      return()
    }
    # change button when working
    shinyjs::disable(id = 'simulate')
    shinyjs::html(id = 'simulate', "Working...")
    tmp <- tempfile('protocol', fileext = '.py')
    write(myprotocol(), file = tmp)
    
    withCallingHandlers({
      processx::run(
        'opentrons_simulate', 
        args = c('-e', tmp),
        stderr_to_stdout = TRUE, 
        error_on_status = FALSE,
        stdout_line_callback = function(line, proc) {message(line)}, 
      )
      shinyjs::enable(id = 'simulate')
      shinyjs::html(id = 'simulate', "Simulate run")
    },
      message = function(m) {
        shinyjs::html(id = "stdout", html = m$message, add = TRUE); 
        runjs("document.getElementById('stdout').scrollTo(0,1e9);") 
        # scroll the page to bottom with each message, 1e9 is just a big number
      }
    )
  })
  
  observeEvent(input$protocol, {
    
    updateSelectizeInput(
      session = session,
      'plex',
      selected = case_when(
        input$protocol == 'flrna' ~ 8,
        input$protocol == 'scrna' ~ 16,
        input$protocol == '16s' ~ 12
      )
    )
    
    mmix_react$water = case_when(
      input$protocol == 'flrna' ~ '88 - X ul',
      input$protocol == 'scrna' ~ '176 - X ul',
      input$protocol == '16s' ~ '132 - X ul'
    )
    mmix_react$kinnexmix = case_when(
      input$protocol == 'flrna' ~ '110 ul',
      input$protocol == 'scrna' ~ '220 ul',
      input$protocol == '16s' ~ '165 ul'
    )
    mmix_react$template = case_when(
      input$protocol == 'flrna' ~ 'X ul (55 ng)',
      input$protocol == 'scrna' ~ 'X ul (55 ng)',
      input$protocol == '16s' ~ 'X ul (35 ng)'
    )
    mmix_react$total = case_when(
      input$protocol == 'flrna' ~ '198 ul',
      input$protocol == 'scrna' ~ '396 ul',
      input$protocol == '16s' ~ '297 ul'
    )
    
  })
  
  
  
  
  ## OUTPUTS
  output$mastermix <- renderReactable({
    pcrvol <- input$MMvol + input$primervol
    df <- tibble(
      Component = c('Water', 'Kinnex PCR mix(2x)', 'Template'),
      Volume = c(mmix_react$water, mmix_react$kinnexmix, mmix_react$template)
    )
    
    reactable(
        df,
        columns = list(
          Component = colDef(footer = 'Total'),
          Volume = colDef(footer = mmix_react$total)
        ),
        defaultColDef = colDef(footerStyle = list(fontWeight = "bold"))
      )
  })
    
  output$rack <- renderReactable({
    reactable(
      rack()$samples,
      highlight = T, wrap = F, 
      bordered = T, compact = F, fullWidth = T, sortable = F,
      defaultColDef = 
        colDef(minWidth = 45,html = TRUE,
               headerStyle = list(background = "#f7f7f8", fontSize = '90%'), 
               style = function(value) {
                 color <- rack_df()$color[ match(value, rack_df()$sample) ]
                 list(background = color)
               }
        )
    )
  })
  
  output$block <- renderReactable({
    reactable(
      block()$primers,
      highlight = T, wrap = F, 
      bordered = T, compact = F, fullWidth = T, sortable = F,
      defaultColDef = 
        colDef(minWidth = 45, html = TRUE,
               headerStyle = list(background = "#f7f7f8", fontSize = '90%'))
    )
  })
  
  output$reaction_plate <- renderReactable({
    reactable(
      reaction_plate()$samples,
      highlight = T, wrap = F, 
      bordered = T, compact = F, fullWidth = T, sortable = F,
      defaultColDef = 
        colDef(minWidth = 40,html = TRUE, headerStyle = list(background = "#f7f7f8", fontSize = '90%'),
               style = function(value) {
                 color <- rack_df()$color[ match(value, rack_df()$sample) ]
                 list(background = color)
               }
        )
    )
  })
  
  output$deck <- renderUI({
    HTML('<img src="deck.png" height="600">')
  })
  
  ## Decide which protocol to show
  
  output$show_protocol <- renderUI({
    actionButton(
      'showpdf', 'Protocol', 
      style = 'margin-top:5px', 
      onclick = paste0("window.open('", input$protocol, ".pdf')")
      )
  })
  
  output$htmlout <- renderUI({
   paste0("Preview for tuberack/primer block and PCR plate")
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