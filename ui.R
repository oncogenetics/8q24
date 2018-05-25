# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# User interface file for shiny

tweaks <- 
  list(
    tags$head(tags$style(HTML("
                              .multicol {
                              -webkit-column-count: 2; /* Chrome, Safari, Opera */
                              -moz-column-count: 2;    /* Firefox */
                              column-count: 2;
                              -moz-column-fill: auto;
                              -column-fill: auto;
                              }
                              "))),
    #push it down 70px to avoid going under navbar
    tags$style(type = "text/css", "body {padding-top: 70px;}"),
    #hide red error messages
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }")
    )

# Define UI ---------------------------------------------------------------
shinyUI(
  navbarPage(
    # Application title
    id = "navBarPageID",
    #title = "8q24 v0.1",
    title = div(h4("8q24 v0.1",
                   style = "margin-top: 0px;"),
                img(src = "icr_logo_white_on_black.PNG", height = "70px",
                    style = "position: relative; top: -60px; right: -800px;")),
    windowTitle = "8q24",
    fluid = FALSE,
    position = "fixed-top",
    inverse = TRUE,
    tabPanel("Plot",
             sidebarPanel(
               tweaks,
               # Type of hits ------------------------------------------------------------
               checkboxGroupInput("methods", "Methods:",
                                  #choices = sort(unique(hitsType$hitType)),
                                  choices = c("final12",
                                              "known16",
                                              "known26",
                                              "credSet38",
                                              "credSet74"),
                                  selected = c("final12", "known16")),
               sliderInput("xRangeZoom", h5("Zoom: start end position"),
                           min = 127330000, max = 129050000,
                           value = c(127600000, 129000000),
                           step = 20000),
               
               hr(),
               #h5("Untick to show other tag SNPs."),
               sliderInput("filterMinLD", h5("LD filter"),
                           min = 0, max = 0.9, value = 0.6, step = 0.05),
               checkboxInput("hitsOnly", label = "Hits only", value = TRUE),
               hr(),
               sliderInput("filterMinBF", h5("BF filter"),
                           min = 0, max = 100, value = 0),
               hr(),
               # Network settings --------------------------------------------------------
               conditionalPanel("input.plotTypePanelSelected=='Network'",
                                selectInput("networkLayout", label = h5("Network layout"),
                                            choices = list("nicely" = "layout_nicely",
                                                           "circle" = "layout_in_circle",
                                                           "star" = "layout_as_star",
                                                           "components" = "layout_components",
                                                           "grid" = "layout_on_grid",
                                                           "sphere" = "layout_on_sphere",
                                                           "random" = "layout_randomly",
                                                           "dh" = "layout_with_dh",
                                                           "drl" = "layout_with_drl",
                                                           "fr" = "layout_with_fr",
                                                           "gem" = "layout_with_gem",
                                                           "graphopt" = "layout_with_graphopt",
                                                           "kk" = "layout_with_kk",
                                                           "lgl" = "layout_with_lgl",
                                                           "mds" = "layout_with_mds",
                                                           "sugiyama" = "layout_with_sugiyama"), 
                                            selected = "layout_nicely"),
                                checkboxInput("nodeSizeBF", label = "BF as node size", value = TRUE)
                                #h5("Untick to have all nodes same size."),
                                ),
               hr(),
               uiOutput("ui_hits"),
               
               width = 3 # width of left sidebar
             ), # END sidebarPanel(
             mainPanel(
               tabsetPanel(
                 id = "plotTypePanelSelected",
                 # Debugging --------------------------------------------------
                 tabPanel("Input",
                          h4("Input"),
                          hr(),
                          h4("hitsOnly"),
                          verbatimTextOutput("I_hitsOnly"),
                          h4("LD filter"),
                          verbatimTextOutput("I_filterLD"),
                          h4("BF filter"),
                          verbatimTextOutput("I_filterBF"),
                          h4("xRangeZoomStart"),
                          verbatimTextOutput("I_xRangeZoomStart"),
                          h4("xRangeZoomEnd"),
                          verbatimTextOutput("I_xRangeZoomEnd"),
                          h4("Selected hits"),
                          verbatimTextOutput("I_hitsSelected"),
                          h4("methods"),
                          verbatimTextOutput("I_methods"),
                          h4("networkLayout"),
                          verbatimTextOutput("I_networkLayout")
                 ),
                 tabPanel("Manhattan",
                          h4("Manhattan"),
                          hr(),
                          plotOutput("PlotManhattan",
                                     click = "plot_click",
                                     dblclick = "plot_dblclick",
                                     hover = "plot_hover",
                                     brush = "plot_brush"),
                          #hr(),
                          #h4("Brush input:"),
                          #verbatimTextOutput("PlotManhattanText"),
                          hr(),
                          h4("Manhattan brush zoom:"),
                          hr(),
                          # zoom with brush
                          plotOutput("PlotManhattanZoom")
                 ),
                 tabPanel("Network",
                          h4("Network"),
                          hr(),
                          h5("Choose layout types: 'nicely', 'circle', 'mds' works best;"),
                          h5("'components' works best to show clusters with tags."),
                          h5("Hover over nodes will display more info about SNPs, and lines(edges) shows LD=R2 value."),
                          h5("To hightlight hits by methods use below dropdown menu:"),
                          visNetworkOutput("networkHits", width = "100%", height = "1000px")
                 ),
                 tabPanel("Arc",
                          h4("Arc"),
                          hr(),
                          plotOutput("PlotArcLD1", width = 1000, height = 500)
                 ),
                 tabPanel("LDmatrix",
                          h4("LD matrix"),
                          hr(),
                          plotOutput("PlotLDmatrix", width = 800, height = 800)
                 ),
                 tabPanel("LD hits vs credSet",
                          h4("LD hits vs credSet"),
                          hr(),
                          DT::dataTableOutput("LD_dt_12_vs_credSet")
                 )
                 
                 
               )
             )
    ), #END tabPanel("Plot",
    tabPanel("Data",
             h4("Data"),
             hr(),
             tabsetPanel(
               tabPanel("Nodes",
                        dataTableOutput("dataNodes")),
               tabPanel("Links",
                        dataTableOutput("dataLinks")),
               tabPanel("LD Subset",
                        dataTableOutput("dataLDSubset")),
               tabPanel("LD Arc",
                        dataTableOutput("dataArcLD")),
               tabPanel("MAP",
                        dataTableOutput("dataMAP")),
               tabPanel("Meta",
                        dataTableOutput("dataMeta")),
               #options
               type = "pills"
             )
    ),
    tabPanel("About")
  ) # END navbarPage(
) # END shinyUI(

