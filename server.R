function(input, output) {
  
  
  # hits ----------------------------------------------------------------------
  hitsSelected <- reactive({ input$hits })
  
  methodsSelected <- reactive({ input$methods })
  
  zoomStart <- reactive({ input$xRangeZoom[1] })
  zoomEnd <- reactive({ input$xRangeZoom[2] })
  
  Xrange <- reactive({ c(zoomStart(), zoomEnd())})
  Yrange <- reactive({c(0, round(max(-log10(metaRegionClean$P)), -1) + 10)})
  
  hitsLDsubset <- reactive({
    hitsLDclean[ (R2 > input$filterMinLD |
                    SNP_A %in% hitsSelected() |
                    SNP_B %in% hitsSelected()) &
                   BP_A >= zoomStart() &
                   BP_A <= zoomEnd() &
                   BP_B >= zoomStart() &
                   BP_B <= zoomEnd(), 
                 ][# filter on BF
                   SNP_A %in% metaRegionClean[ MeanBF > input$filterMinBF, SNP] &
                     SNP_B %in% metaRegionClean[ MeanBF > input$filterMinBF, SNP], ]  
  })
  
  # Dynamic UI --------------------------------------------------------------
  output$ui_hits <- 
    renderUI({
      list(tags$div(align = 'left', 
                    class = 'multicol', 
                    checkboxGroupInput("hits", h5("Hits:"),
                                       choices = sort(
                                         unique(
                                           hitsType[ hitType %in%
                                                       methodsSelected() &
                                                       # filter on BF
                                                       SNP %in% metaRegionClean[ MaxBF > input$filterMinBF, SNP],
                                                     SNP])),
                                       selected = unique(hitsType[ hitType %in%
                                                                     methodsSelected() &
                                                                     BP >= zoomStart() &
                                                                     BP <= zoomEnd(),
                                                                   SNP]))
      ))
    })
  # Network: nodes and links --------------------------------------------------
  # 1. Nodes ------------------------------------------------------------------
  nodes <- reactive({
    nodeHits <- hitsType[ hitType %in% input$methods, 
                          ][ order(hitType), .(SNPid, id = SNP, label = SNP,
                                               hitType = hitType)]
    nodeHits[, group := paste(hitType, collapse = ","),
             by = .(SNPid, id, label) ]
    
    nodeTags <- unique(
      rbind(
        hitsLDsubset()[, .( SNPid = SNPid_A, id = SNP_A, 
                            label = SNP_A, hitType = NA, group = "tag") ],
        hitsLDsubset()[, .( SNPid = SNPid_B, id = SNP_B,
                            label = SNP_B, hitType = NA, group = "tag") ]
        )[ !(id %in% nodeHits$id), ]
      )[ id %in% metaRegionClean[ MaxBF > input$filterMinBF, SNP], ]
    
    #merge hit nodes and other tag nodes
    if(nrow(nodeTags) > 0) {
      res <- rbind(nodeHits, nodeTags)} else {
        res <- nodeHits}
    
    #keep only hits?
    if(input$hitsOnly){
      res <- res[ id %in% hitsSelected(), ]
    } 
    
    # remove nodes that have no links
    res <- res[ id %in% hitsSelected() |
                  id %in% c(links()$from, links()$to), ]
    
    #add meta info: pvalue, infoscore, allele, maf
    res <- merge(res,
                 metaRegionClean[, .(SNPid, MaxBF, title)],
                 by = "SNPid", all.x = TRUE)
    
    # set the size of nodes
    if(input$nodeSizeBF){
      res[, size := as.numeric(
        as.character(cut(MaxBF,
                         breaks = c(-1, 0, 1,  3,  5, 10, 100, Inf),
                         labels = c(1, 4, 10, 12, 15, 20, 25)))) ]
      res <- res[ order(group), .(id, label, group, size, title), ]
    } else {
      res <- res[ order(group), .(id, label, group, title) ]
      }
    #return
    unique(res)
  })
  
  # 2. Links ---------------------------------------------------------------------
  links <- reactive({
    if(input$hitsOnly){
      res <- hitsLDsubset()[, .(from = SNP_A, to = SNP_B, value = R2)
                            ][ from != to, 
                               ][from %in% hitsSelected() & 
                                   to %in% hitsSelected(), ]
      } else {
        res <- hitsLDsubset()[, .(from = SNP_A, to = SNP_B, value = R2) 
                              ][from != to &
                                  from %in% hitsSelected(), ]
      }
    res[, color := "lightblue" ]
    res[, title := paste0("R2:", round(value, 2)) ]
    
    #return
    res[ value >= input$filterMinLD, ]
  })
  
  # Arc plots -----------------------------------------------------------------
  # 1. LD data ----------------------------------------------------------------
  LD <- reactive({
    # data for plotLDarc(), merge to MAP twice to get BP, SNP_A is hits
    # c("BP_A","SNP_A","BP_B","SNP_B","R2")
    merge(
      merge(hitsLD[ R2 > input$filterMinLD, ],
            MAP[ , .(SNPid, rsid, BP)], by.x = "SNP", by.y = "SNPid"),
      MAP[ , .(SNPid, rsid, BP)], by.x = "SNP_hit", by.y = "SNPid"
    )[, .(BP_A = BP.y, SNP_A = SNP_hit, BP_B = BP.x, SNP_B = SNP, R2)]
  })
  
  # Input -------------------------------------------------------------------
  output$I_hitsOnly <- renderPrint({ input$hitsOnly })
  output$I_xRangeZoomStart <- renderPrint({ zoomStart() })
  output$I_xRangeZoomEnd <- renderPrint({ zoomEnd() })
  output$I_hitsSelected <- renderPrint({ hitsSelected() })
  output$I_filterLD <- renderPrint({ input$filterMinLD })
  output$I_filterBF <- renderPrint({ input$filterMinBF })
  output$I_methods <- renderPrint({ input$methods })
  output$I_networkLayout <- renderPrint({ input$networkLayout })
  
  # Data ----------------------------------------------------------------------
  output$dataMAP <- renderDT({datatable(MAP)})
  
  output$dataMeta <- renderDT({
    x <- setdiff(colnames(metaRegionClean), 
                 c("SNPid","MarkerName","V1", "title", "rn"))
    datatable( metaRegionClean[ order(BP), ..x ] )
    })
  
  output$dataLDSubset <- renderDT({ datatable(hitsLDsubset()) })

  output$dataArcLD <- renderDataTable({ LD() })
  output$dataNodes <- renderDataTable({ nodes() })
  output$dataLinks <- renderDataTable({ links() })
  
  # LD Datatable, hits vs credSet -------------------------------------------
  output$LD_dt_12_vs_credSet <- renderDT({
    x <- copy(hitsLDclean)
    x <- x[ (SNP_B %in% hitsType[ hitType == "final12", SNP]), 
       ][, .(SNP = SNP_A, SNP_B, R2 = round(R2, 2)) 
         ][, dcast(.SD, SNP ~ SNP_B, fill = 0, value.var = "R2" )] 

    colNamesOrder <- unique(
        hitsLDclean[ SNP_A %in% hitsType[ hitType == "final12", SNP], 
                     ][order(BP_A) , SNP_A ])
    rowNamesOrder <- unique(
      hitsLDclean[ SNP_B %in% x$SNP, .(SNP_B, BP_B)
                   ][order(BP_B), SNP_B])
    setDF(x)
    rownames(x) <- x$SNP
    x <- x[ match(rowNamesOrder, x$SNP), colNamesOrder ]
    
    datatable(x, 
              options = list(pageLength = 100),
              escape = FALSE) %>%
      formatStyle(
        columns = 1:12,
        backgroundColor = styleInterval(
          cuts = seq(0, 0.8, 0.2),
          values = c("white","#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000")))
  })

  # Plot: network -------------------------------------------------------------
  output$networkHits <- renderVisNetwork({
    visNetwork(nodes(), links()) %>%
      visOptions(highlightNearest = TRUE,
                 selectedBy = "group") %>%
      visLegend() %>%
      visIgraphLayout(layout = input$networkLayout, randomSeed = 12)
  })
  
  # Plot: arc -----------------------------------------------------------------
  plotObjArcLD <- reactive({
    
    statNameSelected <- c(
      input$methods[1],
      ifelse(is.na(input$methods[2]), input$methods[1], input$methods[2]))
    
    plotDat <- hitsLDclean[ (BP_A >= Xrange()[1] & BP_A <= Xrange()[2] &
                               BP_B >= Xrange()[1] & BP_B <= Xrange()[2]), ]

    plotLDarc(plotDat,
              statNames = statNameSelected,
              upper = hitsType[hitType == statNameSelected[1], SNP],
              lower = hitsType[hitType == statNameSelected[2], SNP],
              minR2 = input$filterMinLD,
              hitsOnly = input$hitsOnly,
              title = paste(statNameSelected, collapse = " vs "),
              xStart = Xrange()[1],
              xEnd =   Xrange()[2]
    )
  })
  output$PlotArcLD1 <- renderPlot({
    print(plotObjArcLD())})
  
  # Plot: manhattan ---------------------------------------------------------
  output$PlotManhattan <- renderPlot({
    x <- copy(metaRegionClean)
    x <- x[ BP >= zoomStart() & BP <= zoomEnd(), ]
    x[, Pos := BP ]
    x[, PLog10 := -log10(`P-value`) ]
    x[, col := ifelse(SNP %in% hitsSelected(), "red", "grey") ]
    x[, label := ifelse(SNP %in% hitsSelected(), SNP, NA_character_) ]

    ggplot(x, aes(Pos, PLog10, col = col, label = label)) +
      geom_point(size = 4, alpha = 0.5) +
      geom_text_repel(col = "black", na.rm = TRUE) +
      scale_colour_identity() + 
      coord_cartesian(xlim = Xrange(), ylim = Yrange()) +
      theme_minimal()
  })

  output$PlotManhattanText <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    paste0(
      "click: ", xy_str(input$plot_click),
      "dblclick: ", xy_str(input$plot_dblclick),
      "hover: ", xy_str(input$plot_hover),
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  
  
  output$PlotManhattanZoom <- renderPlot({
    if(is.null(input$plot_brush)){
      xRangeHover <- Xrange()
    } else {
      xRangeHover <- c(input$plot_brush$xmin, input$plot_brush$xmax)  
    }
    
    if(is.null(input$plot_brush)){
      yRangeHover <- Yrange()
    } else {
      yRangeHover <- c(input$plot_brush$ymin, input$plot_brush$ymax)  
    }
    
    x <- copy(metaRegionClean)
    
    x <- x[ BP >=  xRangeHover[1] & BP <= xRangeHover[2], ]
    x[, Pos := BP ]
    x[, PLog10 := -log10(`P-value`) ]
    x[, col := ifelse(SNP %in% hitsSelected(), "red", "grey") ]
    x[, label := ifelse(SNP %in% hitsSelected(), SNP, NA_character_) ]
    x <- x[ PLog10 >= yRangeHover[1] & PLog10 <= yRangeHover[2], ]
    
    ggplot(x, aes(Pos, PLog10, col = col, label = label)) +
      geom_point(size = 4, alpha = 0.5) +
      geom_text_repel(col = "black", na.rm = TRUE) +
      scale_colour_identity() + 
      coord_cartesian(xlim = xRangeHover, ylim = yRangeHover) +
      theme_minimal()
  })
  
  # Plot: LD Matrix ---------------------------------------------------------  
  output$PlotLDmatrix <- renderPlot({
    d <- hitsLDclean[ SNP_A %in% hitsSelected() & SNP_B %in% hitsSelected(), ]
    # get all combo
    x <- data.frame(expand.grid(unique(c(d$SNP_A, d$SNP_B)),
                                unique(c(d$SNP_A, d$SNP_B))))
    colnames(x) <- c("SNP_A", "SNP_B")
    # add R2
    x <- merge(x, d[, .(SNP_A, SNP_B, R2)], by = c("SNP_A", "SNP_B"), all.x = TRUE)
    x[ is.na(x) ] <- 0
    
    # add pos for SNP_A and SNP_B
    res <- x
    d_bp <- merge(res, unique(hitsType[ , .(SNP_A = SNP, BP_A = BP)]), by = "SNP_A")
    d_bp <- merge(d_bp, unique(hitsType[ , .(SNP_B = SNP, BP_B = BP)]), by = "SNP_B")
    
    # order SNPs by pos - as factor levels
    SNPorder <- unique(
      rbind(data.frame(SNP = d_bp$SNP_A, BP = d_bp$BP_A),
            data.frame(SNP = d_bp$SNP_B, BP = d_bp$BP_B)))
    SNPorder <- SNPorder[ order(SNPorder$B), "SNP" ]
    
    res$SNP_A <- factor(res$SNP_A, levels = SNPorder)
    res$SNP_B <- factor(res$SNP_B, levels = SNPorder)
    
    # R2 cut colours
    #res$col <- as.character(cut(res$R2, seq(0, 1, 0.2), labels = grey.colors(5, 0.9, 0)))
    res$col <- as.character(cut(res$R2, seq(0, 1, 0.2),
                                labels = c("#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000")))
    res$col <- ifelse(res$R2 == 0, "white", res$col)
    
    # output plot tile
    ggplot(res, aes(SNP_A, SNP_B, fill = col, col = "grey90")) +
      geom_tile() +
      scale_fill_identity(guide = "legend", name = expression(R^2),
                          breaks = rev(c("white","#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000")),
                          limits = rev(c("white","#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000")),
                          labels = rev(seq(0, 1, 0.2)), drop = FALSE) +
      scale_color_identity() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
            axis.title = element_blank(),
            legend.position = c(0.9, 0.4),
            legend.background = element_rect(colour = "grey90")) 
  })
  
} # END function(input, output) {

