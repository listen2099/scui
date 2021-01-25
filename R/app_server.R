#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import Seurat
#' @import shinythemes
#' @import shinydashboard
#' @import shinycssloaders
#' @import plotly
#' @import ggplot2
#' @import shinyWidgets
#' @import grDevices
#' @import DT
#' @import dplyr
#' @import htmlwidgets
#' @import shinyBS
#' @improt heatmaply
#' @import htmltools
#' @import RColorBrewer
#' @import shinyjs
#' @noRd
app_server <- function( input, output, session ) {
  obj_list <- mod_import_seurat_obj_server("SeuratRdataFile")
  
  values <- reactiveValues(my_gene_list = NULL, # 似乎没有用了
                           seurat_obj = NULL,
                           def_cate_clu = list(), # 自定义的分类
                           def_gene_list = list('Default Feature List' = ''), # 自定义的基因列表
                           category_name = NULL, # 定义的类别名称
                           heatmap_height = 1, # 热图一个基因的高度单位
                           violin_height = 4, # 小提琴一个基因的高度单位
                           dotplot_height = 1,
                           if_show_add_gene_list = F, # 是否展示加基因的框
                           marker_res_save = list(),# maker gene的保存
                           click_res = list(# 点击细胞以后的结果保存
                             click_category_name = 'str',
                             click_cluster_name = 'str',
                             click_cell_num = 0
                           ),
                           new_obj = NULL,
                           new_obj_markers = NULL
                           )
  
  CellBrowserUI_obj1 <- reactive({menuItem("Cell Browser", tabName = "CellBrowser", icon = icon("eye-open", lib = "glyphicon"))})
  CellBrowserUI_obj2 <- reactive({NULL})
  observe({
    if(! is.null(values$seurat_obj) ){
      output$CellBrowserUI <- renderMenu({CellBrowserUI_obj1()})
    }else{
      output$CellBrowserUI <- renderMenu({CellBrowserUI_obj2()})
      }
  })

  FeaturesUI_obj1 <- reactive({menuItem("Features", tabName = "Features" ,icon = icon("list-alt", lib = "glyphicon"))})
  FeaturesUI_obj2 <- reactive({NULL})
  observe({
    if(input$tabs == 'CellBrowser'){
      output$FeaturesUI <- renderMenu({FeaturesUI_obj1()}
      )}else{
        output$FeaturesUI <- renderMenu({FeaturesUI_obj2()})
      }
  })

  HeatmapUI_obj1 <- reactive({menuItem("Heat Map", tabName = "Heatmap", icon = icon("th", lib = "glyphicon"))})
  HeatmapUI_obj2 <- reactive({NULL})
  observe({
    if(input$tabs == 'Features'){output$HeatmapUI <- renderMenu({HeatmapUI_obj1()})}
    else{
      output$HeatmapUI <- renderMenu({HeatmapUI_obj2()})
    }
  })

  Violin_obj1 <- reactive({menuItem("Violin", tabName = "Violin", icon = icon("equalizer", lib = "glyphicon"))})
  Violin_obj2 <- reactive({NULL})
  observe({
    if(input$tabs == 'Heatmap'){output$ViolinUI <- renderMenu({Violin_obj1()})}else{
      output$ViolinUI <- renderMenu({Violin_obj2()})
    }
  })

  DotPlot_obj1 <- reactive({menuItem("Dot Plot", tabName = "DotPlot", icon = icon("braille", lib = "font-awesome"))})
  DotPlot_obj2 <- reactive({NULL})
  observe({
    if(input$tabs == 'Violin'){output$DotPlotUI <- renderMenu({DotPlot_obj1()})}else{
      output$DotPlotUI <- renderMenu({DotPlot_obj2()})
    }
  })

  Re_clustering_obj1 <- reactive({menuItem("Re-clustering", tabName = "Re_clustering", icon = icon("laptop-code", lib = "font-awesome"))})
  Re_clustering_obj2 <- reactive({NULL})
  observe({
    if(length(unlist(values$def_cate_clu)) > 0){
      output$Re_clusteringUI <- renderMenu({Re_clustering_obj1()})
    }else{
      output$Re_clusteringUI <- renderMenu({Re_clustering_obj2()})}
  })
  
  observe({
    values$seurat_obj <- obj_list()$seurat_obj
  })
  
  savevalues <- reactive({obj_list()$savevalues})
  observe({
    if(!is.null(savevalues())){
      values$my_gene_list <- savevalues()$my_gene_list
      values$def_cate_clu <- savevalues()$def_cate_clu
      values$def_gene_list <- savevalues()$def_gene_list
      values$category_name <- savevalues()$category_name
      values$heatmap_height <- savevalues()$heatmap_height
      values$violin_height <- savevalues()$violin_height
      values$dotplot_height <- savevalues()$dotplot_height
      values$if_show_add_gene_list <- savevalues()$if_show_add_gene_list
      values$marker_res_save <- savevalues()$marker_res_save
      values$click_res <- savevalues()$click_res
    }
  })
  
  plotmodel_obj1 <- reactive({
    fluidRow(
      sliderInput("plotheight", "Plot output height (px)",min = 0, max = 2000, value = 500, step = 100),
      sliderInput("plotwidth", "Plot output width (px)",min = 0, max = 2000, value = 500, step = 100),
      radioGroupButtons(inputId = "plotformat",label = "Plot format",choices = c("svg", "png"),
                                          checkIcon = list(yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                           no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))
    )
  })
  plotmodel_obj2 <- reactive({
    fluidRow()
  })
  
  observeEvent(input$tabs,{
    if(input$tabs == 'CellBrowser'){
      output$plotmodel <- renderMenu({
        plotmodel_obj1()
      })
    }else{
      output$plotmodel <- renderMenu({
        plotmodel_obj2()
      })
    }
  })
  

  output$obj_print <- renderPrint({
    print(values$seurat_obj)
  })

  # 按样本分类展示 ***
  cellmap_obj1 <- reactive({
    if(!is.null(input$input_pick_category) && input$input_pick_category == 'Sample'){
      num <- length( unique(values$seurat_obj@meta.data$orig.ident) )
      #c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFE528","#A65628","#F781BF","#999999")
      colours_sample <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- as.character(unique(values$seurat_obj@meta.data$orig.ident))
      myData <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings[,1:2])
      colnames(myData) <- c('D1','D2')
      myData$Sample <- values$seurat_obj@meta.data$orig.ident
      p <- plot_ly()
      pick <- input$input_pick_cluster_for_sample
      for (i in 1:length(pick)) {
        one_cluster <- pick[i]
        temp <- myData[which(values$seurat_obj@meta.data$orig.ident == one_cluster),]
        p <- add_trace(p,data = temp,type = 'scatter', mode = 'markers',
                       size = I(2.5),
                       color = I(colours_sample[one_cluster]),
                       text = rownames(temp),
                       key = paste0(one_cluster,':',rownames(temp)),
                       x = ~D1,y = ~D2,name = I(one_cluster))
      }
      p %>% #hide_legend() %>%
        layout(legend= list(itemsizing='constant'),
               xaxis = list(showgrid=F,zeroline=F,showticklabels=F,title=paste0(toupper('umap'),'_1')),
               yaxis = list(showgrid=F,zeroline=F,showticklabels=F,title=paste0(toupper('umap'),'_2'))) %>%
        onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
        config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth, filename = "sample_cell_cluster"))
    }else{
      return(NULL)
    }
  })

  # 按seurat中已有的分类展示 ***
  cellmap_obj2 <- reactive({
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)){
      num <- length(unique( Idents(values$seurat_obj) ))
      colours_sample<-colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- as.character(levels(Idents(values$seurat_obj)))
      myData <- as.data.frame(cbind(values$seurat_obj@reductions[[input$input_pick_category]]@cell.embeddings[,1:2]))
      colnames(myData) <- c('D1','D2')
      p <- plot_ly()
      pick <- input$input_pick_cluster
      for (i in 1:length(pick)) {
        one_cluster <- pick[i]
        temp <- myData[which(Idents(values$seurat_obj) == one_cluster),]
        p <- add_trace(p,data = temp,size = I(2.5),
                       type = "scatter", mode = "markers",
                       text = rownames(temp),
                       key = paste0(one_cluster,':',rownames(temp)),
                       color = I(colours_sample[one_cluster]),
                       x = ~D1,y = ~D2,name = I(one_cluster)
        )
      }
      p %>% #hide_legend() %>%
        layout(legend= list(itemsizing='constant'),
               xaxis = list(showgrid=F,zeroline=F,showticklabels=F,title=paste0(toupper(input$input_pick_category),'_1')),
               yaxis = list(showgrid=F,zeroline=F,showticklabels=F,title=paste0(toupper(input$input_pick_category),'_2'))) %>%
        onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
        config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth,filename = "cell_cluster"))
        # config(toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth, filename = "sample_cell_cluster"))
    }else{
      return(NULL)
    }
  })

  # 按某基因的表达量展示
  cellmap_obj3 <- reactive({
    mygenelist <- NULL
    if(!is.null( isolate({input$pick_active_feature_list}) )){
      mygenelist <- isolate({values$def_gene_list[[isolate({input$pick_active_feature_list})]]})
      mygenelist <- mygenelist[which(mygenelist != '')]
    }
    if( !is.null(mygenelist[input$selected_features_0_rows_selected]) ){
      pick_gene <- mygenelist[input$selected_features_0_rows_selected]
      if(length(pick_gene)==1){
        plot.data <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings)
        plot.data[['pick_gene']] <- as.numeric(values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]][pick_gene,])
        barcolor <- colorRamp(c("#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b" ))
        if( all(plot.data[['pick_gene']] == 0) ){
          barcolor <- '#e0e0e0'
        }
        p <- plot_ly(plot.data,x = ~UMAP_1,y = ~UMAP_2,text = rownames(plot.data),
                     key = rownames(plot.data)) %>%
          add_markers(color = ~pick_gene,size = I(2.5),
                      colors = barcolor,
                      type = "scatter", mode = "markers") %>%
          colorbar(title = pick_gene) %>%
          hide_legend() %>%
          layout(xaxis = list(showgrid=F,zeroline=F,showticklabels=F),
                 yaxis = list(showgrid=F,zeroline=F,showticklabels=F)) %>%
          onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth,filename = paste0(pick_gene,"_exp")))
      }else{
        return(NULL)
      }
    }else{
      return(plot_ly())
    }
  })

  # 按用户定义的类别展示 ***
  cellmap_obj4 <- reactive({
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
      num <- length(values$def_cate_clu[[input$input_pick_category]])
      colours_sample<-colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- names(values$def_cate_clu[[input$input_pick_category]])
      
      myData <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings[,1:2])
      pick <- input$input_pick_cluster_for_def
      
      p <- plot_ly()
      for (i in 1:length(pick) ) {
        one_cluster <- pick[i]
        cell_in_cluster <- as.character(values$def_cate_clu[[input$input_pick_category]][[one_cluster]])
        
        temp <- myData[cell_in_cluster,]
        p <- add_trace(p,data = temp,size = I(2.5),
                       type = "scatter", mode = "markers",
                       text = rownames(temp),
                       key = paste0(one_cluster,':',rownames(temp)),
                       color = I(colours_sample[one_cluster]),
                       x = ~UMAP_1,y = ~UMAP_2,
                       name = I(one_cluster)
        )
      }
      p %>%
      layout(legend= list(itemsizing='constant'),
             xaxis = list(showgrid=F,zeroline=F,showticklabels=F),
             yaxis = list(showgrid=F,zeroline=F,showticklabels=F)) %>%
      onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
      config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth,filename = 'my_cell_cluster'))
    }else{
      return(NULL)
    }
  })
  
  # 将细胞split以后展示
  cellmap_obj5 <- reactive({
    if( !is.null(input$ifsplit) && input$ifsplit == T){
      mygenelist <- NULL
      if(!is.null( isolate({input$pick_active_feature_list}) )){
        mygenelist <- isolate({values$def_gene_list[[isolate({input$pick_active_feature_list})]]})
        mygenelist <- mygenelist[which(mygenelist != '')]
      }
      if( !is.null(mygenelist[input$selected_features_0_rows_selected]) ){
        pick_gene <- mygenelist[input$selected_features_0_rows_selected]
        if(length(pick_gene)==1){
          plot.data <- NULL  
          splitby <- input$selectsplitcategory
          Cluster <- NULL
          all_cluster <- NULL
          if( splitby %in% names(values$seurat_obj@reductions) ){
            plot.data <- as.data.frame(values$seurat_obj@reductions[[splitby]]@cell.embeddings)
            colnames(plot.data) <- c('D1','D2')
            Cluster <- as.character(Idents(values$seurat_obj))
            plot.data$Cluster <- Cluster
            plot.data[['pick_gene']] <- as.numeric(values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]][pick_gene,])
            all_cluster <- levels(Idents(values$seurat_obj))
          }else if(splitby == 'Sample'){
            plot.data <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings)
            colnames(plot.data) <- c('D1','D2')
            Cluster <- as.character(values$seurat_obj@meta.data$orig.ident)
            plot.data$Cluster <- Cluster
            plot.data[['pick_gene']] <- as.numeric(values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]][pick_gene,])
            all_cluster <- unique(values$seurat_obj@meta.data$orig.ident)
          }else{
            cluster_list <- values$def_cate_clu[[ splitby ]]
            cells <- as.character(unlist(cluster_list))
            res <- NULL
            for (i in 1:length(cluster_list)) {
              res <- c(res,rep(names(cluster_list[i]),length(cluster_list[[i]])))
            }
            Cluster <- res
            plot.data <- as.data.frame(values$seurat_obj[,cells]@reductions$umap@cell.embeddings)
            colnames(plot.data) <- c('D1','D2')
            plot.data$Cluster <- res
            plot.data[['pick_gene']] <- as.numeric(values$seurat_obj[,cells]@assays[[DefaultAssay(values$seurat_obj)]][pick_gene,])
            all_cluster <- names(values$def_cate_clu[[ splitby ]])
          }
          barcolor <- colorRamp(c("#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b" ))
          if( all(plot.data[['pick_gene']] == 0) ){
            barcolor <- '#e0e0e0'
          }
          anno_gen <- function(title){
            f <- list(#family = "Courier New, monospace",
              size = 14,color = "black")
            a <- list(text = title,font = f,xref = "paper",yref = "paper",
                      yanchor = "bottom",xanchor = "center",align = "center",
                      x = 0.5,y = 1,showarrow = FALSE)
            return(a)
          }
          all_G <- list()
          for (i in 1:length(all_cluster)) {
            temp <- plot.data[which(plot.data$Cluster == all_cluster[i]),]
            max_y <- max(plot.data$D2)
            min_y <- min(plot.data$D2)
            max_x <- max(plot.data$D1)
            min_x <- min(plot.data$D1)
            barcolor <- colorRamp(c("#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b" ))
            if( all(temp[['pick_gene']] == 0) ){
              barcolor <- '#e0e0e0'
            }
            p <- plot_ly(temp,x = ~D1,y = ~D2) %>% 
              add_markers(color = ~pick_gene,
                          size = I(2.5),
                          key = rownames(temp),
                          colors = barcolor,
                          type = "scatter", mode = "markers") %>% 
              hide_legend() %>%
              layout(
                annotations = anno_gen(all_cluster[i]),
                xaxis = list(showgrid = F,zeroline = F,showticklabels = F,
                             showline = TRUE,mirror = "ticks",
                             linecolor = toRGB("black"),linewidth = 1,
                             range = c(min_x, max_x)),
                yaxis = list(showgrid = F,zeroline = F,showticklabels = F,
                             showline = TRUE,mirror = "ticks",
                             linecolor = toRGB("black"),linewidth = 1,
                             range = c(min_y, max_y)))
            if(i == 1){
              p <- p %>% colorbar(title = pick_gene,len = 0.5,limits = c(min(plot.data$pick_gene),max(plot.data$pick_gene)))
            }else{
              p <- p %>% colorbar(title = pick_gene,len = 0.5,limits = c(min(plot.data$pick_gene),max(plot.data$pick_gene))) %>% hide_colorbar()
            }
            all_G[[length(all_G) + 1]] <- p
          }
          
          if( length(all_G) <= 3 ){
            subplot(all_G,nrows = 1) %>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
              config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth,filename=paste0(pick_gene,"_split")))
          }else{
            subplot(all_G,nrows = ceiling(sqrt(length(all_G)))) %>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
              config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat, height = input$plotheight, width = input$plotwidth,filename=paste0(pick_gene,"_split")))
          }
        }else{
          return(NULL)
        }
      }else{
        return(plot_ly())
      }
    }else{
      return(NULL)
    }

  })
  
  # 按seurat中已有的分类展示得到表格
  celltable_obj2 <- reactive({
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)){
      myData <- as.data.frame(cbind(values$seurat_obj@reductions[[input$input_pick_category]]@cell.embeddings[,1:2]))
      myData$orig.ident <- values$seurat_obj@meta.data$orig.ident
      res_table <- c()
      pick <- input$input_pick_cluster
      for (i in 1:length(pick)) {
        one_cluster <- pick[i]
        temp <- myData[which(Idents(values$seurat_obj) == one_cluster),]
        temp$cluster <- one_cluster
        res_table <- rbind(res_table,temp)
      }
      return(res_table)
    }else{
      return(NULL)
    }
  })
  
  # 按用户定义的类别得到表格
  celltable_obj4 <- reactive({
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
      myData <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings[,1:2])
      myData$orig.ident <- values$seurat_obj@meta.data$orig.ident
      pick <- input$input_pick_cluster_for_def
      res_table <- c()
      for (i in 1:length(pick) ) {
        one_cluster <- pick[i]
        cell_in_cluster <- as.character(values$def_cate_clu[[input$input_pick_category]][[one_cluster]])
        temp <- myData[cell_in_cluster,]
        temp$cluster <- one_cluster
        res_table <- rbind(res_table,temp)
      }
      return(res_table)
    }else{
      return(NULL)
    }
  })
  
  saveprojectbutton_obj <- reactive({
    downloadButton(
      outputId = "output_scui_obj",
      label = "Save Project",
      style = "color: white; background-color: #1558A0;border-color: #1558A0; height:50px",
      width = "100%"
    )
  })
  
  output$saveprojectbutton <- renderUI({
    if( !is.null(values$seurat_obj) ){
      return(saveprojectbutton_obj())
    }else{
      return(NULL)
    }
  })
  
  output$output_cell_table_buttom <- renderUI({
    if(
      (!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)) ||
      (!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions))
    ){
      if( length(input$selected_features_0_rows_selected) == 0 ){
        downloadButton(
          outputId = "output_cell_table",
          label = "Output",
          style = "unite",
          width = "100%"
        )      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  })
  
  output$output_cell_table <- downloadHandler(
    filename = function(){
      'scui_cell_meta_table.csv'},
    content = function(con){
      r2 <- celltable_obj2()
      r4 <- celltable_obj4()
      res <- NULL
      if( length(input$selected_features_0_rows_selected) == 1 ){
        if( !is.null(input$ifsplit) && input$ifsplit == T ){}else{}
      }else{
        if(!is.null(input$input_pick_category) && input$input_pick_category == 'Sample'){
        }
        if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)){
          res <- r2
        }
        if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
          res <- r4
        }
      }
      write.csv(res,con,quote = FALSE,row.names = TRUE,col.names = TRUE)
    }
  )
  
  output$output_scui_obj <- downloadHandler(
    filename = function(){
      'scui_project_obj.Rdata'
    },
    content = function(con){
      
      seurat_obj <- values$seurat_obj
      seurat_obj_markers <- obj_list()$seurat_obj_markers
      seurat_obj_name <- values$seurat_obj_name
      
      savevalues <- list()
      savevalues[['my_gene_list']] <- values$my_gene_list
      savevalues[['def_cate_clu']] <- values$def_cate_clu
      savevalues[['def_gene_list']] <- values$def_gene_list
      savevalues[['category_name']] <- values$category_name
      savevalues[['heatmap_height']] <- values$heatmap_height
      savevalues[['violin_height']] <- values$violin_height
      savevalues[['dotplot_height']] <- values$dotplot_height
      savevalues[['if_show_add_gene_list']] <- values$if_show_add_gene_list
      savevalues[['marker_res_save']] <- values$marker_res_save
      savevalues[['click_res']] <- values$click_res

      save(
        list = as.character(c(
          'seurat_obj',
          'seurat_obj_markers',
          'seurat_obj_name',
          'savevalues'
        )),
        file = con
      )

    }
  )
  
  output$monitor <- renderPrint({
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
      # 选择自定义类的情况下
      num <- length(values$def_cate_clu[[input$input_pick_category]])
      colours_sample <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- names(values$def_cate_clu[[input$input_pick_category]])
      
      myData <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings[,1:2])
      pick <- input$input_pick_cluster_for_def
      
      res_myData <- c()

      for (i in 1:length(pick) ) {
        one_cluster <- pick[i]
        cell_in_cluster <- as.character(values$def_cate_clu[[input$input_pick_category]][[one_cluster]])
        temp <- myData[cell_in_cluster,]
        temp$text <- rownames(temp)
        temp$key <- paste0(one_cluster,':',rownames(temp))
        temp$color <- colours_sample[one_cluster]
        temp$size <- 2.5
        res_myData <- rbind(res_myData,temp)
      }
      
      clicked_cell <- event_data("plotly_click")
      
      unique(res_myData[which( res_myData$key == clicked_cell$key ),])
    }else{
      
    }
    
  })
  
  output$see_click_cluster <- renderPrint({
    print(values$def_cate_clu)
    print(event_data("plotly_selected",priority = "input"))
  })
  
  ifsplit_obj1 <- reactive({
    fluidRow(
      awesomeCheckbox(
        inputId = "ifsplit",
        label = "Split", 
        value = FALSE,
        width = '100%'
      ),
      style ='background-color:#E8F5FB'
    )
  })
  ifsplit_obj2 <- reactive({
    fluidRow()
  })
  observe({
    if(!is.null(input$selected_features_0_rows_selected) && length(input$selected_features_0_rows_selected) == 1){
      output$ifsplitui <- renderUI({ifsplit_obj1()})
    }else{
      output$ifsplitui <- renderUI({ifsplit_obj2()})
    }
  })
  
  editnamemodel_obj1 <- reactive({
    res <- click_msg_res()
    msg <- sprintf('<h3>You pick cluster <font color="#1558A0">%s</font> in <font color="#1558A0">%s</font> with <font color="#1558A0">%s</font> cells, you can:</h3>', 
                   res$click_cluster_name, res$click_category_name, res$click_cell_num)
    fluidRow(
      column(width = 6,
             h3(HTML(msg)),
             style ='padding-left:1px; padding-top:1px'
      ),
      column(width = 4,
             textInput('editclustername','Input a new cluster name as a modification:',width = '100%'),
             style ='padding-left:1px; padding-top:1px'
             ),
      column(width = 1,
             actionButton('comfirmeditclustername','Confirm!',width = '100%',
                          style = "color: white; background-color: #1f78b4;"),
             style ='padding-left:1px; padding-top:25px'
             ),
      column(width = 1,
             actionButton('closeeditclustername','Close',width = '100%',
                          style = "color: white; background-color: #e31a1c;"),
             style ='padding-left:1px; padding-top:25px'
             ),
    )
  })
  editnamemodel_obj2 <- reactive({fluidRow()})
  editnamemodel_obj3 <- reactive({
    res <- click_msg_res()
    msg <- sprintf('<h3>You pick cluster <font color="#1558A0">%s</font> in <font color="#1558A0">%s</font> with <font color="#1558A0">%s</font> cells, you can:</h3>', 
                   res$click_cluster_name, res$click_category_name, res$click_cell_num)
    fluidRow(
      column(width = 6,
             h3(HTML(msg)),
             style ='padding-left:1px; padding-top:1px'
      ),
      column(width = 6,
             fluidRow(
               column(width = 8,
                      textInput('editclustername','Input a new cluster name as a modification:',width = '100%'),
                      style ='padding-left:1px; padding-top:1px'
                      ),
               column(width = 2,
                      actionButton('comfirmeditclustername','Confirm!',width = '100%',
                                   style = "color: white; background-color: #1f78b4;"),
                      style ='padding-left:1px; padding-top:25px'
                      ),
               column(width = 2,
                      actionButton('closeeditclustername','Close',width = '100%',
                                   style = "color: white; background-color: #e31a1c;"),
                      style ='padding-left:1px; padding-top:25px'
                      )
             ),
             fluidRow(
               actionButton('delpickcluster','Delete this cluster!',width = '100%'),
             )
      )
    )
  })

  observeEvent(event_data("plotly_click",priority = "event"),{
    if(length(input$selected_features_0_rows_selected) == 1){
      output$editnamemodel <- renderUI({editnamemodel_obj2()})
    }else{
      if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
        output$editnamemodel <- renderUI({editnamemodel_obj3()})
      }else{
        output$editnamemodel <- renderUI({editnamemodel_obj1()})
      }
    }
    output$intermodel <- renderUI({intermodel_ob2()})
  })
  
  observeEvent(input$closeeditclustername,{
    output$editnamemodel <- renderUI({editnamemodel_obj2()})
  })
  observeEvent(input$comfirmeditclustername,{
    if( ! is.null(input$editclustername) && input$editclustername != ''){
      # editname时 当选择的是自定义的category
      if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
        ori_name <- names( values$def_cate_clu[[values$click_res$click_category_name]] )
        ori_name[which(ori_name == values$click_res$click_cluster_name)] <- input$editclustername
        names( values$def_cate_clu[[values$click_res$click_category_name]] ) <- ori_name
      }
      
      # 当选择的是 seurat category
      if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)){
        temp <- as.character(Idents(values$seurat_obj))
        temp_levels <- levels(Idents(values$seurat_obj))
        
        temp[which(temp == values$click_res$click_cluster_name)] <- input$editclustername
        temp_levels[which(temp_levels == values$click_res$click_cluster_name)] <- input$editclustername

        new_ident <- factor(temp, levels = temp_levels, order = F)

        names(new_ident) <- names(Idents(values$seurat_obj))
        Idents(values$seurat_obj) <- new_ident
      }
      
      # 当选择的是sample
      if(!is.null(input$input_pick_category) && input$input_pick_category == 'Sample'){
        temp <- as.character(values$seurat_obj@meta.data$orig.ident)
        temp[which(temp == values$click_res$click_cluster_name)] <- input$editclustername
        values$seurat_obj@meta.data$orig.ident <- temp
      }
    }
    output$editnamemodel <- renderUI({editnamemodel_obj2()})
  })
  observeEvent(input$delpickcluster,{
    values$def_cate_clu[[values$click_res$click_category_name]][[values$click_res$click_cluster_name]] <- NULL
    output$editnamemodel <- renderUI({editnamemodel_obj2()})
  })
  
  # 点一下，获得细胞数，cluster名称，category名称
  click_msg_res <- eventReactive(event_data("plotly_click",priority = "event"),{
    res <- list(
      click_category_name = 'str',
      click_cluster_name = 'str',
      click_cell_num = 0
    )
    
    # 当选择的是自定义的category
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
      res$click_category_name <- input$input_pick_category
      values$click_res$click_category_name <- input$input_pick_category
      
      num <- length(values$def_cate_clu[[input$input_pick_category]]) # cluster数量
      colours_sample <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- names(values$def_cate_clu[[input$input_pick_category]])
      myData <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings[,1:2])
      pick <- input$input_pick_cluster_for_def

      res_myData <- c()

      for (i in 1:length(pick) ) {
        one_cluster <- pick[i]
        cell_in_cluster <- as.character(values$def_cate_clu[[input$input_pick_category]][[one_cluster]])
        temp <- myData[cell_in_cluster,]
        temp$text <- rownames(temp)
        temp$key <- one_cluster#paste0(one_cluster,':',rownames(temp))
        temp$color <- colours_sample[one_cluster]
        temp$size <- 2.5
        res_myData <- rbind(res_myData,temp)
      }

      clicked_cell <- event_data("plotly_click")
      clicked_cluster <- strsplit(clicked_cell$key,':')[[1]][1]
      
      res$click_cluster_name <- unique(res_myData[which( res_myData$key == clicked_cluster ),'key'])
      values$click_res$click_cluster_name <- unique(res_myData[which( res_myData$key == clicked_cluster ),'key'])
      
      res$click_cell_num <- length(res_myData[which( res_myData$key == clicked_cluster ),'key'])
      values$click_res$click_cell_num <- length(res_myData[which( res_myData$key == clicked_cluster ),'key'])
    }
    
    # 当选择的是seurat category
    if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)){
      res$click_category_name <- input$input_pick_category
      values$click_res$click_category_name <- input$input_pick_category
      
      num <- length(unique( Idents(values$seurat_obj) ))
      colours_sample<-colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- as.character(levels(Idents(values$seurat_obj)))
      myData <- as.data.frame(cbind(values$seurat_obj@reductions[[input$input_pick_category]]@cell.embeddings[,1:2]))
      colnames(myData) <- c('D1','D2')
      pick <- input$input_pick_cluster
      
      res_myData <- c()
      for (i in 1:length(pick)) {
        one_cluster <- pick[i]
        temp <- myData[which(Idents(values$seurat_obj) == one_cluster),]
        temp$text <- rownames(temp)
        temp$key <- one_cluster #paste0(one_cluster,':',rownames(temp))
        temp$color <- colours_sample[one_cluster]
        temp$size <- 2.5
        res_myData <- rbind(res_myData,temp)
      }
      
      clicked_cell <- event_data("plotly_click")
      clicked_cluster <- strsplit(clicked_cell$key,':')[[1]][1]
      
      res$click_cluster_name <- (unique(res_myData[which( res_myData$key == clicked_cluster ),'key']))
      values$click_res$click_cluster_name <- (unique(res_myData[which( res_myData$key == clicked_cluster ),'key']))
      
      res$click_cell_num <- length(res_myData[which( res_myData$key == clicked_cluster ),'key'])
      values$click_res$click_cell_num <- length(res_myData[which( res_myData$key == clicked_cluster ),'key'])
    }
    
    # 当选择的是sample
    if(!is.null(input$input_pick_category) && input$input_pick_category == 'Sample'){
      res$click_category_name <- input$input_pick_category
      values$click_res$click_category_name <- input$input_pick_category
      
      num <- length( unique(values$seurat_obj@meta.data$orig.ident) )
      colours_sample <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- as.character(unique(values$seurat_obj@meta.data$orig.ident))
      myData <- as.data.frame(values$seurat_obj@reductions$umap@cell.embeddings[,1:2])
      colnames(myData) <- c('D1','D2')
      
      myData$Sample <- values$seurat_obj@meta.data$orig.ident
      pick <- input$input_pick_cluster_for_sample
      
      res_myData <- c()
      for (i in 1:length(pick)) {
        one_cluster <- pick[i]
        temp <- myData[which(values$seurat_obj@meta.data$orig.ident == one_cluster),]
        temp$text <- rownames(temp)
        temp$key <- one_cluster # paste0(one_cluster,':',rownames(temp))
        temp$color <- colours_sample[one_cluster]
        temp$size <- 2.5
        res_myData <- rbind(res_myData,temp)
      }
      clicked_cell <- event_data("plotly_click")
      clicked_cluster <- strsplit(clicked_cell$key,':')[[1]][1]
      
      res$click_cluster_name <- (unique(res_myData[which( res_myData$key == clicked_cluster ),'key']))
      values$click_res$click_cluster_name <- (unique(res_myData[which( res_myData$key == clicked_cluster ),'key']))
      res$click_cell_num <- length(res_myData[which( res_myData$key == clicked_cluster ),'key'])
      values$click_res$click_cell_num <- length(res_myData[which( res_myData$key == clicked_cluster ),'key'])
    }
    
    return(res)
  })
  
  output$cellmap <- renderPlotly({
    r1 <- cellmap_obj1()
    r2 <- cellmap_obj2()
    r3 <- cellmap_obj3()
    r4 <- cellmap_obj4()
    r5 <- cellmap_obj5()
    
    if( length(input$selected_features_0_rows_selected) == 1 ){
      if( !is.null(input$ifsplit) && input$ifsplit == T ){
        return(r5)
      }else{
        return(r3)
      }
    }else{
      if(!is.null(input$input_pick_category) && input$input_pick_category == 'Sample'){
        return(r1)
      }
      if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && input$input_pick_category %in% names(values$seurat_obj@reductions)){
        return(r2)
      }
      if(!is.null(input$input_pick_category) && input$input_pick_category != "Sample" && ! input$input_pick_category %in% names(values$seurat_obj@reductions)){
        return(r4)
      }
    }
  })
  
  output$selectsplitcategory <- renderUI({
    if( !is.null(input$ifsplit) && input$ifsplit == T){
      selectInput(
        inputId = 'selectsplitcategory',
        label = 'Split by:',
        choices = values$category_name,
        selected = NULL,
        selectize = TRUE,
        width = '100%'
      )
    }else{
      return(NULL)
    }
  })
  
  output$inSelect <- renderUI({
    x <- rownames(values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]]@data)
    if(is.null(x)){
      x <- character(0)
    }
    selectInput(inputId = "inSelect",
                label = "Search and select genes",
                choices = x,
                selected = NULL,
                selectize = T
    )
  })

  output$monitor1 <- renderText({
    select_points <- event_data("plotly_selected",priority = "input")
    msg <- NULL
    if( length(select_points) == 0 ){
      msg <- "no cell selected.\n"
    }else{
      cell_number <- length(select_points[,1])
      if(cell_number > 1){
        msg <- paste0(cell_number," cells were selected.\n",collapse = '')
      }else if(cell_number == 1){
        msg <- paste0(cell_number," cell was selected.\n",collapse = '')
      }else{
        msg <- "The number of cells in this cluster is 1.\n"
      }
    }
    return(msg)
  })

  intermodel_ob1 <- reactive({
    fluidRow(
      h2("Defind your category and cluster"),
      fluidRow(
        h3(textOutput("monitor1")),
        h3("You can:"),
        style ='padding-left:20px'
      ),
      fluidRow(
        column(width = 5,
               textInput("categorynametxt", "Enter a new category name:",width = '100%'),
               style ='padding-right:1px'),
        column(width = 2,
               actionButton(inputId = 'addcategory',
                            label = 'Add category',
                            width = '100%',
                            style = "color: white; background-color: #a6cee3;"
               ),
               style ='padding-left:1px; padding-top:25px'),
        column(width = 5,
               uiOutput('selectcategory'))
      ),
      fluidRow(
        column(width = 5,
               textInput("clusternametxt", "Enter a new cluster name:",width = '100%'),
               style ='padding-right:1px'),
        column(width = 2,
               actionButton(inputId = 'addcluster',
                            label = 'Add cluster',
                            width = '100%',
                            style = "color: white; background-color: #a6cee3;"
               ),
               style ='padding-left:1px; padding-top:25px'),
        column(width = 5,
               uiOutput('selectcluster'))
      ),
      fluidRow(
        column(
          width = 6,
          actionButton(
            inputId = "confirmdef",
            label = "Confirm your definition of selected cells",
            width = '100%',
            style = "color: white; background-color: #1f78b4;"
          )
        ),
        column(
          width = 6,
          actionButton(
            inputId = "closedef",
            label = "Close",
            width = '100%',
            style = "color: white; background-color: #e31a1c;"
          )
        )
      )
    )
  })
  intermodel_ob2 <- reactive({
    fluidRow()
  })
  
  observeEvent(event_data("plotly_selected",priority = "event"),{
    output$intermodel <- renderUI({intermodel_ob1()})
    output$editnamemodel <- renderUI({editnamemodel_obj2()})
  })
  
  observeEvent(input$confirmdef,{
    
    if(!is.null(input$selectcategory) && !is.null(input$selectcluster)){
      
      split_res <- unlist( strsplit(as.character(event_data("plotly_selected")$key),':'))
      pick_cell <- split_res[seq(2,length(split_res),2)]
      
      cate_name <- names(values$def_cate_clu[[input$selectcategory]])
      for (one in cate_name) {
        values$def_cate_clu[[input$selectcategory]][[one]] <-  setdiff(values$def_cate_clu[[input$selectcategory]][[one]],pick_cell)
      }
      
      values$def_cate_clu[[input$selectcategory]][[input$selectcluster]] <- pick_cell
      output$intermodel <- renderUI({intermodel_ob2()})
      need_to_add <- names(values$def_cate_clu)[which(!names(values$def_cate_clu) %in% values$category_name)]
      if(length(which(!names(values$def_cate_clu) %in% values$category_name)) > 0){
        for (add_one in need_to_add) {
          if( length(values$def_cate_clu[[add_one]]) > 0 ){
            values$category_name <- c(values$category_name,add_one)
          }
        }
      }
      
    }
    
  })
  
  observeEvent(input$closedef,{
    output$intermodel <- renderUI({intermodel_ob2()})
  })
  
  output$genelistpickmodel <- renderUI({
    fluidRow(
      column(width = 7,
             pickerInput(
               inputId = "pick_active_feature_list",
               label = "Default Feature List",
               choices = names(values$def_gene_list),
               selected = NULL,#names(values$def_gene_list)[length(names(values$def_gene_list))],
               options = list(style = "btn-primary"),
               width = "100%"
             ),
             style ='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px'
      ),
      column(width = 2,
             actionButton(
               inputId = "addgenelistname",
               label = NULL,
               status = "primary",
               style = "unite",
               icon = icon("plus-sign",lib="glyphicon"),
               width = "100%"
             ),
             style ='padding-left:0px; padding-right:0px; padding-top:25px; padding-bottom:0px'
      ),
      column(width = 2,
             actionButton(
               inputId = "delgenelistname",
               label = NULL,
               status = "primary",
               style = "unite",
               icon = icon("minus-sign",lib="glyphicon"),
               width = "100%"
             ),
             style ='padding-left:0px; padding-right:0px; padding-top:25px; padding-bottom:0px'
      ),
      style ='padding-left:1px; padding-right:1px; padding-top:1px; padding-bottom:1px'
    )
  })
  
  addgenelisenamemodel_obj1 <- reactive({
    fluidRow(
      column(width = 7,
             textInput(inputId = "inputnewgenelistname",label = NULL,value = "New list name", width = "100%"),
             style ='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px'
      ),
      column(width = 4,
             actionButton(inputId = "confirmaddgenelist",label = "Add!",width = "100%"),
             style ='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px'
      )
    )
  })
  addgenelisenamemodel_obj2 <- reactive({
    fluidRow()
  })
  # values$def_gene_list[[input$inputnewgenelistname]]
  # values$def_gene_list[[input$pick_active_feature_list]]
  
  pick_gene_list <-  reactive({
    input$pick_active_feature_list
  })
  
  observeEvent(input$addgenelistname,{
    if(!values$if_show_add_gene_list){
      output$addgenelisenamemodel <- renderUI({addgenelisenamemodel_obj1()})
      values$if_show_add_gene_list <- T
    }else{
      output$addgenelisenamemodel <- renderUI({addgenelisenamemodel_obj2()})
      values$if_show_add_gene_list <- F
    }
  })
  observeEvent(input$confirmaddgenelist,{
    output$addgenelisenamemodel <- renderUI({addgenelisenamemodel_obj2()})
    values$if_show_add_gene_list <- F
  })
  observeEvent(input$confirmaddgenelist,{
    if(!input$inputnewgenelistname %in% names(values$def_gene_list)){
      values$def_gene_list[[input$inputnewgenelistname]] <- ''
      updatePickerInput(session,
                        inputId = "pick_active_feature_list",
                        selected = input$inputnewgenelistname
                        )
    }
  })
  observeEvent(input$delgenelistname,{
    if( pick_gene_list() != 'Default Feature List' ){
      values$def_gene_list[[pick_gene_list()]] <- NULL
    }
  })
  
  output$selectcategory <- renderUI({
    selectInput(inputId = "selectcategory",
                label = "Select a category:",
                width = '100%',
                choices = names(values$def_cate_clu),
                selected = names(values$def_cate_clu)[length(names(values$def_cate_clu))],
                selectize = T)
  })

  output$selectcluster <- renderUI({
    selectInput(inputId = "selectcluster",
                label = "Select a cluster:",
                width = '100%',
                choices = names(values$def_cate_clu[[input$selectcategory]]),
                selected = names(values$def_cate_clu[[input$selectcategory]])[length(names(values$def_cate_clu[[input$selectcategory]]))],
                selectize = T)
  })

  observeEvent(input$addcategory,{
    if(!is.null(input$categorynametxt)){
      values$def_cate_clu[[input$categorynametxt]] <- list()
    }
  })

  observeEvent(input$addcluster,{
    if(!is.null(input$clusternametxt)){
      values$def_cate_clu[[input$selectcategory]][[input$clusternametxt]] <- 'N'
    }
  })

  # output$monitor2 <- renderPrint({
  #   print( event_data("plotly_selected",priority = "input"))
  # })

  observe({
    if(is.null(savevalues())){
      if("integrated" %in% names(values$seurat_obj@assays)){
        values$category_name <- c(names(values$seurat_obj@reductions),"Sample")
      }else{
        values$category_name <- names(values$seurat_obj@reductions)
      }   
    }
  })

  output$pick_category <- renderUI({
    radioGroupButtons(
      inputId = "input_pick_category",
      label = "Choice Category:",
      choices = values$category_name,
      selected =  'umap',
      justified = TRUE,
      checkIcon = list( yes = icon("ok",lib = "glyphicon")) )
  })
  
  pick_category_all_obj1 <- reactive({actionButton("delpickcategory","Delete this category!",width = '100%')})
  pick_category_all_obj2 <- reactive({fluidRow()})
  
  observeEvent(input$input_pick_category,{
    if( ! input$input_pick_category %in% c(names(values$seurat_obj@reductions),"Sample")){
      output$pick_category_all <- renderUI(pick_category_all_obj1())
    }else{
      output$pick_category_all <- renderUI(pick_category_all_obj2())
    }
  })
  
  observeEvent(input$delpickcategory,{
    temp <- input$input_pick_category
    values$def_cate_clu[[input$input_pick_category]] <- NULL
    values$category_name <- values$category_name[-which(values$category_name == temp)]
    output$pick_category_all <- renderUI({fluidRow()})
  })

  output$pick_cluster <- renderUI({
    if(input$input_pick_category  == "Sample"){
      cluster_name <- unique(values$seurat_obj@meta.data$orig.ident)
      selected <- unique(values$seurat_obj@meta.data$orig.ident)
      checkboxGroupButtons(
        inputId = "input_pick_cluster_for_sample",
        label = "Pick Cluster:",
        choices = cluster_name,
        selected = selected,
        checkIcon = list(
          yes = tags$i(class = "fa fa-check-square",
                       style = "color: steelblue"),
          no = tags$i(class = "fa fa-square-o",
                      style = "color: steelblue"))
      )
    }else if(input$input_pick_category %in% names(values$seurat_obj@reductions)){
      cluster_name <- levels(Idents(values$seurat_obj))
      selected <- levels(Idents(values$seurat_obj))
      checkboxGroupButtons(
        inputId = "input_pick_cluster",
        label = "Pick Cluster:",
        choices = cluster_name,
        selected = selected,
        checkIcon = list(
          yes = tags$i(class = "fa fa-check-square",
                       style = "color: steelblue"),
          no = tags$i(class = "fa fa-square-o",
                      style = "color: steelblue"))
      )
    }else{
      cluster_name <- names(values$def_cate_clu[[input$input_pick_category]])
      selected <- names(values$def_cate_clu[[input$input_pick_category]])
      checkboxGroupButtons(
        inputId = "input_pick_cluster_for_def",
        label = "Pick Cluster:",
        choices = cluster_name,
        selected = selected,
        checkIcon = list(
          yes = tags$i(class = "fa fa-check-square",
                       style = "color: steelblue"),
          no = tags$i(class = "fa fa-square-o",
                      style = "color: steelblue"))
      )
    }
  })
  
  output$selectcolorforpoints <- renderUI({
    selectInput(
      inputId = 'selectcolorforpoints',
      label = 'Select color for points:',
      choices =  rownames(brewer.pal.info[brewer.pal.info$category == 'qual',]),
      selected = NULL,
      selectize = TRUE,
      width = '100%'
    )
  })
  
  observeEvent(input$inSelect,{
    mygenelist <- isolate({values$def_gene_list[[pick_gene_list()]]})
    mygenelist <- mygenelist[which(mygenelist != '')]
    if(input$inSelect != "" && !input$inSelect %in% mygenelist){
      values$def_gene_list[[pick_gene_list()]] <- c(mygenelist,input$inSelect)
      updatePickerInput(session,
                        inputId = "pick_active_feature_list",
                        selected = pick_gene_list()
      )
    }
  })

  observeEvent(input$clean_my_gene_list,{
    if (!is.null(input$selected_features_0_rows_selected)) {
      values$def_gene_list[[pick_gene_list()]] <- values$def_gene_list[[pick_gene_list()]][-as.numeric(input$selected_features_0_rows_selected)]
    }
  })

  output$selected_features_0 <- renderDT({
    res <- values$def_gene_list[[pick_gene_list()]][which(values$def_gene_list[[pick_gene_list()]] != '')]
    DT::datatable(as.data.frame(res),
                  rownames = F,
                  selection = 'single',
                  escape = FALSE,
                  colnames = 'My gene list',
                  options = list(sDom  = 't<"bottom"><"clear">', #<"top">
                                 pageLength = 20))
  })

  output$output_my_gene_list <- downloadHandler(
    filename = function(){
      'scui_gene_list.txt'},
    content = function(con){
      write.table(as.matrix(values$def_gene_list[[pick_gene_list()]]),
                 con,quote = FALSE,sep = '\t',
                 row.names = FALSE,col.names = FALSE)
    }
  )
  
  output$output_feature_table <- downloadHandler(
    filename = function(){
      'scui_feature_table.csv'},
    content = function(con){
      res <- NULL
      if(!is.null(input$input_ifusemycategory) && input$input_ifusemycategory){
        res <- findmarkers_res()
      }else{
        res <- ori_findmarkers_res()
      }
      write.csv(res,con,quote = FALSE,row.names = FALSE,col.names = FALSE)
    }
  )
  
  output$ifusemycategory <- renderUI({
    if( is.list(values$def_cate_clu) && length(values$def_cate_clu) == 0 ){
      NULL
    }else if(is.list(values$def_cate_clu) && length(values$def_cate_clu) > 0){
      fluidRow(
        box(width = 12,
          h3("Check if find marker genes according to 'My Category':"),
          awesomeCheckbox(
            inputId = "input_ifusemycategory",
            label = "Use my category", 
            value = FALSE,
            width = '100%'
          ))
        )
    }else{
      NULL
    }
  })
  
  output$categoryforgetmarker <- renderUI({
    if(! is.null(input$input_ifusemycategory) && input$input_ifusemycategory){
      fluidRow(
        uiOutput('monitor_for_findmarkers'),
        column(width = 6,
               selectInput(
                 inputId = "input_selectcategoryforgetmarker",
                 label = "Select Category",
                 choices = names(values$def_cate_clu)
               ))
      )
    }else{
      fluidRow()
    }
  })
  
  findmarkers_res <- reactive({
    req(input$input_ifusemycategory)
    cluster_list <- values$def_cate_clu[[input$input_selectcategoryforgetmarker]]
    # 拿出细胞
    cells <- as.character(unlist(cluster_list))
    temp_obj <- values$seurat_obj[VariableFeatures(object = values$seurat_obj),cells]
    
    # 配置group
    res <- NULL
    for (i in 1:length(cluster_list)) {
      res <- c(res,rep(names(cluster_list[i]),length(cluster_list[[i]])))
    }
    temp_obj@meta.data[['marker_group']] <- res
    
    # 给每个
    if(length(cluster_list) >= 2){
      myfindmarker_res <- c()
      for (i in 1:length(cluster_list)) {
        one <- FindMarkers(temp_obj,ident.1 = names(cluster_list[i]),
                           group.by = 'marker_group', only.pos = T)
        one$marker_group <- names(cluster_list[i])
        myfindmarker_res <- rbind(myfindmarker_res,one)
      }
      
      myfindmarker_res <- as.data.frame(myfindmarker_res)
      is.num <- sapply(myfindmarker_res, is.numeric)
      myfindmarker_res[is.num] <- lapply(myfindmarker_res[is.num], signif, 3)
      myfindmarker_res$gene <- rownames(myfindmarker_res)
      
      myfindmarker_res$marker_group <- as.factor(myfindmarker_res$marker_group)
      colnames(myfindmarker_res)[2] <- 'avg_logFC'
      values$marker_res_save[['findmarkers_res']] <- myfindmarker_res
      return(myfindmarker_res)
    }else{
      NULL
    }
  })
  ori_findmarkers_res <- reactive({
    pick <- obj_list()$seurat_obj_markers %>% group_by(cluster) %>% arrange(desc(avg_logFC))
    pick <- as.data.frame(pick)
    is.num <- sapply(pick, is.numeric)
    pick[is.num] <- lapply(pick[is.num], signif, 3)
    values$marker_res_save[['ori_findmarkers_res']] <- pick
    return(pick)
  })
  
  output$monitor_for_findmarkers <- renderUI({
    if(is.null(values$def_cate_clu)){
      return(box(width = 12,h3('')))
    }else{
      cluster_list <- values$def_cate_clu[[input$input_selectcategoryforgetmarker]]
      if(input$input_ifusemycategory){
        if(length(cluster_list) < 2){
          return(box(width = 12,h3('Need more then 2 clusters for my category!')))
        }else{
          return(box(width = 12,h3('Run FindMarkers Complete!')))
        }
      }else{
        return(box(width = 12,h3('')))
      }    
    }
  })
  
  output$features_datatable <- renderDT({
    if(!is.null(input$input_ifusemycategory) && input$input_ifusemycategory){
      DT::datatable(findmarkers_res(),rownames= F,filter = 'top',
                    options = list(pageLength = 15,lengthMenu = c(5,10,15,25,50)))
    }else{
      DT::datatable(ori_findmarkers_res(),rownames= F,filter = 'top',options = list(
        pageLength = 15,
        lengthMenu = c(5,10,15,25,50)
      ))
    }
  })

  add_one_line <- reactive({
    pick <- obj_list()$seurat_obj_markers %>% group_by(cluster) %>% arrange(desc(avg_logFC))
    return(pick$gene[input$features_datatable_rows_selected])
  })

  output$selected_features <- renderTable({
    add_one_line()
  })

  output$monitor3 <- renderPrint({
    cat(add_one_line())
  })

  DTproxy <- dataTableProxy("features_datatable")
  observeEvent(input$clean_list_switch,{
    DTproxy %>% selectRows(NULL)
  })

  observeEvent(input$add_list_switch,{
    now_list <- add_one_line()
    if( length(which(!now_list %in% values$def_gene_list[[pick_gene_list()]])) > 0 ){
      newtab <- switch(input$tabs,"CellBrowser" = "Features","Features"="CellBrowser")
      updateTabItems(session, "tabs", newtab)
      temp <- pick_gene_list()
      updatePickerInput(session,
                        inputId = "pick_active_feature_list",
                        selected = temp
      )
      values$def_gene_list[[temp]] <- c(values$def_gene_list[[temp]],now_list[which(!now_list %in% values$def_gene_list[[temp]])])
      if(length(which(values$def_gene_list[[temp]] == ''))>0){
        values$def_gene_list[[temp]] <- values$def_gene_list[[temp]][-which(values$def_gene_list[[temp]] == '')]
      }
    }
  })

  output$selecttopormygenelist <- renderUI({
    options <- NULL
    if(length(values$def_gene_list[[pick_gene_list()]]) == 0 || values$def_gene_list[[pick_gene_list()]] == ''){
      options <- 'Top gene'
    }else{
      options <- c( 'Top gene','My gene list' )
    }
    radioGroupButtons(
      inputId = "selecttopormygenelist",
      label = "Features come from:",
      choices = options,
      justified = TRUE,
      selected = 'Top gene',
      checkIcon = list(
        yes = icon("ok", lib = "glyphicon"))
    )
  })
  
  output$topnforheatmap <- renderUI({
    if(input$selecttopormygenelist == 'Top gene'){
      selectInput(inputId = "topnforheatmap",
                  label = "Top n for heatmap:",
                  width = '100%',
                  choices = 1:15,
                  selectize = T)
    }else{
      return(NULL)
    }
  })
  
  
  
  output$categoryforheatmap <- renderUI({
    options <- NULL
    if( length(names(values$def_cate_clu)) == 0 ){
      options <- 'Seurat cluster'
    }else{
      if(input$input_ifusemycategory){
        options <- c('Seurat cluster',names(values$def_cate_clu))
      }else{
        options <- 'Seurat cluster'
      }
    }
    radioGroupButtons(
      inputId = "categoryforheatmap",
      label = "Cells come from:",
      choices = options,
      justified = TRUE,
      selected = 'Seurat cluster',
      checkIcon = list(
        yes = icon("ok", lib = "glyphicon"))
    )
  })
  
  output$clusterforheatmap <- renderUI({
    options <- NA
    selected <- NA
    if(input$categoryforheatmap == 'Seurat cluster'){
      options <- as.character(levels(Idents(values$seurat_obj)))
      selected <- as.character(levels(Idents(values$seurat_obj)))
    }else{
      options <- names(values$def_cate_clu[[input$categoryforheatmap]])
      selected <- names(values$def_cate_clu[[input$categoryforheatmap]])
    }
    selectInput(inputId = "clusterforheatmap",
                label = "Cluster for heatmap:",
                width = '100%',
                choices = options,
                selected = selected,
                multiple = T,
                selectize = T)
  })
  
  output$checkboxforheatmapcluster <- renderUI({
    awesomeCheckbox(
      inputId = "if_cluster_for_heatmap",
      label = "If cluster for heatmap", 
      value = T
    )
  })
  
  output$pickcolorforheatmap <- renderUI({
    selectInput(
      inputId = 'selectcolorforheatmap',
      label = 'Select color for heatmap:',
      choices =  c('viridis',rownames(brewer.pal.info[brewer.pal.info$category != 'qual',])),
      selected = NULL,
      selectize = TRUE,
      width = '100%'
    )
  })
  
  plotmodel_stat_obj1 <- reactive({
    fluidRow(
      column(width = 4, sliderInput("plotheight_stat", "Plot output height (px)",min = 0, max = 2000, value = 500, step = 100)),
      column(width = 4, sliderInput("plotwidth_stat", "Plot output width (px)",min = 0, max = 2000, value = 500, step = 100)),
      column(width = 4, radioGroupButtons(inputId = "plotformat_stat",label = "Plot format",choices = c("svg", "png"),
                               checkIcon = list(yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))))
    )
  })
  plotmodel_stat_obj2 <- reactive({fluidRow()})
  observeEvent(input$tabs,{
    if(input$tabs == 'Heatmap'){
      output$plotmodel_stat <- renderUI({plotmodel_stat_obj1()})
    }else{
      output$plotmodel_stat <- renderUI({plotmodel_stat_obj2()})
    }
  })
  
  plotmodel_stat_v_obj1 <- reactive({
    fluidRow(
      column(width = 4, sliderInput("plotheight_stat_v", "Plot output height (px)",min = 0, max = 2000, value = 500, step = 100)),
      column(width = 4, sliderInput("plotwidth_stat_v", "Plot output width (px)",min = 0, max = 2000, value = 500, step = 100)),
      column(width = 4, radioGroupButtons(inputId = "plotformat_stat_v",label = "Plot format",choices = c("svg", "png"),
                                          checkIcon = list(yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                           no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))))
    )
  })
  plotmodel_stat_v_obj2 <- reactive({fluidRow()})
  observeEvent(input$tabs,{
    if(input$tabs == 'Violin'){
      output$plotmodel_stat_v <- renderUI({plotmodel_stat_v_obj1()})
    }else{
      output$plotmodel_stat_v <- renderUI({plotmodel_stat_v_obj2()})
    }
  })
  
  observeEvent(c(input$plotheight_stat,input$plotwidth_stat,input$plotformat_stat),{
    click("heat_map_run")
  })
  observeEvent(c(input$plotheight_stat_v,input$plotwidth_stat_v,input$plotformat_stat_v),{
    click("violin_map_run")
  })

  heatmapdata_obj1 <- eventReactive(input$heat_map_run,{
    # req(input$plotformat, input$plotheight, input$plotwidth)
    data <- values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]]@data
    features <- NULL # 'Top gene','My gene list'
    if( input$selecttopormygenelist == 'Top gene' ){ # 如果选择的是top gene
      if( input$categoryforheatmap == 'Seurat cluster' ){ # 自己的top
        features <- values$marker_res_save[['ori_findmarkers_res']] %>% group_by(cluster) %>% top_n(n = as.numeric(input$topnforheatmap), wt = avg_logFC)
        features <- features$gene
      }else{ # 人为的top
        features <- values$marker_res_save[['findmarkers_res']] %>% group_by(marker_group) %>% top_n(n = as.numeric(input$topnforheatmap), wt = avg_logFC)
        features <- features$gene
      }
    }else{ # 如果选择的不是top gene
      features <- values$def_gene_list[[pick_gene_list()]]
    }
    values$heatmap_height <- length(features)
    
    data <- data[features,]
    
    mat <- matrix(0)
    if(input$categoryforheatmap == 'Seurat cluster'){ # 按照原版的分类
      mat <- matrix(0,
                    nrow = length(features),
                    ncol = length(input$clusterforheatmap))
      rownames(mat) <- features
      colnames(mat) <- input$clusterforheatmap
      
      for (i in 1:length(mat[1,])) {
        one_cluster <- input$clusterforheatmap[i]
        one_cluser_data <- data[,Idents(values$seurat_obj) == one_cluster]
        mat[,i] <- rowMeans(as.matrix(one_cluser_data))
      }
    }else{ # 按照新的分类
      cluster_list <- values$def_cate_clu[[input$categoryforheatmap]]
      mat <- matrix(0,
                    nrow = length(features),
                    ncol = length(cluster_list))
      rownames(mat) <- features
      colnames(mat) <- names(cluster_list)
      
      for (i in 1:length(mat[1,])) {
        one_cluster <- names(cluster_list)[i]
        one_cluser_data <- data[,values$def_cate_clu[[input$categoryforheatmap]][[one_cluster]] ]
        mat[,i] <- rowMeans(as.matrix(one_cluser_data))
      }
    }
    
    pick_color <- NULL
    if(input$selectcolorforheatmap == 'viridis'){
      pick_color <- NULL
      if( length(mat[1,]) == 1 ){
        heatmaply(mat,
                  Rowv = NULL,
                  Colv = NULL,
                  fontsize_row = 6,
                  show_dendrogram = c(F,T)
        ) %>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat_stat, height = input$plotheight_stat, width = input$plotwidth_stat,filename=paste0("heatmap")))
      }else{
        heatmaply(mat,
                  Rowv=NULL,
                  Colv = input$if_cluster_for_heatmap,
                  fontsize_row = 6,
                  show_dendrogram = c(F,T)
        )%>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat_stat, height = input$plotheight_stat, width = input$plotwidth_stat,filename=paste0("heatmap")))
      }
    }else{
      pick_color <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforheatmap,'maxcolors'], input$selectcolorforheatmap))(25)
      if( length(mat[1,]) == 1 ){
        heatmaply(mat,
                  colors = pick_color,
                  Rowv = NULL,
                  Colv = NULL,
                  fontsize_row = 6,
                  show_dendrogram = c(F,T)
        ) %>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat_stat, height = input$plotheight_stat, width = input$plotwidth_stat,filename=paste0("heatmap")))
      }else{
        heatmaply(mat,
                  colors = pick_color,
                  Rowv=NULL,
                  Colv = input$if_cluster_for_heatmap,
                  fontsize_row = 6,
                  show_dendrogram = c(F,T)
        )%>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
          config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat_stat, height = input$plotheight_stat, width = input$plotwidth_stat,filename=paste0("heatmap")))
      }
    }
    

  })
  heatmapdata_obj2 <- reactive({fluidRow()})
  output$heat_map <- renderPlotly({
    if(!is.null(input$heat_map_run) && input$heat_map_run != 0){
      ggplotly(heatmapdata_obj1())
    }else{
      heatmapdata_obj2()
    }
  })

  output$all_heat_map_ui <- renderUI({
    if(!is.null(input$heat_map_run) && input$heat_map_run != 0){
      fact <- 300
      if(values$heatmap_height < 30){
        fact <- 300
      }else{
        fact <- values$heatmap_height*10
      }
      fluidRow(
        shinycssloaders::withSpinner(
          plotlyOutput('heat_map',height = paste0(fact,'px'))
        )
      )
     }else{
       fluidRow()
      }

  })
  
  output$selecttopormygenelist_for_violin <- renderUI({
    options <- NULL
    if(length(values$def_gene_list[[pick_gene_list()]]) == 0 || values$def_gene_list[[pick_gene_list()]] == ''){
      options <- 'Top gene'
    }else{
      options <- c( 'Top gene','My gene list' )
    }
    radioGroupButtons(
      inputId = "selecttopormygenelist_for_violin",
      label = "Features come from:",
      choices = options,
      justified = TRUE,
      selected = 'My gene list',
      checkIcon = list(yes = icon("ok", lib = "glyphicon"))
    )
  })
  
  output$categoryforviolin <- renderUI({
    options <- NULL
    if( length(names(values$def_cate_clu)) == 0 ){
      options <- 'Seurat cluster'
    }else{
      if(input$input_ifusemycategory){
        options <- c('Seurat cluster',names(values$def_cate_clu))
      }else{
        options <- 'Seurat cluster'
      }
    }
    radioGroupButtons(
      inputId = "categoryforviolin",
      label = "Cells come from:",
      choices = options,
      justified = TRUE,
      selected = 'Seurat cluster',
      checkIcon = list(
        yes = icon("ok", lib = "glyphicon"))
    )
  })
  
  output$clusterforviolin <- renderUI({
    options <- NA
    selected <- NA
    if(input$categoryforviolin == 'Seurat cluster'){
      options <- as.character(levels(Idents(values$seurat_obj)))
      selected <- as.character(levels(Idents(values$seurat_obj)))
    }else{
      options <- names(values$def_cate_clu[[input$categoryforviolin]])
      selected <- names(values$def_cate_clu[[input$categoryforviolin]])
    }
    selectInput(inputId = "clusterforviolin",
                label = "Cluster for violin:",
                width = '100%',
                choices = options,
                selected = selected,
                multiple = T,
                selectize = T)
  })
  
  output$selectviolinfeature <- renderUI({
    if(input$selecttopormygenelist_for_violin == "My gene list"){
      selectInput(
        inputId = 'input_selectviolinfeature',
        label = "Violin Features",
        choices = values$def_gene_list[[pick_gene_list()]],
        selected = values$def_gene_list[[pick_gene_list()]][1],
        multiple = T
      )
    }else if(input$selecttopormygenelist_for_violin == 'Top gene'){
      selectInput(inputId = "topnforviolin",
                  label = "Top n for violin:",
                  width = '100%',
                  choices = 1:15,
                  selectize = T)
    }else{
      return(NULL)
    }
  })

  output$split_sample_check_box <- renderUI({
    if("integrated" %in% names(values$seurat_obj@assays) ){
      awesomeCheckbox(
        inputId = "if_split_sample",
        label = "Split sample", 
        value = F
      )
    }else{
      return(NULL)
    }
  })

  violin_output <- eventReactive(input$violin_map_run,{
    anno_gen <- function(title){
      f <- list( # family = " Courier New , monospace ",
        size = 14,color = "black")
      a <- list(text = title,font = f,xref = "paper",yref = "paper",
        yanchor = "bottom",xanchor = "center",align = "center",
        x = 0.5,y = 1,showarrow = FALSE)
      return(a)
    }
    
    pick_features <- NULL
    if(input$selecttopormygenelist_for_violin == "My gene list"){
      pick_features <- input$input_selectviolinfeature
    }else{ # top gene
      if( input$categoryforviolin == 'Seurat cluster' ){ # 原来类的top基因
        pick_features <- values$marker_res_save[['ori_findmarkers_res']] %>% group_by(cluster) %>% top_n(n = as.numeric(input$topnforviolin), wt = avg_logFC)
        pick_features <- pick_features$gene
      }else{ # 自定义类的top基因
        pick_features <- values$marker_res_save[['findmarkers_res']] %>% group_by(marker_group) %>% top_n(n = as.numeric(input$topnforviolin), wt = avg_logFC)
        pick_features <- pick_features$gene
      }
    }
    values$violin_height <- length(pick_features)
    
    colours_sample <- NULL
    if(input$categoryforviolin == 'Seurat cluster'){
      num <- length(unique(Idents(values$seurat_obj)))
      colours_sample <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- as.character(levels(Idents(values$seurat_obj)))
    }else{
      cluster_list <- values$def_cate_clu[[input$categoryforviolin]]
      num <- length(cluster_list)
      colours_sample <- colorRampPalette(brewer.pal(brewer.pal.info[input$selectcolorforpoints,'maxcolors'],input$selectcolorforpoints))(num)[num:1]
      names(colours_sample) <- names(cluster_list)
    }
    
    all_G <- list()
    if("integrated" %in% names(values$seurat_obj@assays)){
      # 有多样本
      if(input$if_split_sample == F){
        # 不分割
        for (i in 1:length(pick_features)) {
          feature <- pick_features[i]
          data <- NULL
          if(input$categoryforviolin == 'Seurat cluster'){
            data <- data.frame(
              exp = values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]]@data[feature,],
              cluster = as.character(Idents(values$seurat_obj)),
              Sample = values$seurat_obj@meta.data$orig.ident)
            
              data <- data[which(as.character(data$cluster) %in% input$clusterforviolin),]
              #data$cluster <- factor(as.character(data$cluster),levels = levels(Idents(values$seurat_obj))[which(levels(Idents(values$seurat_obj)) %in% as.character(data$cluster))])
              data$cluster <- factor(as.character(data$cluster),levels = input$clusterforviolin)
              
          }else{ # 如果使用的是自己定义的category
            #gene在细胞群内的表达量  #每个细胞的cluster情况 #每个细胞的sample
            
            # 得到cluster列表
            cluster_list <- values$def_cate_clu[[input$categoryforviolin]]
            # 得到所有选择的细胞
            cells <- as.character(unlist(cluster_list))
            # 得到这些细胞的类别
            res <- NULL
            for (i in 1:length(cluster_list)) {
              res <- c(res,rep(names(cluster_list[i]),length(cluster_list[[i]])))
            }
            temp_obj <- values$seurat_obj[,cells]
            exp <- temp_obj@assays[[DefaultAssay(values$seurat_obj)]]@data[feature,]
            # 得到这些细胞来自哪个sample
            Sample <- temp_obj@meta.data$orig.ident
            # 得到这些细胞的exp
            data <- data.frame(
              exp = exp,
              cluster = res,
              Sample = Sample
            )
          }
          p <- plot_ly() %>%
            add_trace(
              data = data,
              x = ~cluster,
              y = ~exp,
              split = ~cluster,
              type = 'violin',
              color = ~cluster,
              colors = colours_sample,
              box = list(visible = T),
              meanline = list(visible = T)
            )  %>% 
            layout(
              annotations = anno_gen(feature),
              xaxis = list(title = 'Cluster',tickangle = -45),
              yaxis = list(title = 'Expression Level',zeroline = T)
            ) %>% hide_legend()
          all_G[[(length(all_G)+1)]] <- p
        }
      }else if(input$if_split_sample == T){
        # 分割
        samples <- unique(values$seurat_obj@meta.data$orig.ident)
        all_G <- list()
        for (i in 1:length(pick_features)) {
          feature <- pick_features[i]
          showlegend <- NULL
          if(i == length(pick_features)){
            showlegend <- T
          }else{
            showlegend <- F
          }
          data <- NULL
          if(input$categoryforviolin == 'Seurat cluster'){
            data <- data.frame(
              exp = values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]]@data[feature,],
              cluster = as.character(Idents(values$seurat_obj)),
              Sample = values$seurat_obj@meta.data$orig.ident)
            data <- data[which(as.character(data$cluster) %in% input$clusterforviolin),]
            #data$cluster <- factor(as.character(data$cluster),levels = levels(Idents(values$seurat_obj))[which(levels(Idents(values$seurat_obj)) %in% as.character(data$cluster))])
            data$cluster <- factor(as.character(data$cluster),levels = input$clusterforviolin)
          }else{
            cluster_list <- values$def_cate_clu[[input$categoryforviolin]]
            cells <- as.character(unlist(cluster_list))
            res <- NULL
            for (i in 1:length(cluster_list)) {
              res <- c(res,rep(names(cluster_list[i]),length(cluster_list[[i]])))
            }
            temp_obj <- values$seurat_obj[,cells]
            exp <- temp_obj@assays[[DefaultAssay(values$seurat_obj)]]@data[feature,]
            # 得到这些细胞来自哪个sample
            Sample <- temp_obj@meta.data$orig.ident
            # 得到这些细胞的exp
            data <- data.frame(exp = exp,cluster = res,Sample = Sample)
          }
          
          p <- plot_ly() %>%
            add_trace(
              data = data,
              x = ~cluster[ data$Sample == samples[1] ],
              y = ~exp[ data$Sample == samples[1] ],
              legendgroup = samples[1],
              scalegroup = samples[1],
              name = samples[1],
              #split = ~cluster[ data$Sample == samples[1] ],
              type = 'violin',
              side = 'negative',
              color = I("#8dd3c7"),
              showlegend = showlegend,
              box = list(visible = T),
              meanline = list(visible = T)
            ) %>%
            add_trace(
              data = data,
              x = ~cluster[ data$Sample == samples[2] ],
              y = ~exp[ data$Sample == samples[2] ],
              legendgroup = samples[2],
              scalegroup = samples[2],
              name = samples[2],
              #split = ~cluster[ data$Sample == samples[2] ],
              type = 'violin',
              side = 'positive',
              color = I("#bebada"),
              showlegend = showlegend,
              box = list(visible = T),
              meanline = list(visible = T)
            ) %>%
            layout(
              annotations = anno_gen(feature),
              xaxis = list(title = 'Cluster', tickangle = -45),
              yaxis = list(title = 'Expression Level',zeroline = T)
            )
          all_G[[ (length(all_G)+1) ]] <- p
        }
      }
    }else{
      # 根本就不是多样本
      all_G <- list()
      for (i in 1:length(pick_features)) {
        feature <- pick_features[i]
        data <- NULL
        if(input$categoryforviolin == 'Seurat cluster'){
          data <- data.frame(
            exp = values$seurat_obj@assays[[DefaultAssay(values$seurat_obj)]]@data[feature,],
            cluster = as.character(Idents(values$seurat_obj)),
            Sample = values$seurat_obj@meta.data$orig.ident)
          data <- data[which(as.character(data$cluster) %in% input$clusterforviolin),]
          #data$cluster <- factor(as.character(data$cluster),levels = levels(Idents(values$seurat_obj))[which(levels(Idents(values$seurat_obj)) %in% as.character(data$cluster))])
          data$cluster <- factor(as.character(data$cluster),levels = input$clusterforviolin)
        }else{
          cluster_list <- values$def_cate_clu[[input$categoryforviolin]]
          cells <- as.character(unlist(cluster_list))
          res <- NULL
          for (i in 1:length(cluster_list)) {
            res <- c(res,rep(names(cluster_list[i]),length(cluster_list[[i]])))
          }
          
          temp_obj <- values$seurat_obj[,cells]
          exp <- temp_obj@assays[[DefaultAssay(values$seurat_obj)]]@data[feature,]
          # 得到这些细胞来自哪个sample
          Sample <- temp_obj@meta.data$orig.ident
          # 得到这些细胞的exp
          data <- data.frame(
            exp = exp,
            cluster = res,
            Sample = Sample
          )
        }
        
        p <- plot_ly() %>%
          add_trace(
            data = data,
            x = ~cluster,
            y = ~exp,
            split = ~cluster,
            type = 'violin',
            color = ~cluster,
            colors = colours_sample,
            box = list(visible = T),
            meanline = list(visible = T)
          ) %>%
          layout(
            annotations = anno_gen(feature),
            xaxis = list(title = 'Cluster',tickangle = -45),
            yaxis = list(title = 'Expression Level',zeroline = T)
          ) %>% hide_legend()
        all_G[[ (length(all_G)+1) ]] <- p
      }
    }
    if(length(all_G) > 1){
      subplot( all_G, nrows = length(pick_features), shareY = T, shareX = T, titleX = F, titleY = F) %>% layout(
        title = list(y = 0.1),
        xaxis = list(title = 'Cluster'),
        yaxis = list(title = 'Expression Level')
      )%>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
        config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat_stat_v, height = input$plotheight_stat_v, width = input$plotwidth_stat_v,filename=paste0("violin")))
    }else if(length(all_G) == 1){
      subplot( all_G[[1]], nrows = length(pick_features), shareY = T, shareX = T, titleX = F, titleY = F)%>% layout(
        title = list(y = 0.1),
        xaxis = list(title = 'Cluster'),
        yaxis = list(title = 'Expression Level')
      )%>% onRender("function(el,x){el.on('plotly_legendclick', function(){ return false; })}") %>%
        config(modeBarButtonsToRemove = c("autoScale2d","zoomIn2d","zoomOut2d","hoverClosestCartesian", "hoverCompareCartesian"),displaylogo = FALSE,toImageButtonOptions = list(format = input$plotformat_stat_v, height = input$plotheight_stat_v, width = input$plotwidth_stat_v,filename=paste0("violin")))
    }
  })
  violin_output_obj2 <- reactive({NULL})
  
  output$all_violin_plot <- renderPlotly({
    if(!is.null(input$violin_map_run) && input$violin_map_run != 0 ){
      violin_output()
    }else{
      violin_output_obj2()
    }})
  
  output$all_violin_plot_ui <- renderUI({
    if(!is.null(input$violin_map_run) && input$violin_map_run != 0 ){
      fluidRow(
        shinycssloaders::withSpinner(
          plotlyOutput('all_violin_plot',height =  paste0(values$violin_height*200,'px'))
        )
      )
      }else{
        fluidRow()
    }
  })
  
  
  output$selecttopormygenelist_for_dotplot <- renderUI({
    options <- NULL
    if(length(values$def_gene_list[[pick_gene_list()]]) == 0 || values$def_gene_list[[pick_gene_list()]] == ''){
      options <- 'Top gene'
    }else{
      options <- c( 'Top gene','My gene list' )
    }
    radioGroupButtons(
      inputId = "selecttopormygenelist_for_dotplot",
      label = "Features come from:",
      choices = options,
      justified = TRUE,
      selected = 'Top gene',
      checkIcon = list(yes = icon("ok", lib = "glyphicon"))
    )
  })
  
  output$selectdotplotfeature <- renderUI({
    if(input$selecttopormygenelist_for_dotplot == "My gene list"){
      selectInput(
        inputId = 'input_selectdotplotfeature',
        label = "DotPlot Features",
        choices = values$def_gene_list[[pick_gene_list()]],
        selected = values$def_gene_list[[pick_gene_list()]][1],
        multiple = T
      )
    }else if(input$selecttopormygenelist_for_dotplot == 'Top gene'){
      selectInput(inputId = "topnfordotplot",
                  label = "Top n for DotPlot:",
                  width = '100%',
                  choices = 1:15,
                  selectize = T)
    }else{
      return(NULL)
    }
  })
  
  output$split_sample_check_box_for_dotplot <- renderUI({
    if("integrated" %in% names(values$seurat_obj@assays) ){
      awesomeCheckbox(
        inputId = "if_split_sample_for_dotplot",
        label = "Split sample", 
        value = F
      )
    }else{
      return(NULL)
    }
  })
  
  output$categoryfordotplot <- renderUI({
    options <- NULL
    if( length(names(values$def_cate_clu)) == 0 ){
      options <- 'Seurat cluster'
    }else{
      if(input$input_ifusemycategory){
        options <- c('Seurat cluster',names(values$def_cate_clu))
      }else{
        options <- 'Seurat cluster'
      }
    }
    radioGroupButtons(
      inputId = "categoryfordotplot",
      label = "Cells come from:",
      choices = options,
      justified = TRUE,
      selected = 'Seurat cluster',
      checkIcon = list(
        yes = icon("ok", lib = "glyphicon"))
    )
  })
  
  output$clusterfordotplot <- renderUI({
    options <- NA
    selected <- NA
    if(input$categoryfordotplot == 'Seurat cluster'){
      options <- as.character(levels(Idents(values$seurat_obj)))
      selected <- as.character(levels(Idents(values$seurat_obj)))
    }else{
      options <- names(values$def_cate_clu[[input$categoryfordotplot]])
      selected <- names(values$def_cate_clu[[input$categoryfordotplot]])
    }
    selectInput(inputId = "clusterfordotplot",
                label = "Cluster for DotPlot:",
                width = '100%',
                choices = options,
                selected = selected,
                multiple = T,
                selectize = T)
  })
  
  output$pickcolorfordotplot <- renderUI({
    selectInput(
      inputId = "selectcolorfordotplot",
      label = "Select color for DotPlot",
      choices = brewer.pal(9,'Set1'),
      selected = NULL,
      selectize = TRUE,
      width = '100%'    
    )
  })
  
  output$plotmodel_dotplot <- renderUI({
    fluidRow(
      column(width = 3, sliderInput("plotheight_stat_dot", "Plot output height (px)",min = 0, max = 2000, value = 500, step = 100)),
      column(width = 3, sliderInput("plotwidth_stat_dot", "Plot output width (px)",min = 0, max = 2000, value = 500, step = 100)),
      column(width = 3, radioGroupButtons(inputId = "plotformat_stat_dot",label = "Plot format",choices = c("svg", "png"),
                                          checkIcon = list(yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                           no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))),
      column(width = 3, downloadButton(
        outputId = "save_dotplot",
        label = "Save DotPlot",
        style = "color: white; background-color: #1558A0;border-color: #1558A0",
        width = "100%"
      ))
    )
  })
  
  output$save_dotplot <- downloadHandler(
    filename = function(){
      paste0('scui_DotPlot.',input$plotformat_stat_dot)
    },
    content = function(con){
      if(input$plotformat_stat_dot == 'svg'){
        svg(file = con,width = input$plotwidth_stat_dot,height = input$plotheight_stat_dot)
        print(dotplot_output())
        dev.off()
      }else if(input$plotformat_stat_dot == 'png'){
        png(file = con,width = input$plotwidth_stat_dot,height = input$plotheight_stat_dot)
        print(dotplot_output())
        dev.off()
      }
    }
  )
  
  dotplot_output <- eventReactive(input$dot_plot_run,{
    pick_features <- NULL
    if(input$selecttopormygenelist_for_dotplot == "My gene list"){
      pick_features <- input$input_selectdotplotfeature
    }else{ # top gene
      if( input$categoryfordotplot == 'Seurat cluster' ){ # 原来类的top基因
        pick_features <- values$marker_res_save[['ori_findmarkers_res']] %>% group_by(cluster) %>% top_n(n = as.numeric(input$topnfordotplot), wt = avg_logFC)
        pick_features <- pick_features$gene
      }else{ # 自定义类的top基因
        pick_features <- values$marker_res_save[['findmarkers_res']] %>% group_by(marker_group) %>% top_n(n = as.numeric(input$topnfordotplot), wt = avg_logFC)
        pick_features <- pick_features$gene
      }
    }
    
    values$dotplot_height <- length(levels(Idents(values$seurat_obj)))
    if( ! is.null(input$if_split_sample_for_dotplot) && input$if_split_sample_for_dotplot == T){
      values$dotplot_height <- values$dotplot_height * 1.5
    }
    
    values$seurat_obj@meta.data$seurat_clusters <- Idents(values$seurat_obj)
    group.by <- 'seurat_clusters'
    split.by <- NULL
    
    if( ! is.null(input$if_split_sample_for_dotplot) && input$if_split_sample_for_dotplot == T){
      split.by <- 'orig.ident'
    }
    
    if(input$categoryfordotplot == 'Seurat cluster'){ # 按照原版的分类
      p <- DotPlot(values$seurat_obj, 
                   assay = DefaultAssay(values$seurat_obj),
                   features = unique(pick_features), 
                   cols = c("lightgrey", input$selectcolorfordotplot),#split.by = "orig.ident",
                   group.by = group.by,
                   split.by = split.by,
                   idents = input$clusterfordotplot,
                   dot.scale = 4) + RotatedAxis()
      return(p)
    }else{ # 按照自定义的分类
      # 1按自定义的细胞类型取子集
      # 2为新的子集设置my——def
      # groupby mydef
      # 得到cluster列表
      cluster_list <- values$def_cate_clu[[input$categoryfordotplot]]
      # 得到所有选择的细胞
      cells <- as.character(unlist(cluster_list))
      # 得到这些细胞的类别
      res <- NULL
      for (i in 1:length(cluster_list)) {
        res <- c(res,rep(names(cluster_list[i]),length(cluster_list[[i]])))
      }
      temp_obj <- values$seurat_obj[,cells]
      
      temp_obj@meta.data$my_clusters <- res
      
      p <- DotPlot(temp_obj, 
                   assay = DefaultAssay(values$seurat_obj),
                   features = unique(pick_features), 
                   cols = c("lightgrey", input$selectcolorfordotplot),#split.by = "orig.ident",
                   group.by = 'my_clusters',
                   split.by = split.by,
                   idents = input$clusterfordotplot,
                   dot.scale = 4) + RotatedAxis()
      return(p)
    }
  })
  
  output$dotplot <- renderPlot({
    # ggplotly(dotplot_output())
    dotplot_output()
  })
  
  output$all_dot_plot_ui <- renderUI({
    fact <- 250
    if(values$dotplot_height < 10){
      fact <- 250
    }else{
      fact <- values$dotplot_height*25
    }
    fluidRow(
      shinycssloaders::withSpinner(
        plotOutput('dotplot',height = paste0(fact,'px'))
      )
    )
  })
  
  # ------------------------------------------- #
  output$pick_category_for_recluster <- renderUI({
    all_def_cateory <- names(values$def_cate_clu)
    if(length(all_def_cateory) >= 1){
      selectInput(
        inputId = 'pickcategoryforrecluster',
        label = "Category for Re-clustering",
        choices = all_def_cateory,
        selected = all_def_cateory[1],
        multiple = FALSE,
        width = NULL,
      ) 
    }
  })
  output$def_category_msg <- renderUI({
    pick_cell <- unlist(values$def_cate_clu[[input$pickcategoryforrecluster]])
    msg <- sprintf('<h3>There are <font color="#1558A0">%s</font> cells in category <font color="#1558A0">%s</font> defined by you.</h3>',length(pick_cell),input$pickcategoryforrecluster)
    h3(HTML(msg))
  })
  
  running_msg_obj1 <- reactive({
    fluidRow(
      h3('Running Seurat cluster process...')
    )
  })
  running_msg_obj2 <- reactive({
    fluidRow(
      h3('Process Complete!')
    )
  })
  running_msg_obj3 <- reactive({
    fluidRow(
      h3('Please select more then 100 cell for re-clustering!')
    )
  })
  
  observeEvent(input$recluster_run,{
    pick_cell <- unlist(values$def_cate_clu[[input$pickcategoryforrecluster]])
    if(length(pick_cell) > 100){
      sub_obj <- values$seurat_obj
      sub_obj <- values$seurat_obj[,pick_cell]
      sub_obj <- SCTransform(sub_obj, verbose = T)
      sub_obj <- RunPCA(sub_obj, assay ='SCT', verbose = FALSE)
      sub_obj <- RunUMAP(sub_obj, assay ='SCT', dims = 1:30, verbose = FALSE)
      sub_obj <- FindNeighbors(sub_obj,assay ='SCT', dims = 1:30, verbose = FALSE)
      sub_obj <- FindClusters(sub_obj,assay ='SCT', verbose = FALSE)
      sub_obj_markers <- FindAllMarkers(
        sub_obj,
        assay = 'SCT',
        only.pos = T,
        features = VariableFeatures(sub_obj)[1:2000]
      )
      values$new_obj = sub_obj
      values$new_obj_markers = sub_obj_markers
    }else{
    }
  })
  
  # observeEvent(input$recluster_run,{
  #   pick_cell <- unlist(values$def_cate_clu[[isolate({input$pickcategoryforrecluster})]])
  #   if(length(pick_cell) > 100){
  #     output$running_msg <- renderUI({running_msg_obj2()})
  #   }else{
  #     output$running_msg <- renderUI({running_msg_obj3()})
  #   }
  # })
  
  output$reclustermsg <- renderPrint({
    msg <- ''
    if(input$recluster_run){
      pick_cell <- unlist(values$def_cate_clu[[input$pickcategoryforrecluster]])
      if(length(pick_cell) > 100){
        msg <- 'Running... may take several minutes..'
      }else{
        msg <- 'Please select more then 100 cell for re-clustering!'
      }
    }else{
      msg <- "Click \'Run!\' to start. It may takes several minutes.."
    }
    if(is.null(values$new_obj)){
      print(msg)
    }else{
      print(values$new_obj)
    }
  })

  recluster_result_obj1 <- reactive({
    fluidRow(
      downloadButton(
        outputId = "downloadnewobj",
        label = "Download New Seurat Object!",
        style = "color: white; background-color: #1558A0;border-color: #1558A0",
        width = "100%"
      )
    )
  })
  recluster_result_obj2 <- reactive({fluidRow()})

  observe({
    if( !is.null(values$new_obj) ){
      output$recluster_result <- renderUI({
        recluster_result_obj1()
        })
    }else{
      output$recluster_result <- renderUI({
        recluster_result_obj2()
        })
    }
  })
  
  output$downloadnewobj <- downloadHandler(
    filename = function(){
      'scui_new_obj.Rdata'
    },
    content = function(con){
      seurat_obj <- values$new_obj
      seurat_obj.markers <- values$new_obj_markers
      save(
        list = as.character(c(
          'seurat_obj',
          'seurat_obj.markers'
        )),
        file = con
      )
    }
  )
}








