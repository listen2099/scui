# 读取seurat和marker
# 单样本和混合分情况

# combined <- RunTSNE(combined)
# save(combined,combined.markers,file = 'C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/raw_combined.markers.Rdata')

# load('C:/BGIworkdir/BGIDOC/monocle2/Cluster_Biomarkers.Rdata')

load('C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/raw_combined.markers.Rdata')
set.seed(1)
pickcell <- sample(colnames(combined),5000)
combined <- combined[,pickcell]
save(combined,combined.markers,file = 'C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/raw_combined.markers_2.Rdata')

# ------------准备数据------------- #
load('C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/raw_combined.markers_1.Rdata')
# load('C:/Users/weikun1/Downloads/scui_project_obj2.Rdata')
# load('C:/Users/weikun1/Downloads/scui_new_obj.Rdata')

#load('C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/Cluster_Biomarkers.Rdata')
read_env <- ls(environment())

seurat_obj <- NULL
seurat_obj_name <- NULL

for (obj in read_env) {
  if(class(environment()[[obj]]) == "Seurat"){
    seurat_obj <- environment()[[obj]]
    seurat_obj_name <- obj
  }
}

#table_cell <- seurat_obj@meta.data
seurat_obj_markers <- environment()[[ paste0(seurat_obj_name,'.markers') ]]

res <- list(
  seurat_obj = seurat_obj,
  seurat_obj_markers = seurat_obj_markers
)

res$seurat_obj

library(Seurat)
set.seed(1)
pick_id <- sample(colnames(res$seurat_obj),100)


sub_obj <- res$seurat_obj[,pick_id]
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

DimPlot(sub_obj, label = TRUE) + NoLegend()
print(DimPlot(sub_obj, label = TRUE))



sub_obj@assays$RNA <- sub_obj@assays$SCT
DefaultAssay(sub_obj) <- 'RNA'
sub_obj@assays$SCT <- NULL


library(dplyr)
res$seurat_obj
marker <- res$seurat_obj_markers %>% group_by(cluster) %>% top_n(2,avg_logFC)
marker$gene

DefaultAssay(res$seurat_obj)
p <- DotPlot(res$seurat_obj,
             assay = "RNA",
             group.by = 'seurat_clusters',
        features = c("Xkr4","Tcea1"),# unique(marker$gene), 
        cols = c("lightgrey", "blue"),split.by = "stim",
        dot.scale = 4) + RotatedAxis()
p
ggplotly(p)

FindAllMarkers(seurat_obj)

seurat_obj@meta.data[['marker_group']] 

msg <- print(a <- FindMarkers(seurat_obj,ident.1 = "control", group.by = 'stim'))

# ------------准备数据------------- #


# ------------画点图---------- #
library(Seurat)
library(ggplot2)
library(plotly)

plot.data <- as.data.frame(cbind(seurat_obj@reductions$umap@cell.embeddings[,1:2],
                   Cluster = Idents(seurat_obj)))
plot.data$Sample = seurat_obj@meta.data$orig.ident

p <- plot_ly(plot.data,x = ~UMAP_1,y = ~UMAP_2) %>% 
  add_markers(color = ~factor(Cluster),
              symbol = ~factor(Sample), 
              symbols = c("circle","square"),
              type = "scatter", mode = "markers") %>% 
  layout(xaxis = list(showgrid=F,zeroline=F,showticklabels=F),
         yaxis = list(showgrid=F,zeroline=F,showticklabels=F)) 

pdf(file = 'C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/test.pdf')
export(p, file = "C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/test.pdf")
export(p, file = "C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/test.svg")
orca(p, file = "C:/BGIworkdir/BGIDOC/Single cell cell type study/SingleR/test/test.pdf")
toWebGL(p)

plotly_json(p)
schema(p)


plotly_example("shiny", "event_data",edit = TRUE)

num <- length(unique( Idents(seurat_obj) ))
colours_sample<-colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFE528","#A65628","#F781BF","#999999"))(num)[num:1]
names(colours_sample)<-levels(seurat_obj)
# ------------画点图---------- #

# ------------画点图，cluster单独对象---------- #
library(Seurat)
library(ggplot2)
library(plotly)

num <- length(unique( Idents(seurat_obj) ))
colours_sample<-colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFE528","#A65628","#F781BF","#999999"))(num)[num:1]
names(colours_sample) <- levels(Idents(seurat_obj))

plot.data <- as.data.frame(cbind(seurat_obj@reductions$umap@cell.embeddings[,1:2]),)
colnames(plot.data) <- c('D1','D2')
plot.data.list <- list()
p <- plot_ly()
for (i in 1:length(levels(Idents(seurat_obj)))) {
  one_cluster <- levels(Idents(seurat_obj))[i]
  plot.data.list[[ one_cluster ]] <- plot.data[which(Idents(seurat_obj) == one_cluster),]
  p <- add_trace(p,data = plot.data.list[[i]], type = "scatter", mode = "markers",
                 size = I(2.5),
                 color = I(colours_sample[i]),
                 x = ~D1,y = ~D2)#,name = one_cluster)
}
p %>% layout(showlegend = FALSE)
p %>% hide_legend()

# 创建多个对应的check box
library(shiny)
library(plotly)

hoverinfo='skip'

ui <- fluidPage(
  sliderInput("marker", "Marker size", min = 0, max = 20, value = 8),
  sliderInput("path", "Path size", min = 0, max = 30, value = 2),
  plotlyOutput("p")
)

server <- function(input, output, session) {
  
  output$p <- renderPlotly({
    plot_ly(
      economics, x = ~pce, y = ~psavert, z = ~unemploy, 
      color = ~as.numeric(date), mode = "markers+lines"
    )
  })
  
  observeEvent(input$marker, {
    plotlyProxy("p", session) %>%
      plotlyProxyInvoke(
        "restyle", 
        # could also do list(marker = list(size = input$marker))
        # but that overwrites the existing marker definition
        # https://github.com/plotly/plotly.js/issues/1866#issuecomment-314115744
        list(marker.size = input$marker)
      )
  })
  
  observeEvent(input$path, {
    plotlyProxy("p", session) %>%
      plotlyProxyInvoke(
        "restyle", list(line.width = input$path)
      )
  })
  
}

shinyApp(ui, server)




# ------------画点图，cluster单独对象---------- #


# --------小提琴图-------- #
library(Seurat)
library(ggplot2)
library(plotly)
query <- seurat_obj_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) %>% select(gene)
vln <- VlnPlot(seurat_obj, features = query$gene[1:2],
               combine = T,pt.size=0) + NoLegend() + theme(text = element_text(size = 10),axis.text = element_text(size = 8))
p <-plot_ly()
subplot(
  nrows = 2,
  ggplotly(vln[[1]]) %>% layout(annotations = list(text='A')),
  ggplotly(vln[[2]]) %>% layout(annotations = list(text='B')),
  shareX = F,
  shareY = F,
  titleX = T,
  titleY = T
)

toWebGL(ggplotly(vln))
p <- VlnPlot(seurat_obj, features = query$gene[1:2],
        combine = F,pt.size=0)
subplot(p[[1]],p[[2]])
subplot(vln[[1]],vln[[2]])
# --------小提琴图-------- #

# ------plot_ly小提琴，实现多样本比较------- #
# 小提琴的本质：
# 某个feature在不同cluster内的表达量的统计
# 左半边是sampleA，右半边是sampleB
library(plotly)

df <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv")
head(df)

num <- length(unique( Idents(seurat_obj) ))
colours_sample<-colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFE528","#A65628","#F781BF","#999999"))(num)[num:1]
names(colours_sample) <- as.character(levels(Idents(seurat_obj)))

df <- df[which(df$day %in% c('Fri','Thur')),]


plot_ly() %>%
  add_trace(
    data = df,
    x = ~day,
    y = ~total_bill,
    split = ~day,
    type = 'violin',
    colors = colours_sample,
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
) %>% layout(
  xaxis = list(title = "Day"),
  yaxis = list(title = "Total Bill",zeroline = F)
) 

# ------plot_ly小提琴，实现多样本比较------- #
query <- seurat_obj_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) %>% select(gene)
feature <- query$gene[1]
data <- data.frame(
  exp = seurat_obj@assays$RNA@data[feature,],
  cluster = Idents(seurat_obj),
  Sample = seurat_obj@meta.data$orig.ident
)

f <- list(
  #family = "Courier New, monospace",
  size = 14,
  color = "black")

a <- list(
  text = "SUBPLOT TITLE A",
  font = f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "center",
  x = 0.5,
  y = 1,
  showarrow = FALSE
)

b <- list(
  text = "SUBPLOT TITLE B",
  font = f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "center",
  x = 0.5,
  y = 1,
  showarrow = FALSE
)

p1 <- plot_ly() %>%
  add_trace(
    data = data,
    x = ~cluster[ data$Sample == 'AD' ],
    y = ~exp[ data$Sample == 'AD' ],
    split = ~cluster[ data$Sample == 'AD' ],
    type = 'violin',
    side = 'negative',
    color = ~cluster[ data$Sample == 'AD' ],
    colors = colours_sample,
    box = list(visible = T),
    meanline = list(visible = T)
  ) %>% 
  add_trace(
    data = data,
    x = ~cluster[ data$Sample == 'control' ],
    y = ~exp[ data$Sample == 'control' ],
    split = ~cluster[ data$Sample == 'control' ],
    type = 'violin',
    side = 'positive',
    color = ~cluster[ data$Sample == 'control' ],
    colors = colours_sample,
    box = list(visible = T),
    meanline = list(visible = T)
  ) %>%
  layout(
    annotations = a,
    xaxis = list(title = 'Cluster'),
    yaxis = list(title = 'Expression Level',zeroline = T)
  ) %>% hide_legend()

query <- seurat_obj_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) %>% select(gene)
feature <- query$gene[2]
data <- data.frame(
  exp = seurat_obj@assays$RNA@data[feature,],
  cluster = Idents(seurat_obj),
  Sample = seurat_obj@meta.data$orig.ident
)

p2 <- plot_ly() %>%
  add_trace(
    data = data,
    x = ~cluster[ data$Sample == 'AD' ],
    y = ~exp[ data$Sample == 'AD' ],
    split = ~cluster[ data$Sample == 'AD' ],
    type = 'violin',
    side = 'negative',
    color = ~cluster[ data$Sample == 'AD' ],
    colors = colours_sample,
    box = list(visible = T),
    meanline = list(visible = T)
  ) %>% 
  add_trace(
    data = data,
    x = ~cluster[ data$Sample == 'control' ],
    y = ~exp[ data$Sample == 'control' ],
    split = ~cluster[ data$Sample == 'control' ],
    type = 'violin',
    side = 'positive',
    color = ~cluster[ data$Sample == 'control' ],
    colors = colours_sample,
    box = list(visible = T),
    meanline = list(visible = T)
  ) %>%
  layout(
    annotations = b,
    xaxis = list(title = 'Cluster'),
    yaxis = list(title = 'Expression Level',zeroline = T)
  ) %>% hide_legend()

subplot( list(p1,p2), nrows = 2,shareY = T ,shareX =T,
         titleX = F, titleY = F
         ) %>% layout(
           xaxis = list(title = 'Cluster'),
           yaxis = list(title = 'Expression Level')
         ) 

# ----------绘制基因表达量------------ #
library(dplyr)
gene_test_list <- seurat_obj_markers %>% group_by(cluster) %>% top_n(1,desc(avg_logFC )) %>% select(gene)
pick_gene <- gene_test_list$gene[1]


plot.data <- as.data.frame(cbind(seurat_obj@reductions$umap@cell.embeddings[,1:2]))
plot.data$Cluster <- Idents(seurat_obj)

plot.data[['pick_gene']] <- as.numeric(seurat_obj@assays$RNA[pick_gene,])

p <- plot_ly(plot.data,x = ~UMAP_1,y = ~UMAP_2) %>% 
  add_markers(color = ~pick_gene,
              colors = colorRamp(c("#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b" )),
              type = "scatter", mode = "markers") %>% 
  colorbar(title = pick_gene) %>%
  hide_legend() %>%
  layout(xaxis = list(showgrid=F,zeroline=F,showticklabels=F),
         yaxis = list(showgrid=F,zeroline=F,showticklabels=F)) 

anno_gen <- function(title){
  f <- list(#family = "Courier New, monospace",
    size = 14,color = "black")
  a <- list(text = title,font = f,xref = "paper",yref = "paper",
            yanchor = "bottom",xanchor = "center",align = "center",
            x = 0.5,y = 1,showarrow = FALSE)
  return(a)
}

all_G <- list()
all_cluster <- levels(Idents(seurat_obj))
for (i in 1:length(all_cluster)) {
  temp <- plot.data[which(plot.data$Cluster == all_cluster[i]),]
  max_y<-ceiling(max(temp$pick_gene)+max(temp$pick_gene)*0.1)
  min_y<-floor(min(temp$pick_gene)+max(temp$pick_gene)*0.1)
  p <- plot_ly(temp,x = ~UMAP_1,y = ~UMAP_2) %>% 
    add_markers(color = ~pick_gene,
                colors = colorRamp(c("#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b" )),
                coloraxis = 'coloraxis',
                type = "scatter", mode = "markers") %>% 
    colorbar(title = pick_gene,limits = c(min(plot.data$pick_gene),max(plot.data$pick_gene))) %>% 
    hide_legend() %>% hide_colorbar() %>%
    layout(
      annotations = anno_gen(all_cluster[i]),
      xaxis = list(showgrid=F,zeroline=F,showticklabels=F,
                        showline = TRUE,
                        mirror = "ticks",
                        linecolor = toRGB("black"),
                        linewidth = 2,
                        range = c(min_y, max_y)),
           yaxis = list(showgrid=F,zeroline=F,showticklabels=F,
                        showline = TRUE,
                        mirror = "ticks",
                        linecolor = toRGB("black"),
                        linewidth = 2,
                        range = c(-1, 10)))
  all_G[[ length(all_G) + 1 ]] <- p
}

subplot(all_G,
        nrows = round(sqrt(length(all_G)))
        ) %>% layout(coloraxis=list(colorscale='Jet'))


# ----------热图绘制2------------ #
# 热图有3种方法：
# 1.seurat自己的热图
library(Seurat)
library(dplyr)
# feature 可以选择my list，还可以选择top几
#  |-首先选择topn还是my list
#     |-如果是topn就选top几
#     |-如果是mylist直接画
# cell 可以选择哪些cluster
#  |-首先选择category，再选择这个category里的某些类
# run触发画图的操作
top10 <- seurat_obj_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
p <- DoHeatmap(seurat_obj, features = top10$gene,
               label = F)
p <- p + theme(
  axis.text.y = element_text(size = 5)
)
ggplotly(p)

# 2.我们设计的seurat热图
library(ComplexHeatmap)
library(ggplotify)
DoHeatmapPlot_num <- 20 #展示所有cluster中最显著10个基因的DoHeatmapPlot
DoHeatmapPlot <- NA #NA表示不自定义基因'char'
seurat_obj <- seurat_obj
markers_obj <- seurat_obj_markers
fontsize <- 5
DoHeatMap <- function(seurat_obj,markers_obj,DoHeatmapPlot_num,DoHeatmapPlot,fontsize=5){
  seurat_cluster <- Idents(seurat_obj)
  
  num = length(levels(seurat_cluster))
  colours_sample=colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFE528","#A65628","#F781BF","#999999"))(num)[num:1]
  names(colours_sample) <- levels(seurat_cluster)
  scale_mat <- GetAssayData(seurat_obj, slot = "scale.data")
  
  top_plot <- markers_obj %>% group_by(cluster) %>% top_n(n = DoHeatmapPlot_num, wt = avg_logFC)
  mat <- as.matrix(scale_mat[top_plot$gene, ])
  
  ha = HeatmapAnnotation(
    Cluster = as.character(seurat_cluster),
    col = list(Cluster = colours_sample),
    show_legend = F
  )
  
  la = rowAnnotation(
    Cluster = top_plot$cluster,
    col = list(Cluster = colours_sample)
  )
  
  if(is.na(DoHeatmapPlot)){# 没有指定
    if_show_row_names <- T
    if_right_annotation <- NULL
  }else if(is.character(DoHeatmapPlot)){ # 指定了
    DoHeatmap_gene = strsplit(DoHeatmapPlot,',')[[1]]
    if_show_row_names <- F
    pick_res <- markers_obj %>% filter( gene %in% DoHeatmap_gene )
    s1 <- paste0(pick_res$gene,'-',pick_res$cluster)
    s2 <- paste0(top_plot$gene,'-',top_plot$cluster)
    gene_pos <- which( s2 %in% s1 )
    if_right_annotation = rowAnnotation(
      mark_gene = anno_mark(at = gene_pos, labels = rownames(mat)[gene_pos])
    )
  }
  ran <- quantile(as.numeric(mat), probs = c(0.01,0.99))
  Heatmap(mat,
          cluster_rows = FALSE,
          show_column_names = FALSE,
          show_row_names = if_show_row_names,
          cluster_columns = cluster_within_group(mat, seurat_cluster),
          column_dend_side = 'top',
          column_dend_reorder = F,
          column_split = length(levels(seurat_cluster)),
          column_title = NULL,
          top_annotation = ha,
          left_annotation = la,
          right_annotation = if_right_annotation,
          heatmap_legend_param = list(title = 'Gene Exp'),
          row_names_gp = gpar(fontsize = fontsize),
          col = circlize::colorRamp2(c(ran[1],
                                       mean( as.numeric(mat)[which(as.numeric(mat) > 0 & as.numeric(mat)<=ran[2])] ),
                                       ran[2]),
                                     c("#4575b4","#f7f7f7","#d73027"))
  )
}
pdf(file = 'C:/BGIworkdir/DoHeatMap/test3.pdf')
p <- DoHeatMap(seurat_obj,seurat_obj_markers,
          DoHeatmapPlot_num,DoHeatmapPlot
          )
dev.off()
p2 <- as.ggplot(p)
p3 <- as.grob(p)
class(p3)
ggplotly(p3)

# 3.模仿loupe的热图
library(dplyr)
library(ggplot2)
library(Seurat)
library(heatmaply)
top <- seurat_obj_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top$gene

seurat_obj@assays$RNA@data



heatmap <- DoHeatmap(seurat_obj, features = top$gene) + NoLegend() + theme(text = element_text(size = 10),axis.text = element_text(size = 8))
hdata <- heatmap$data
rowN <- names(table(hdata$Identity))
colN <- names(table(hdata$Feature))
res_met <- matrix(
  NA,
  nrow = length(rowN),
  ncol = length(colN)
)
rownames(res_met) = rowN
colnames(res_met) = colN
for (i in 1:length(colN) ) {
  gene <- colN[i]
  tdata <- hdata[which( hdata$Feature == gene),c(3,4)]
  res_met[,i] <- aggregate(tdata[,1],list(tdata$Identity),mean,na.rm = T  )[,2]
}
res_met <- t(scale( t(res_met)))

set.seed(1)
e.dist <- dist(t(res_met),method = "euclidean")
h.cluster <- hclust(e.dist, method="ward.D2")
cut.h.cluster <- cutree(h.cluster, k=13)
cut.h.cluster <- sort(cut.h.cluster)
res_met <- res_met[,names(cut.h.cluster)]

e.dist <- dist((res_met),method = "euclidean")
h.cluster <- hclust(e.dist, method="ward.D2")
cut.h.cluster <- cutree(h.cluster, k=4)
cut.h.cluster <- sort(cut.h.cluster)
res_met <- res_met[names(cut.h.cluster),]

p <- ggplot()
p <- p + geom_tile(data = reshape2::melt(res_met), aes(Var2, Var1, fill = value))
p <- p + scale_fill_gradient2(name = 'Gene Exp',high = 'red',low = 'blue',mid = "white")
p <- p + scale_y_continuous(breaks = 0:(length(rowN)-1),labels = paste0('cluster ',0:(length(rowN)-1)))
p <- p + theme_minimal()
p <- p + theme(
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  panel.grid = element_blank()
)
ggplotly(p)

# -------plot_ly热图测试---------- # 
library(dplyr)
library(ggplot2)
library(Seurat)
library(plotly)
library(heatmaply)
top <- seurat_obj_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
features <- top$gene

data <- seurat_obj@assays$RNA@data

mat <- matrix(0,
  nrow = length(features),
  ncol = length(levels(Idents(seurat_obj)))
)

rownames(mat) <- features
colnames(mat) <- levels(Idents(seurat_obj))

for (i in 1:length(mat[1,])) {
  one_cluster <- levels(Idents(seurat_obj))[i]
  one_cluser_data <- data[features,Idents(seurat_obj) == one_cluster]
  mat[,i] <- rowMeans(as.matrix(one_cluser_data))
}

heatmaply(mat,Rowv=NULL,fontsize_row = 5)

#基因是横，cluster是列



# ----------热图绘制2------------ #


# ------------------------------ #
library(shiny)
library(plotly)

ui <- fluidPage(
  radioButtons("plotType", "Plot Type:", choices = c("ggplotly", "plotly")),
  plotlyOutput("plot"),
  # verbatimTextOutput("hover"),
  verbatimTextOutput("click"),
  # verbatimTextOutput("brushing"),
  # verbatimTextOutput("selecting"),
  # verbatimTextOutput("brushed"),
  verbatimTextOutput("selected")
)

server <- function(input, output, session) {
  
  nms <- row.names(mtcars)
  
  output$plot <- renderPlotly({
    p <- if (identical(input$plotType, "ggplotly")) {
      ggplotly(ggplot(mtcars, aes(x = mpg, y = wt, customdata = nms)) + geom_point())
    } else {
      plot_ly(mtcars, x = ~mpg, y = ~wt, customdata = nms)
    }
    p %>% 
      layout(dragmode = "select") %>%
      event_register("plotly_selecting")
  })
  
  # output$hover <- renderPrint({
  #   d <- event_data("plotly_hover")
  #   if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  # })
  # 
  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) "Click events appear here (double-click to clear)" else d
  })
  
  # output$brushing <- renderPrint({
  #   d <- event_data("plotly_brushing")
  #   if (is.null(d)) "Brush extents appear here (double-click to clear)" else d
  # })
  # 
  # output$selecting <- renderPrint({
  #   d <- event_data("plotly_selecting")
  #   if (is.null(d)) "Brush points appear here (double-click to clear)" else d
  # })
  # 
  # output$brushed <- renderPrint({
  #   d <- event_data("plotly_brushed")
  #   if (is.null(d)) "Brush extents appear here (double-click to clear)" else d
  # })
  # 
  output$selected <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Brushed points appear here (double-click to clear)" else d
  })
  
}

shinyApp(ui, server, options = list(display.mode = "showcase"))

# 16 ------------------------------------ # 
library(shiny)
library(plotly)

# Generate 100,000 observations from 2 correlated random variables
s <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
d <- MASS::mvrnorm(1e5, mu = c(0, 0), Sigma = s)
d <- setNames(as.data.frame(d), c("x", "y"))

# fit a simple linear model
m <- lm(y ~ x, data = d)

# generate y predictions over a grid of 10 x values
dpred <- data.frame(
  x = seq(min(d$x), max(d$x), length.out = 10)
)
dpred$yhat <- predict(m, newdata = dpred)

ui <- fluidPage(
  plotlyOutput("scatterplot"),
  checkboxInput(
    "smooth", 
    label = "Overlay fitted line?", 
    value = FALSE
  )
)

server <- function(input, output, session) {
  
  output$scatterplot <- renderPlotly({
    
    p <- plot_ly(d, x   = ~x, y = ~y) %>%
      add_markers(color = I("black"), alpha = 0.05) %>%
      toWebGL() %>%
      layout(showlegend = FALSE)
    
    if (!input$smooth) return(p)
    
    add_lines(p, data = dpred, x = ~x, y = ~yhat, color = I("red"))
  })
  
}

shinyApp(ui, server)

# 17------------------------------------------- #
library(shiny)
library(plotly)

# Generate 100,000 observations from 2 correlated random variables
s <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
d <- MASS::mvrnorm(1e5, mu = c(0, 0), Sigma = s)
d <- setNames(as.data.frame(d), c("x", "y"))

# fit a simple linear model
m <- lm(y ~ x, data = d)

# generate y predictions over a grid of 10 x values
dpred <- data.frame(
  x = seq(min(d$x), max(d$x), length.out = 10)
)
dpred$yhat <- predict(m, newdata = dpred)

ui <- fluidPage(
  plotlyOutput("scatterplot"),
  checkboxInput(
    "smooth", 
    label = "Overlay fitted line?", 
    value = FALSE
  )
)


server <- function(input, output, session) {
  
  output$scatterplot <- renderPlotly({
    plot_ly(d, x = ~x, y = ~y) %>%
      add_markers(color = I("black"), alpha = 0.05) %>%
      toWebGL()
  })
  
  observe({
    if (input$smooth) {
      # this is essentially the plotly.js way of doing
      # `p %>% add_lines(x = ~x, y = ~yhat) %>% toWebGL()`
      # without having to redraw the entire plot
      plotlyProxy("scatterplot", session) %>%
        plotlyProxyInvoke(
          "addTraces", 
          list(
            x = dpred$x,
            y = dpred$yhat,
            type = "scattergl",
            mode = "lines",
            line = list(color = "red")
          )
        )
    } else {
      # JavaScript index starts at 0, so the '1' here really means
      # "delete the second traces (i.e., the fitted line)"
      plotlyProxy("scatterplot", session) %>%
        plotlyProxyInvoke("deleteTraces", 1)
    }
  })
}

shinyApp(ui, server)


# -------Modifying traces------------ #
library(shiny)
library(plotly)

ui <- fluidPage(
  sliderInput("marker", "Marker size", min = 0, max = 20, value = 8),
  sliderInput("path", "Path size", min = 0, max = 30, value = 2),
  plotlyOutput("p")
)

server <- function(input, output, session) {
  
  output$p <- renderPlotly({
    plot_ly(
      economics, x = ~pce, y = ~psavert, z = ~unemploy, 
      color = ~as.numeric(date), mode = "markers+lines"
    )
  })
  
  observeEvent(input$marker, {
    plotlyProxy("p", session) %>%
      plotlyProxyInvoke(
        "restyle", 
        list(marker.size = input$marker)
      )
  })
  
  observeEvent(input$path, {
    plotlyProxy("p", session) %>%
      plotlyProxyInvoke(
        "restyle", list(line.width = input$path)
      )
  })
  
}

shinyApp(ui, server)

# ----------Updating the layout------------ #
plotly_example("shiny", "proxy_mapbox",edit = T)

library(shiny)
library(plotly)

# get all the available mapbox styles
mapStyles <- schema()$layout$layoutAttributes$mapbox$style$values

ui <- fluidPage(
  selectInput("style", "Select a mapbox style", mapStyles),
  plotlyOutput("map")
)

server <- function(input, output, session) {
  
  output$map <- renderPlotly({
    plot_mapbox()
  })
  
  observeEvent(input$style, {
    plotlyProxy("map", session) %>%
      plotlyProxyInvoke(
        "relayout",
        list(mapbox = list(style = input$style))
      )
  })
  
}

shinyApp(ui, server)

# -------------- #
plotly_example("shiny", "proxy_relayout")

# --------Streaming data---------- #
plotly_example("shiny", "stream")


# ----------  #
library(shiny)
library(plotly)

# Generate 100,000 observations from 2 correlated random variables
s <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
d <- MASS::mvrnorm(1e5, mu = c(0, 0), Sigma = s)
d <- setNames(as.data.frame(d), c("x", "y"))

# fit a simple linear model
m <- lm(y ~ x, data = d)

# generate y predictions over a grid of 10 x values
dpred <- data.frame(
  x = seq(min(d$x), max(d$x), length.out = 10)
)
dpred$yhat <- predict(m, newdata = dpred)

ui <- fluidPage(
  plotlyOutput("scatterplot"),
  checkboxInput(
    "smooth", 
    label = "Overlay fitted line?", 
    value = FALSE
  )
)

# --------------- #
demo("sf-dt", package = "plotly")

# --------------- #
library(shiny)
library(dplyr)
library(nycflights13)
library(ggstatsplot)
arr_time <- flights$arr_time
dep_time <- flights$dep_time
arr_bins <- bin_fixed(arr_time, bins = 250)
dep_bins <- bin_fixed(dep_time, bins = 250)
arr_stats <- compute_stat(arr_bins, arr_time) %>%
  filter(!is.na(xmin_))
dep_stats <- compute_stat(dep_bins, dep_time) %>%
  filter(!is.na(xmin_))
ui <- fluidPage(
  plotlyOutput("arr_time", height = 250),
  plotlyOutput("dep_time", height = 250)
)
server <- function(input, output, session) {
  output$arr_time <- renderPlotly({
    plot_ly(arr_stats, source = "arr_time") %>%
      add_bars(x = ~xmin_, y = ~count_)
  })
  output$dep_time <- renderPlotly({
    plot_ly(dep_stats, source = "dep_time") %>%
      add_bars(x = ~xmin_, y = ~count_)
  })
  # arr_time brush updates dep_time view
  observe({
    brush <- event_data("plotly_brushing", source = "arr_time")
    p <- plotlyProxy("dep_time", session)
    # if brush is empty, restore default
    if (is.null(brush)) {
      plotlyProxyInvoke(p, "restyle", "y", list(dep_stats$count_), 0)
    } else {
      in_filter <- between(dep_time, brush$x[1], brush$x[2])
      dep_count <- dep_bins %>%
        compute_stat(dep_time[in_filter]) %>%
        filter(!is.na(xmin_)) %>%
        pull(count_)
      plotlyProxyInvoke(p, "restyle", "y", list(dep_count), 0)
    }
  })
  observe({
    brush <- event_data("plotly_brushing", source = "dep_time")
    p <- plotlyProxy("arr_time", session)
    # if brush is empty, restore default
    if (is.null(brush)) {
      plotlyProxyInvoke(p, "restyle", "y", list(arr_stats$count_), 0)
    } else {
      in_filter <- between(arr_time, brush$x[1], brush$x[2])
      arr_count <- arr_bins %>%
        compute_stat(arr_time[in_filter]) %>%
        filter(!is.na(xmin_)) %>%
        pull(count_)
      plotlyProxyInvoke(p, "restyle", "y", list(arr_count), 0)
    }
  })
}
shinyApp(ui, server)

# -------------- #
plotly_example("shiny", "event_data_persist",edit = T)

library(shiny)
library(plotly)

ui <- fluidPage(
  plotlyOutput("p"),
  tableOutput("table")
)

server <- function(input, output, session) {
  
  # keep track of which cars have been hovered on
  cars <- reactiveVal()
  
  # On hover, the key field of the event data contains the car name
  # Add that name to the set of all "selected" cars
  observeEvent(event_data("plotly_hover"), {
    car <- event_data("plotly_hover")$customdata
    cars_old_new <- c(cars(), car)
    cars(unique(cars_old_new))
  })
  
  # clear the set of cars when a double-click occurs
  observeEvent(event_data("plotly_doubleclick"), {
    cars(NULL)
  })
  
  output$p <- renderPlotly({
    
    # if the car is selected, paint it red
    cols <- ifelse(row.names(mtcars) %in% cars(), "red", "black")
    
    mtcars %>%
      plot_ly(
        x = ~wt, y = ~mpg, 
        customdata = row.names(mtcars), 
        marker = list(color = cols)
      ) %>%
      add_markers()
  })
  
  output$table <- renderTable({
    filter(mtcars, row.names(mtcars) %in% cars())
  })
  
}

shinyApp(ui, server)

# dash bar --------#
library(shinydashboard)
library(shinyWidgets)

renderUI()

ui <- dashboardPage(
  dashboardHeader(title = "Dynamic sidebar"),
  dashboardSidebar(
    sidebarMenu(
      menuItemOutput("menuitem"),
      menuItem("Menu Item 0", tabName = "menu_0"),
      menuItem("Menu Item 1", tabName = "menu_1",
               menuItemOutput("mymenu")),
      menuItem("Menu Item 2", tabName = "menu_2",
               menuItemOutput("mysubmenu"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "menu_0", 
              fluidRow(
                h1("Homepage 000")
              )
      ),
      tabItem(tabName = "menu_1", 
              fluidRow(
                h1("Homepage 1111")
              )
      ),
      tabItem(tabName = "menu_2",
              fluidRow(
                h1("Homepage 222")
              )
      )
    )
  )
)

server <- function(input, output) {
  output$menuitem <- renderMenu({
    menuItem("Menu item",tabName = "menu_0", icon = icon("calendar"))
  })
  output$mysubmenu <- renderMenu({
    menuSubItem('subtext2',tabName = "menu_2")
  })
  output$mymenu <- renderMenu({
    menuItem("my menu",tabName = "menu_1",
             awesomeRadio(
               inputId = "id1", label = "Make a choice:",
               choices = c("graphics", "ggplot2")
             ),
             awesomeCheckboxGroup(
               inputId = "id2", label = "Make a choice:",
               choices = c("graphics", "ggplot2")
             ),
             radioGroupButtons(
               inputId = "id3",
               label = "Label",
               choices = c("A", 
                           "B"),
               justified = TRUE
             ),
             menuSubItem('subtext1',tabName = "menu_1"))
  })
}

shinyApp(ui, server)
# -------------------------- #

ui <- fluidPage(
  uiOutput("moreControls")
)

server <- function(input, output) {
  output$moreControls <- renderUI({
    tagList(
      sliderInput("n", "N", 1, 1000, 500),
      textInput("label", "Label")
    )
  })
}
shinyApp(ui, server)

# ---------------- #
library(shiny)
library(plotly)
ui <- fluidPage(
  selectizeInput(
    inputId = "cities",
    label = "Select a city",
    choices = unique(txhousing$city),
    selected = "Abilene",
    multiple = TRUE
  ),
  plotlyOutput(outputId = "p")
)
server <- function(input, output, ...) {
  output$p <- renderPlotly({
    plot_ly(txhousing, x = ~date, y = ~median) %>% 
      filter(city %in% input$cities) %>%
      group_by(city) %>%
      add_lines()
  })
}
shinyApp(ui, server)

# ----------------------- #
library(shiny)
library(plotly)
library(htmlwidgets)

js <- "function(el, x, inputName){
  var id = el.getAttribute('id');
  var d3 = Plotly.d3;
  $(document).on('shiny:inputchanged', function(event) {
    if (event.name === 'Remove') {
      var out = [];
      d3.select('#' + id + ' g.legend').selectAll('.traces').each(function(){
        var trace = d3.select(this)[0][0].__data__[0].trace;
        out.push([name=trace.name, index=trace.index]);
      });
      Shiny.setInputValue(inputName, out);
    }
  });
}"

ui <- fluidPage(
  textInput("TraceName", "Trace Name"),
  verbatimTextOutput("PrintTraceMapping"),
  actionButton("Add", "Add Trace"),
  actionButton("Remove", "Remove Trace"),
  plotlyOutput("MyPlot")
)

server <- function(input, output, session) {
  
  output$MyPlot <- renderPlotly({
    plot_ly(type = "scatter", mode = "markers") %>%
      layout(showlegend  = TRUE) %>% onRender(js, data = "TraceMapping") 
  })
  
  output$PrintTraceMapping <- renderPrint({unlist(input$TraceMapping)})
  
  observeEvent(input$Add, {
    req(input$TraceName)
    plotlyProxy("MyPlot", session) %>%
      plotlyProxyInvoke("addTraces", list(x = rnorm(10),y = rnorm(10),
                                          type = "scatter",mode = "markers",
                                          name = input$TraceName))
  })
  
  observeEvent(input$Remove, {
    req(input$TraceName, input$TraceMapping)
    traces <- matrix(input$TraceMapping, ncol = 2, byrow = TRUE)
    indices <- as.integer(traces[traces[, 1] == input$TraceName, 2])
    plotlyProxy("MyPlot", session) %>%
      plotlyProxyInvoke("deleteTraces", indices)
  })
  
}

shinyApp(ui, server)

# ---------------------- #
library(shiny)
library(plotly)

ui <- fluidPage(
  selectizeInput(inputId="myTraces", label="Trace names", choices = NULL, multiple = TRUE, options = list('plugins' = list('remove_button'), 'create' = TRUE, 'persist' = TRUE, placeholder = "...add or remove traces")),
  plotlyOutput("MyPlot")
)

server <- function(input, output, session){
  
  myData <- reactiveVal()
  
  observeEvent(input$myTraces, {
    tmpList <- list()
    for(myTrace in input$myTraces){
      tmpList[[myTrace]] <- data.frame(name = myTrace, x = rnorm(10),y = rnorm(10))
    }
    myData(do.call("rbind", tmpList))
    
    return(NULL)
  }, ignoreNULL = FALSE)
  
  output$MyPlot <- renderPlotly({
    if(is.null(myData())){
      plot_ly(type = "scatter", mode = "markers")
    } else {
      plot_ly(myData(), x = ~x, y = ~y, color = ~name, type = "scatter", mode = "markers") %>%
        layout(showlegend  = TRUE)
    }
  })
}

shinyApp(ui, server)


# --------------------------------------------------- #
# 从服务器出发下拉菜单
library("shiny")
library("shinyWidgets")
ui <- fluidPage(
  tags$h2("Toggle Dropdown Button"),
  br(),
  fluidRow(
    column(
      width = 6,
      dropdownButton(
        tags$h3("List of Inputs"),
        selectInput(inputId = 'xcol',
                    label = 'X Variable',choices = names(iris)),
        sliderInput(inputId = 'clusters',
                    label = 'Cluster count',
                    value = 3,min = 1,max = 9),
        actionButton(inputId = "toggle2",
                     label = "Close dropdown"),
        circle = TRUE, status = "danger",
        inputId = "mydropdown",
        icon = icon("gear"), width = "300px"
      )
    ),
    column(
      width = 6,
      actionButton(inputId = "toggle1",
                   label = "Open dropdown")
    )
  )
)
server <- function(input, output, session) {

}
shinyApp(ui = ui, server = server)


# 下拉菜单本身的功能
library("shiny")
library("shinyWidgets")

ui <- fluidPage(
  tags$h2("Dropdown Button"),
  br(),
  dropdown(
    tags$h3("List of Input"),
    fluidRow(
      column(width = 6,
             pickerInput(inputId = 'xcol2',
                         label = 'X Variable',
                         choices = names(iris),
                         options = list(`style` = "btn-info"))
             ),
      column(width = 6,
             pickerInput(inputId = 'ycol2',
                         label = 'Y Variable',
                         choices = names(iris),
                         selected = names(iris)[[2]],
                         options = list(`style` = "btn-warning"))
             )
    ),
    

    

    
    sliderInput(inputId = 'clusters2',
                label = 'Cluster count',
                value = 3,
                min = 1, max = 9),
    
    style = "unite", icon = icon("gear"),
    status = "danger", width = "300px",
    animate = animateOptions(
      enter = animations$fading_entrances$fadeInLeftBig,
      exit = animations$fading_exits$fadeOutRightBig
    )
  ),
  
  plotOutput(outputId = 'plot2')
)

server <- function(input, output, session) {
  
  selectedData2 <- reactive({
    iris[, c(input$xcol2, input$ycol2)]
  })
  
  clusters2 <- reactive({
    kmeans(selectedData2(), input$clusters2)
  })
  
  output$plot2 <- renderPlot({
    palette(c("#E41A1C", "#377EB8", "#4DAF4A",
              "#984EA3", "#FF7F00", "#FFFF33",
              "#A65628", "#F781BF", "#999999"))
    
    par(mar = c(5.1, 4.1, 0, 1))
    plot(selectedData2(),
         col = clusters2()$cluster,
         pch = 20, cex = 3)
    points(clusters2()$centers, pch = 4, cex = 4, lwd = 4)
  })
  
}

shinyApp(ui = ui, server = server)
























