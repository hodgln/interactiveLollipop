

# server <- 
library(shiny)
library(g3vizModified)
library(glue)
library(httr)
library(jsonlite)
library(tidyverse)
library(purrr)

# change the file path to your own
geneSymbols <- read.csv("/Users/danhodgin/Documents/University/Intercalation/annotation/allData.csv")["Hugo_Symbol"]

listGenes <- unique(geneSymbols$Hugo_Symbol)

index <- which(listGenes == ".")

listGenes <- listGenes[-index]


# change the file path to your own
mutation.csv <- "/Users/danhodgin/Documents/University/Intercalation/annotation/allData.csv"

mutation.dat <- readMAF(mutation.csv,
                        gene.symbol.col = "Hugo_Symbol",
                        variant.class.col = "Variant_Classification",
                        protein.change.col = "amino_acid_change",
                        start.pos.col = "Start_Position",
                        chr.col = "Chromosome",
                        sep = ",",
                        quote = "\"")  # column-separator of csv file



geneSelected <- "TP53" 

# set up chart options



# shinyApp(ui = ui, server = server)

shinyApp(
  ui = fluidPage(
    
     sidebarPanel(selectInput(
      "gene", 
      "Gene:",
                # select the gene names where the column name = hugo...
      c(listGenes)),
    textOutput("result")),
    mainPanel(
      g3LollipopOutput("lolliplot")
    
 )
 ),
  server = function(input, output) {
    geneSelected <- reactive(as.character(input$gene)) 
    
    
    
    output$lolliplot <- renderG3Lollipop({
      
        server <- "https://rest.ensembl.org"
        ext <- "/lookup/symbol/homo_sapiens/"
        
      
        r <- GET(paste(server, ext, geneSelected(), "?expand=1", sep = ""), content_type("application/json"))
        

        parsedData <- head(data.frame(t(sapply(content(r),c))))
        
        
        dataone <- parsedData$Transcript$Transcript
        
        exons <- sapply(dataone, function(i) i["Exon"])
        
        
        newExons <- unique(exons)
        
        
        listExons <- purrr::flatten(newExons)
    
        
        end <- (lapply(listExons, `[[`, "end"))
        start <- (lapply(listExons, `[[`, "start"))
        name <- (lapply(listExons, `[[`, "object_type"))
        
        transcripts <- data.frame(cbind(end, start, name))
        
        #select chromosome using gene. 
        geneData <- subset(mutation.dat, Hugo_Symbol == geneSelected())
        
        #startAndStop <- parsedData$gene
      
        g3Lollipop(mutation.dat,
                   # gene symbol here can be chosen from a dropdown menu
                   gene.symbol = geneSelected(),
                   protein.change.col = "amino_acid_change",
                   start.pos.col = "Start_Position",
                   chr.col = "Chromosome",
                   transcript.values = transcripts,
                   gene.pos.values = c(parsedData$start$start, parsedData$end$end),
                   chr.value = geneData$Chromosome[1],
                   btn.style = "blue", # blue-style chart download buttons
                   plot.options = g3Lollipop.options(
                     # Chart settings
                     chart.width = 900,
                     chart.type = "pie",
                     chart.margin = list(left = 30, right = 20, top = 20, bottom = 30),
                     chart.background = "white",
                     transition.time = 300,
                     # Lollipop track settings
                     lollipop.track.height = 200,
                     lollipop.track.background = "white",
                     lollipop.pop.min.size = 1,
                     lollipop.pop.max.size = 8,
                     lollipop.pop.info.limit = 5.5,
                     lollipop.pop.info.dy = "0.24em",
                     lollipop.pop.info.color = "#d3d3d3",
                     lollipop.line.color = "#a9A9A9",
                     lollipop.line.width = 3,
                     lollipop.circle.color = "#ffdead",
                     lollipop.circle.width = 0.4,
                     lollipop.label.ratio = 2,
                     lollipop.label.min.font.size = 12,
                     lollipop.color.scheme = "dark2",
                     highlight.text.angle = 60,
                     # Domain annotation track settings
                     anno.height = 16,
                     anno.margin = list(top = 0, bottom = 0),
                     anno.background = "white",
                     anno.bar.fill = "#a9a9a9",
                     anno.bar.margin = list(top = 4, bottom = 4),
                     domain.color.scheme = "pie5",
                     domain.margin = list(top = 2, bottom = 2),
                     domain.text.color = "#d3d3d3",
                     domain.text.font = "italic 8px Serif",
                     # Y-axis label
                     y.axis.label = paste("Number of ", geneSelected(), "mutations"),
                     axis.label.color = "#303030",
                     axis.label.alignment = "end",
                     axis.label.font = "italic 12px Serif",
                     axis.label.dy = "-1.5em",
                     y.axis.line.color = "#303030",
                     y.axis.line.width = 0.5,
                     y.axis.line.style = "line",
                     y.max.range.ratio = 1.1,
                     # Chart title settings
                     title.color = "#303030",
                     title.text = paste(geneSelected(), "gene"),
                     title.font = "bold 12px monospace",
                     title.alignment = "start",
                     # Chart legend settings
                     legend = TRUE,
                     legend.margin = list(left=20, right = 0, top = 10, bottom = 5),
                     legend.interactive = TRUE,
                     legend.title = "Variant classification",
                     # Brush selection tool
                     brush = TRUE,
                     brush.selection.background = "#F8F8FF",
                     brush.selection.opacity = 0.3,
                     brush.border.color = "#a9a9a9",
                     brush.border.width = 1,
                     brush.handler.color = "#303030",
                     # tooltip and zoom
                     tooltip = TRUE,
                     zoom = TRUE
                   ),
                   output.filename = "customized_plot")
      })
    
  }

)


