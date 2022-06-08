# Metaborcoding fragments Apps 

library(shiny)
library(ggplot2)
library(shinydashboard)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Amplicon Sequencing"),
  dashboardSidebar(
    numericInput("frag", "Amplicon length", value = 420, step = 1),
    numericInput("fseq", "Forward Sequence length", value = 250, step = 1),
    numericInput("rseq", "Forward Sequence length", value = 250, step = 1),
    textInput("fprim", "Forward primer", value = "Enter primer sequence"),
    textInput("rprim", "Reverse primer", value = "Enter primer sequence")
  ),
  dashboardBody(
    fluidRow(
    box(
      title = "Sequences visualization", #background = "teal", solidHeader = TRUE,
      width = 12, status="primary",
      plotOutput("plot", height = 500)
    )
  )
)
)

server <- function(input, output, session) {
  output$plot <- renderPlot({
    if (input$frag-nchar(input$rprim) > input$fseq & input$frag-nchar(input$fprim) > input$rseq){
      over_forw <- -input$frag/2+input$fseq
      over_rev <- input$frag/2-input$rseq
    }
    else if (input$frag-nchar(input$rprim) > input$fseq & input$frag-nchar(input$fprim) <= input$rseq) {
      over_forw <- -input$frag/2+input$fseq
      over_rev <- -input$frag/2+nchar(input$fprim)
    }
    
    else if (input$frag-nchar(input$rprim) <= input$fseq & input$frag-nchar(input$fprim) > input$rseq){
      over_forw <- input$frag/2-nchar(input$rprim)
      over_rev <- input$frag/2-input$rseq
    }
    else {
      over_forw <- input$frag/2-nchar(input$rprim)
      over_rev <- -input$frag/2+nchar(input$fprim)
    }
    
    over_length <- over_forw-(over_rev)
    
    ggplot()+
      #Amplicon
      geom_segment(aes(x=-input$frag/2+nchar(input$fprim), y=0, xend=input$frag/2-nchar(input$rprim), yend=0, color = "#7F7F7F"),size=4)+
      
      #Primers on amplicon
      geom_segment(aes(x=-input$frag/2, y=0, xend=-input$frag/2+nchar(input$fprim), yend=0),size=4, color = "#F39C12")+
      geom_segment(aes(x=input$frag/2-nchar(input$rprim), y=0, xend=0+input$frag/2, yend=0),size=4, color = "#F39C12")+
      
      #Primers on sequences fragments
      #geom_segment(aes(x=0-input$frag/2-nchar(input$fprim), y=1, xend=0-input$frag/2, yend=1),size=2, color = "#E1B12D")+
      #geom_segment(aes(x=0+input$frag/2, y=-1, xend=0+input$frag/2+nchar(input$rprim), yend=-1),size=2, color = "#E1B12D")+
      
      #Primers on trimmed sequences fragments
      geom_segment(aes(x=-input$frag/2, y=1, xend=-input$frag/2+nchar(input$fprim), yend=1),size=4, color = "#F39C12")+
      geom_segment(aes(x=input$frag/2-nchar(input$rprim), y=-1, xend=+input$frag/2, yend=-1),size=4, color = "#F39C12")+
      
      #Sequenced fragments
      #geom_segment(aes(x=0-input$frag/2, y=1, xend=0-input$frag/2+input$seq-nchar(input$fprim), yend=1),size=2, color = "#357CA5")+
      #geom_segment(aes(x=0+input$frag/2-input$seq+nchar(input$rprim), y=-1, xend=0+input$frag/2, yend=-1),size=2, color = "#357CA5")+
      
      #Trimmed fragments
      geom_segment(aes(x=-input$frag/2+nchar(input$fprim), y=1, xend=-input$frag/2+input$fseq, yend=1),size=4, color = "#0073B7")+
      geom_segment(aes(x=input$frag/2-input$rseq, y=-1, xend=input$frag/2-nchar(input$rprim), yend=-1),size=4, color = "#0073B7")+
      
      #Overlap without trimming 
      #geom_segment(aes(x=0-input$frag/2+input$seq-nchar(input$fprim), y=1, xend=0-input$frag/2+input$seq-nchar(input$fprim), yend=-10),linetype = "dashed", color = "#A1141E")+
      #geom_segment(aes(x=0+input$frag/2-input$seq+nchar(input$rprim), y=1, xend=0+input$frag/2-input$seq+nchar(input$rprim), yend=-10),linetype = "dashed", color = "#A1141E")+
      #geom_segment(aes(x=0+input$frag/2-input$seq+nchar(input$rprim), y=-10, xend=0-input$frag/2+input$seq-nchar(input$fprim), yend=-10),size=2, color = "#A1141E")+
      #geom_text(aes(x=0, y=-11, label = as.character(paste(input$seq-nchar(input$fprim)-input$frag/2+input$seq-nchar(input$rprim)-input$frag/2," bps"))), color = "#A1141E")+
      
      #Overlap with trimming 
      geom_segment(aes(x=over_forw, y=-1, xend=over_forw, yend=10),linetype = "dashed", color = "#A1141E")+
      geom_segment(aes(x=over_rev, y=-1, xend=over_rev, yend=10),linetype = "dashed", color = "#A1141E")+
      geom_segment(aes(x=over_rev, y=10, xend=over_forw, yend=10),size=4, color = "#A1141E")+
      geom_text(aes(x=(over_forw+over_rev)/2, y=11, label = as.character(paste(over_length," bps"))), color = "#A1141E")+
      
      
      ylim(-12,12)+
      scale_color_manual(name='',
                         breaks=c('Amplicon', 'Sequence', 'Primer','Overlap'),
                         values=c('Amplicon'='#7F7F7F', 'Sequence'='#0073B7', 'Primer'='#F39C12','Overlap'='#A1141E'))+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="bottom",
            legend.title=element_text(size=20),
            legend.text=element_text(size=14),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
  }, res = 96)
}

shinyApp(ui, server)
