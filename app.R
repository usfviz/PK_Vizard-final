library(shiny)
library(dplyr)
library(googleVis)
library(ggvis)
library('networkD3')
library(cluster)


drug <- read.csv('data/drug_1May.txt', sep='$')
reac <- read.csv('data/reac_1May.txt', sep = '$')

#lists of drug product names and reaction pts, in order of frequency
drugCount <- drug %>% subset(select= c(prod_ai, caseid)) %>% distinct() %>% group_by(prod_ai) %>% summarize(count = n()) %>% arrange(-count)
reacCount <- reac %>% subset(select= c(pt, caseid)) %>% distinct() %>% group_by(pt) %>% summarize(count = n()) %>% arrange(-count)

top_ten <- drugCount$prod_ai[1:10]
children <- list()
for(j in seq(1,length(top_ten))){
  d <- top_ten[j]
  s <- summary(reac$pt[reac$caseid %in% drug$caseid[drug$prod_ai==d]])[1:10]
  l <- list()
  for(i in seq(1,length(s))) l[[i]] <-list(name=names(s[i]), value=s[[i]])
  children[[j]] <- list(name=d, children=l)
}
radial_list <- list(name="Top Ten Drugs", children=children)

RR_func <- function(a,b,c,d){
  return((a/(a+b))/(c/(c+d)))
}
LLR_func <- function(a,b,c,d){
  return(a*(log(a)-log(a+b))+c*(log(c)-log(c+d))-(a+c)*(log(a+c)-log(a+b+c+d)))
}

PRRdf_func <- function(target, drug, reac) {
    #find caseids that involve target drug (and fit other criteria)
    ids <- drug %>% 
        subset(prod_ai == target, select = caseid) %>%
        distinct()
    ids <- ids$caseid
    
    #calculate pt reaction counts and frequencies among target caseids
    drugAEs <- reac %>%
        subset(caseid %in% ids, select = pt) %>%
        group_by(pt) %>%
        summarize(counts = n())
    drugAEs$freq <- drugAEs$counts / sum(drugAEs$counts)
    
    #calculate comparator frequencies among all other caseids (that fit other criteria)
    otherAEs <- reac %>%
        subset(!(caseid %in% ids), select = pt) %>%
        group_by(pt) %>%
        summarize(counts = n())
    otherAEs$freq <- otherAEs$counts / sum(otherAEs$counts)
    
    #inner join target and comparators and calculate PRR
    PRRdf<- merge(drugAEs, otherAEs, by = "pt", suffixes = c('_target', '_comparator'))
    PRRdf$PRR <- PRRdf$freq_target / PRRdf$freq_comparator
    
    #final results df as pt, drug, count, frequency, and PRR; restrict to at least 3 cases in each
    results <- subset(PRRdf, (counts_target > 9 & counts_comparator > 9), c(pt, counts_target, freq_target, PRR))
    colnames(results) <- c("AE", "Count", "Frequency", "PRR")
    results$Drug <- target
    
    return(results)
}

######################
##     UI           ##
######################
ui <- shinyUI(fluidPage(
    tags$style(type="text/css",".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; }"),
    navbarPage("PK Vizard - FAERS Data Dashboard",
               tabPanel("By Drug",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("drugname", "Select Drug Product(s)", c("All drugs", as.character(drugCount$prod_ai)),
                                            selected = 'All drugs', multiple = T, selectize = F)
                            ),
                            mainPanel(
                                htmlOutput("drug_plot")
                            )
                        )
               ),
               tabPanel("By Event",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("ae", "Select Adverse Event(s)", c("All AEs", as.character(reacCount$pt)),
                                            selected = "All AEs", multiple = T, selectize = F)
                            ),
                            mainPanel(
                                htmlOutput("ae_plot")
                            )
                        )
               ),
               tabPanel("Radial Network of Top 10 Drugs",
                        radialNetworkOutput("radial")),
               tabPanel("By LLR",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("drugname2", "Select Drug Product(s)", c("All drugs", as.character(drugCount$prod_ai)),
                                        selected = 'All drugs', multiple = T, selectize = F)
                          ),
                          mainPanel(
                            dataTableOutput("llr_plot")
                          )
                        )),
               tabPanel("Drug-Event Heatmap",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("hmap_drug", "Select Drug for Sorting", c(as.character(head(drugCount$prod_ai), 20)),
                                            selected = 'ETANERCEPT', multiple = F, selectize = F)
                            ),
                            mainPanel(
                                uiOutput("ggvis_ui"),
                                ggvisOutput("ggvis")
                            )
                        ))
               )
))


###################
##    Server     ##
###################
server <- function(input, output) {
  #all_values <- function(x) {
  #  if(is.null(x)) return(NULL)
  #  row <- final_df[final_df$id == x$id, ]
  #  paste0(row$country, ", ", row$year, collapse = "<br />")
  #}

  #final_df$id <- 1:nrow(final_df)
  #multiplier <- reactive(c(min(final_df$population, na.rm=T),  
  #                         max(final_df$population, na.rm=T) / input$size_scale))
  #df <- reactive({subset(final_df, year==input$year & Region %in% input$region_select)})
  ae_df <- reactive({
      ifelse(input$ae == "All AEs",
             ids <- unique(reac$caseid),
             ids <- unique(subset(reac, pt %in% input$ae)$caseid))
      data <- subset(drug, caseid %in% ids) %>% 
          group_by(prod_ai) %>%
          summarise(count = n()) %>%
          arrange(-count)
      return(head(as.data.frame(data), 10))})
  
  drug_df <- reactive({
      ifelse(input$drugname == "All drugs",
             ids <- unique(drug$caseid),
             ids <- unique(subset(drug, prod_ai %in% input$drugname)$caseid))
      data <- subset(reac, caseid %in% ids) %>% 
          group_by(pt) %>%
          summarise(count = n()) %>%
          arrange(-count)
      return(head(as.data.frame(data), 10))})
  
  llr_df <- reactive({
    ifelse(input$drugname2 == "All drugs",
           ids <- unique(drug$caseid),
           ids <- unique(subset(drug, prod_ai %in% input$drugname2)$caseid))
    data <- subset(reac, caseid %in% ids) %>% 
      group_by(pt) %>%
      summarise(count = n()) %>%
      arrange(-count)
    RR <- c()
    LLR <- c()
    sum_drug <- sum(data$count)
    d <- sum(drugCount$count)
    for(ae in data$pt){
      a <- data$count[data$pt==ae]
      c <- sum_drug-a
      b <- reacCount$count[reacCount$pt==ae]-a
      RR <- append(RR, RR_func(a,b,c,d))
      LLR <- append(LLR, LLR_func(a,b,c,d))
    }
    sorted <- sort(LLR, decreasing=TRUE, index.return=TRUE)
    llr <- sorted$x
    rr <- RR[sorted$ix]
    tops <- data$pt[sorted$ix]
    return(data.frame(tops, llr, rr))
  })
  
  heatmap_df <- PRRdf_func(as.character(drugCount$prod_ai[1]), drug, reac)
  for (target in drugCount$prod_ai[2:20]) {
      PRRdf <- PRRdf_func(as.character(target), drug, reac)
      heatmap_df <-rbind(subset(heatmap_df, AE %in% PRRdf$AE), subset(PRRdf, AE %in% heatmap_df$AE))
  }
  hmap_df <- reactive({
        heatmap_df$AE <- factor(heatmap_df$AE, levels = heatmap_df$AE[order(-subset(heatmap_df, Drug == input$hmap_drug, PRR))], ordered = T)
        return(droplevels(heatmap_df))
    })

  output$drug_plot <- renderGvis({
      gvisBarChart(drug_df())
  })
  output$radial <- renderRadialNetwork({
    radialNetwork(List = radial_list, fontSize = 10, opacity = 0.9)
  })
  
  output$ae_plot <- renderGvis({
      gvisBarChart(ae_df())
  })
  output$llr_plot <- renderDataTable(llr_df())
  
    hmap_df %>%
      ggvis(~Drug, ~AE, fill = ~PRR) %>% 
      layer_rects(width = band(), height = band()) %>%
      scale_nominal("x", padding = 0, points = FALSE) %>%
      scale_ordinal("y", padding = 0, points = FALSE) %>%
      scale_numeric("fill", reverse = T) %>%
      add_axis("x", title = "", properties = axis_props(
          labels = list(angle = -75, align = 'right'))) %>%
      add_axis("y", title = "", properties = axis_props(
          labels = list(fontSize = 8))) %>%
        bind_shiny("ggvis", "ggvis_ui")

}


shinyApp(ui = ui, server = server)
