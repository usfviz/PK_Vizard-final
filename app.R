library(shiny)
library(dplyr)
library(googleVis)
library(ggvis)
library(networkD3)


if (!exists('drug')) drug <- read.csv('data/drug_1May.txt', sep='$')
if (!exists('reac')) reac <- read.csv('data/reac_1May.txt', sep = '$')
if (!exists('drugCount'))drugCount <- drug %>% subset(select= c(prod_ai, caseid)) %>% distinct() %>% group_by(prod_ai) %>% summarize(count = n()) %>% arrange(-count)
if (!exists('reacCount'))reacCount <- reac %>% subset(select= c(pt, caseid)) %>% distinct() %>% group_by(pt) %>% summarize(count = n()) %>% arrange(-count) 
if (!exists('counts_df'))counts_df <- read.csv('data/counts_df.csv', row.names = 1)
if (!exists('heatmap_df'))heatmap_df <- read.csv('data/heatmap_df.csv')

if (!(exists('drug') & exists('children') & exists('radial_list'))) {
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
}

RR_func <- function(a,b,c,d){
    return((a/(a+b))/(c/(c+d)))
}
LLR_func <- function(a,b,c,d){
    return(a*(log(a)-log(a+b))+c*(log(c)-log(c+d))-(a+c)*(log(a+c)-log(a+b+c+d)))
}


################ UI ###########################
ui <- shinyUI(fluidPage(
    navbarPage("PK Vizard - FAERS Data Dashboard",
               tabPanel("AE Distributions by Drug",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("drugname", "Select Drug Product(s)", c("ALL DRUGS", as.character(drugCount$prod_ai)),
                                            selected = c('ALL DRUGS', "ETANERCEPT", "ASPIRIN"), multiple = T, selectize = T),
                                sliderInput("p1_n", "Limit to Top X AEs for Each Drug", min = 1, max = 50,
                                            value = 10, step = 1)
                            ),
                            mainPanel(
                                uiOutput("drug_ui"),
                                ggvisOutput("drug_plot")
                            )
                        )
               ),
               tabPanel("Top Medications by AE",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("ae", "Select Adverse Event (AE)", c("All AEs", as.character(reacCount$pt)),
                                            selected = "All AEs", multiple = F, size = 10, selectize = F),
                                sliderInput("p2_n", "Number of Medications to Show", min = 0, max = 1000,
                                            value = 25, step = 1)
                            ),
                            mainPanel(
                                htmlOutput("ae_plot")
                            )
                        )
               ),
               tabPanel("Log-Likelihood Ratios for AE Signals",
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("drugname2", "Select Drug Product(s)", c("ALL DRUGS", levels(drugCount$prod_ai)),
                                            selected = '.ALPHA.-LIPOIC ACID', multiple = T, selectize = F, size = 25)
                            ),
                            mainPanel(
                                dataTableOutput("llr_plot")
                            )
                        )),
               tabPanel("Drug-Event Heatmap",
                        mainPanel(
                            selectInput("hmap_drug", "Select Drug for Sorting", as.character(unique(heatmap_df$Drug)),
                                        selected = 'ETANERCEPT', multiple = F, selectize = F, size = 1),
                            uiOutput("hmap_ui"),
                            ggvisOutput("hmap_plot")
                        )
               )
    )
))


##################### SERVER #############################     
server <- function(input, output) {
    ############ Plot 1 ###################
    drug_df <- reactive({
        df1 <- head(counts_df[counts_df$type == "AE", c('counts', 'freq')], input$p1_n)
        df1 <- data.frame('pt' = row.names(df1), 'counts' = df1$counts, 
                          'freq' = df1$freq, 'alldrugfreq' = df1$freq, 'Drug' = as.character("ALL DRUGS"))
        
        for(drugname in input$drugname) {
            if (!drugname %in% df1$Drug) {
                ids <- unique(drug[drug$prod_ai == drugname, 'caseid'])
                temp_df <- subset(reac, caseid %in% ids) %>% 
                    group_by(pt) %>%
                    summarise(counts = n()) %>%
                    arrange(-counts) %>% head(input$p1_n)
                temp_df$freq <- temp_df$counts / counts_df[drugname,'counts']
                temp_df$alldrugfreq <- sapply(temp_df$pt, function(x) counts_df[x, 'freq'])
                temp_df$Drug <- drugname
                df1 <- rbind(df1, temp_df)
            }
        }
        return (df1[df1$Drug %in% input$drugname, ])
    })
    
    drug_tooltip <- function(x) {
        if(is.null(x)) return(NULL)
        paste0("<strong>",format(x)[1],"</strong><br />",
               format(x)[2]," reports: ", format(as.numeric(format(x)[3])*2, big.mark=",", scientific=FALSE),
               " (", round(as.numeric(format(x)[4]), 2),"%)<br />All reports frequency: ",
               round(as.numeric(format(x)[5]),2), "%<br />")
    }
    
    drug_df %>% arrange(-counts) %>%
        ggvis(~freq*100, ~alldrugfreq*100, size := ~counts*.5, stroke = ~ pt,
              fill = ~ Drug, fillOpacity := .6, fillOpacity.hover := 1, 
              stroke.hover := 'silver') %>%
        add_tooltip(drug_tooltip, on = "hover") %>%
        layer_points() %>% 
        scale_numeric("x", domain = c(0, NA), zero = T, nice = T) %>%
        scale_numeric("y", domain = c(0, NA), zero = T, nice = T) %>%
        scale_ordinal("stroke", range = c('#66666C', '#66666C')) %>%
        add_axis("x", title = "Percent of Target Drug Reports", title_offset = 50,
                 properties = axis_props(title = list(fontSize=15), 
                                         grid = list(stroke = 'white'),
                                         labels = list(align = "center", fontSize = 12))) %>%
        add_axis("y", title = "Percent of All Drug Reports", title_offset = 50,
                 properties = axis_props(title = list(fontSize=15), 
                                         grid = list(stroke = 'white'),
                                         labels = list(align = "right", fontSize = 12))) %>% 
        add_legend('fill', title = 'Drug Product',
                   properties = legend_props(labels = list(fontSize=12),
                                             title = list(fontSize=15))) %>%
        hide_legend('size') %>% 
        hide_legend('stroke') %>% 
        set_options(height = 600, width = 650) %>%
        bind_shiny("drug_plot", "drug_ui")
    
    ############ Plot 2 ###################
    ae_df <- reactive({
        ifelse(input$ae == "All AEs",
               ids <- unique(reac$caseid),
               ids <- unique(subset(reac, pt %in% input$ae)$caseid))
        data <- subset(drug, caseid %in% ids) %>% 
            group_by(prod_ai) %>%
            summarise(count = n()) %>%
            arrange(-count)
        return(head(as.data.frame(data), input$p2_n))})
    
    output$ae_plot <- renderGvis({
        gvisBarChart(ae_df(), 
                     options=list(height = 30 * input$p2_n,
                                  fontSize = 14,
                                  chartArea="{left:200,top:50,width:400,height:\"100%\"}",
                                  legend="none", 
                                  title=paste0("Top ",input$p2_n," Medications Among Reports Involving Selected AE: ", input$ae)))
    })
    
    ############ Plot 3 ###################
    llr_df <- reactive({
        ifelse(input$drugname2 == "ALL DRUGS",
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
    output$llr_plot <- renderDataTable(llr_df())
    
    ########### Plot 4 ##################
    hmap_tooltip <- function(x) {
        if(is.null(x)) return(NULL)
        row <- heatmap_df[heatmap_df$id == x$id, ]
        paste0("AE: ",format(row)[1],"<br />",
               "Drug: ",format(row)[5],"<br />",
               "PRR: ", round(as.numeric(format(row)[4]), 2), "<br />",
               "Count: ", format(as.numeric(format(row)[2]), big.mark=",", scientific=FALSE),
               " (", round(as.numeric(format(row)[3])*100, 2),"%)<br />")
    }
    
    hmap_df <- reactive({
        heatmap_df$AE <- factor(heatmap_df$AE, levels = heatmap_df$AE[order(-subset(heatmap_df, Drug == input$hmap_drug, PRR))], ordered = T)
        return(droplevels(heatmap_df))
    })
    hmap_df %>%
        ggvis(~Drug, ~AE, fill := ~fill, key := ~id,stroke.hover := 'grey') %>%
        add_tooltip(hmap_tooltip, on = 'hover') %>%
        layer_rects(width = band(), height = band()) %>%
        scale_nominal("x", padding = 0, points = FALSE) %>%
        scale_ordinal("y", padding = 0, points = FALSE) %>%
        add_axis("x", title = "", properties = axis_props(
            labels = list(angle = -75, align = 'right'))) %>%
        add_axis("y", title = "", properties = axis_props(
            labels = list(fontSize = 8))) %>%
        set_options(height = 800, width = 900) %>%
        bind_shiny("hmap_plot", "hmap_ui")
}
# Run the application 
shinyApp(ui = ui, server = server)

