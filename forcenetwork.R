library(networkD3)
library(dplyr)
setwd('/Users/PhilRotella/msan622/Project/')
# Load data
data(MisLinks)
data(MisNodes)

# Plot
forceNetwork(Links = MisLinks, Nodes = MisNodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.8)


drug <- read.csv('data/drug_1May.txt', sep='$')
reac <- read.csv('data/reac_1May.txt', sep = '$')

drugCount <- drug %>% subset(select= c(prod_ai, caseid)) %>% distinct() %>% group_by(prod_ai) %>% summarize(count = n()) %>% arrange(-count)
reacCount <- reac %>% subset(select= c(pt, caseid)) %>% distinct() %>% group_by(pt) %>% summarize(count = n()) %>% arrange(-count)



connection_count <- function(source, target, drugdb, drugindex){
    drug1 <- drugindex[[source + 1, 1]]
    drug2 <- drugindex[[target + 1, 1]]
    return()
}

interactingIDs <- distinct(subset(drug, role_cod == "I", caseid))$caseid
interacting_df <- subset(drug, 
                         caseid %in% interactingIDs & role_cod %in% c('PS', 'I', 'SS'), 
                         c(caseid, prod_ai))

getDrugIDs <-function(target, df) {
    df %>% 
    subset(prod_ai == target, select = caseid) %>%
    distinct() -> id_list
    id_list$caseid
}

idIntersection <- function(source, target, drugIDs) {
    length(intersect(drugIDs$ids[source + 1][[1]], drugIDs$ids[target + 1][[1]]))
}

n_drugs = 100
druglinks <- data.frame('source' = rep(0:(n_drugs - 1), each = n_drugs), 
                        'target' = rep(0:(n_drugs - 1), times = n_drugs),
                        'value' = numeric(n_drugs^2))
druglinks <- subset(druglinks, source != target)
drugIDs <- data.frame("target" = drugCount$prod_ai[1:n_drugs])
drugIDs$ids <- lapply(drugIDs$target, function(x) getDrugIDs(x, interacting_df))
druglinks$value <- lapply(1:nrow(druglinks), function(x) idIntersection(druglinks$source[x], druglinks$target[x], drugIDs))
druglinks <- droplevels(subset(druglinks, value > 0 & source != target))
drugNodes <- as.data.frame(droplevels(drugCount[1:n_drugs, ]))
drugNodes$size <- sapply(drugIDs$ids, length)
drugNodes$group <- numeric(nrow(drugNodes))
#drugNodes <- subset(drugNodes, size > 0)

druglinks$logvalue <- log(as.numeric(druglinks$value))
druglinks <- droplevels(subset(druglinks, logvalue > 0 & source != target))

forceNetwork(Links = druglinks, Nodes = drugNodes, Nodesize = 'size',
             Source = "source", Target = "target",
             Value = "logvalue", NodeID = "prod_ai", Group = 'group' ,
             opacity = 0.6, fontSize = 14, fontFamily = 'Impact,fantasy',
             zoom = T)