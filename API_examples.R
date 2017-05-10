# library(ggplot2)
# library(ggvis)
# library(dplyr)
# library(jsonlite)
# library(networkD3)
# library(urltools)
# 
# #key fields:
# # DRUG NAMES/CLASSES
# # patient.drug.openfda.substance_name.exact
# # patient.drug.openfda.brand_name.exact
# # patient.drug.openfda.pharm_class_epc.exact
# # REACTIONS
# # patient.reaction.reactionmeddrapt
# # REPORT CHARACTERISTICS
# # receivedate
# # DRUG CHARACTERISTICS
# # patient.drug.drugindication.exact
# # patient.drug.drugadministrationroute.exact
# # patient.drug.drugcharacterization (1 = suspect, 2 = concomitant, 3 = interacting)
# # 
# # DEMOGRAPHICS
# # "patientonsetage: [5+TO+17]
# # "patientsex": "2", (female = 2, male = 1)
# 


api_key <- 'kFbeLzfFDAAgerR1hbhz65IdE81QbsxlLAlGRvIl'
base_url <- 'https://api.fda.gov/drug/event.json?api_key=%s&search=%s&count=%s'

# 
# # get vector of top 1000 drugs since 2016
search_query <- 'receiptdate:([2016-01-01+TO+2017-04-30])&limit=1000'
url <- sprintf(base_url, api_key, search_query, 'patient.drug.openfda.substance_name.exact')
top_1000_drug_response <- fromJSON(url)
top1000drug_df <- top_1000_drug_response$results
top1000drugs <- top1000drug_df$term

ggplot(head(top1000drug_df), aes(x = term, y = count)) + geom_col()
# get vector of top 1000 AEs since 2016
search_query <- 'receiptdate:([2016-01-01+TO+2017-04-30])&limit=1000'
url <- sprintf(base_url, api_key, search_query, 'patient.reaction.reactionmeddrapt.exact')
top_1000_AEs_response <- fromJSON(url)
top1000AEs_df <- top_1000_AEs_response$results
 top1000AEs <- top1000AEs_df$term

test_DRUG_AEs <- fromJSON(url_encode(url))
url <- sprintf(base_url, api_key, paste0(search_query, "+AND+patient.reaction.reactionmeddrapt.exact:",top1000AEs[1]), 'patient.drug.openfda.substance_name.exact')



