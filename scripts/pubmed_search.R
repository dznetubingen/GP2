#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
INPUT = args[1]
FAMID = args[2]


require(openxlsx)
require(dplyr)
require(rentrez)
require(rlang)

#https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
##########################################
keywords <- c("","parkinson","substianal nigra","dopaminergic","movement disorders")
##########################################

#######################
max_publications <- 5
sort_by <- "relevance"
#sort_by <- "pub+date"
#######################

df <- read_excel(INPUT, sheet = "tier_1.1")
pubmed_genes <- unique(df$gene)

set_entrez_key("1c960f401d0d943d2a42324f8201a0be8308")
Sys.getenv("ENTREZ_KEY")

find_abstract <- function(pmid) {
  abstract_search <- entrez_search(term = paste(pmid, "[UID]", sep=""), db="pubmed")
  abstract_fetch <- parse_pubmed_xml(entrez_fetch(db="pubmed", id=abstract_search$ids, rettype="xml"))
  return(abstract_fetch$abstract)}

dfList <- list()
for (key in keywords){
  pubmed_table <- tibble()
  for (i in pubmed_genes) {
    if(key == "") {keyword_query <- ""}else{keyword_query <- paste(" AND", key)}
    query <- paste0(i, keyword_query)
    pubmed_part <- tibble(to_delete = "to_delete")
    pubmed_search <- entrez_search(db="pubmed", term = query, retmax = max_publications, sort = sort_by)
    Sys.sleep(0.1) 
    if(length(pubmed_search$ids) != 0) {
      for (t in pubmed_search$ids) {
        pubmed_summary <- entrez_summary(db = "pubmed", id = t)
        Sys.sleep(0.1) 
        title <- pubmed_summary$title
        if (is_empty(find_abstract(t) ==  TRUE)) {abstract <- "NO ABSTRACT"}else{abstract <- find_abstract(t)}
        for (a in 1:length(abstract)) {abstract <- glue::glue_collapse(abstract)}
        pubmed_part <- bind_cols(pubmed_part, tibble(title = title, abstract = abstract, link = paste0("https://pubmed.ncbi.nlm.nih.gov/", t, "/")))
        if (pubmed_part[1] == "to_delete") {pubmed_part <- select(pubmed_part, -to_delete)}
        colnames(pubmed_part) <- c(1:ncol(pubmed_part))
        gene_part <- tibble(gene = i)
        all_parts <- bind_cols(gene_part, pubmed_part)
      }
    }else{all_parts <- tibble(gene = i, "1" = paste("No association with", key))}
    pubmed_table <- bind_rows(pubmed_table, all_parts)
  }
    n_articles <- (ncol(pubmed_table) -1) / 3
    columns <- rep(c("article", "abstract", "link"), n_articles)
    colnames(pubmed_table) <- c("gene", columns)
    if(length(columns) != 0){colnames(pubmed_table) <- c("gene", columns)}else{colnames(pubmed_table) <- c("gene", "article")}
    if(key == "") {dfList[['gene']] <-pubmed_table}else{dfList[[paste0("gene_",key)]] <-pubmed_table}
}

style <- createStyle(textDecoration = "bold")
write.xlsx(dfList, paste0("reports/",FAMID,"_pubmed_table.xlsx"),zoom = 115, withFilter = T, headerStyle = style, colWidths = "18", overwrite = T)
