############################
# setwd("D:/Box Sync/Lab/Exomes/Families")
# setwd("D:/Box Sync/Lab/GP2/Analysis")
# dir.create("Rentrez")
# setwd("Rentrez")

# file.create("OMIM_genes_dominant.txt")
# file.create("OMIM_genes_recessive.txt")
############################
# library(rentrez)
# library(tidyverse)
# library(openxlsx)
############################
inheritance <- "dominant"
#inheritance <- "recessive"

disease <- "chorea"
############################

if(file.exists("OMIM_genes_dominant.txt")) {
  dominant_genes <- read_tsv("OMIM_genes_dominant.txt", col_names = FALSE)
  dominant_genes <- pull(dominant_genes)
  dominant_genes <- unique(dominant_genes)}

if(file.exists("OMIM_genes_recessive.txt")) {
  recessive_genes <- read_tsv("OMIM_genes_recessive.txt", col_names = FALSE)
  recessive_genes <- pull(recessive_genes)
  recessive_genes <- unique(recessive_genes)}

set_entrez_key("")
Sys.getenv("ENTREZ_KEY")

find_abstract <- function(pmid) {
  abstract_search <- entrez_search(term = paste(pmid, "[UID]", sep=""), db="pubmed")
  abstract_fetch <- parse_pubmed_xml(entrez_fetch(db="pubmed", id=abstract_search$ids, rettype="xml"))
  return(abstract_fetch$abstract)
}

table <- tibble()
if(inheritance == "dominant") {gene_list <- dominant_genes}else{gene_list <- recessive_genes}
for (i in gene_list) {
  pubmed_part <- tibble(to_delete = "to_delete")
  if(disease == "") {disease_query <- ""}else{disease_query <- paste(" AND", disease)}
  query <- paste0(i, "[TITL] AND ", inheritance, disease_query)
  omim_search <- entrez_search(db="omim", term=query, retmax=10)
  
  if (length(omim_search$ids) != 0) {
    omim_summary <- entrez_summary(db = "omim", id = omim_search$ids)
    medgen_link <- entrez_link(dbfrom="omim", id = omim_search$ids, db="medgen")
    if(!is.null(medgen_link$links$omim_medgen)) {medgen_summary <- entrez_summary(db = "medgen", id = medgen_link$links$omim_medgen)}else{medgen_summary <- tibble(definition = "")}
    gene_part <- tibble(gene = i, disease = if(is.null(omim_summary$title)){disease = "NA"}else{disease = omim_summary$title}, definition = if(length(medgen_summary$definition) == 0 | is.null(medgen_summary$definition)){definition = "NA"}else{definition = medgen_summary$definition[[1]]})
  
    pubmed_link <- entrez_link(dbfrom="omim", id = omim_search$ids, db="pubmed")
    pubmed_cited <- head(pubmed_link$links$omim_pubmed_cited, 5)
  
    for (c in pubmed_cited) {
      pubmed_summary <- entrez_summary(db="pubmed", id = c)
      if (is_empty(find_abstract(c) ==  TRUE)) {abstract <- "NO ABSTRACT"}else{abstract <- find_abstract(c)}
      for (a in 1:length(abstract)) {abstract <- glue::glue_collapse(abstract)}
      pubmed_part <- bind_cols(pubmed_part, tibble(article_title = pubmed_summary$title, article_abstract = abstract, link = paste0("https://pubmed-ncbi-nlm-nih-gov.ezproxy.galter.northwestern.edu/", c, "/")))
    }
  
    pubmed_all <- pubmed_part
    all_parts <- bind_cols(gene_part, pubmed_all)
    
  }else{all_parts <- tibble(gene = i, disease = paste("No association with", inheritance, disease))}
    table <- bind_rows(table, all_parts)
}
if (!is.na(colnames(table)[4])) {table <- select(table, -to_delete)}

style <- createStyle(textDecoration = "bold")
write.xlsx(table, file = paste0("./", "OMIM_table_", inheritance, "_", disease, ".xlsx"), zoom = 115, withFilter = T, headerStyle = style, colWidths = "18", overwrite = T)
