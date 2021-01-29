set1 <- read.table(".txt", header=T, stringsAsFactors = FALSE)
set2 <- read.table(".txt", header=T, stringsAsFactors = FALSE)

set1_v <- set1$gene_name
set2_v <- set2$gene_name

setEnrichment <- function(set1, set2, universe = 20000){
  a = sum(set1 %in% set2)
  c = length(set1) - a
  
  b = length(set2) - a
  d = universe - length(set2) - c
  
  contingency_table = matrix(c(a, c, b, d), nrow = 2)
  # one-tailed test for enrichment only
  fisher_results = fisher.test(contingency_table, alternative = "greater")
  # returns data.frame containing the lengths of the two sets, the overlap, the enrichment ratio (odds ratio) and P value
  df <- tibble::tibble( set1_length = length(set1), set2_length = length(set2), overlap = a, ratio = fisher_results$estimate, p.value = fisher_results$p.value)
  return(df)
}

setEnrichment( set1 = set1_v, set2 = set2_v)