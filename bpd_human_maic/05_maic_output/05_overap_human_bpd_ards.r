human_bpd_priority <- ARDSMAICR::inflection_point(human_bpd)

human_bpd_slice <- human_bpd_priority$gene_number

clinical_maic_priority <- ARDSMAICR::inflection_point(clinical_maic)

clinical_maic_slice <- clinical_maic_priority$gene_number

human_bpd_priority_genes <- human_bpd |>
  dplyr::slice(1:human_bpd_slice) |>
  dplyr::pull(gene)

clinical_maic_priority_genes <- clinical_maic |>
  dplyr::slice(1:clinical_maic_slice) |>
  dplyr::pull(gene)

overlap_prioritised <- human_bpd_priority_genes[human_bpd_priority_genes %in% clinical_maic_priority_genes]

num_overlap_prioritised <- length(overlap_prioritised)

percent_overlap_prioritised <- (length(overlap_prioritised)/length(human_bpd_priority_genes))*100

#   ____________________________________________________________________________
#   plots                                                                   ####

euler_plot <- function(eu_obj) {
  plot <- plot(eu_obj,
               quantities = list(type = "counts"),
               shape = "ellipse")
  return(plot)
}

gross_list_genes <- list(BPD = human_bpd_genes, ARDS = clinical_maic_genes)

gross_eu_genes <- eulerr::euler(gross_list_genes)

euler_plot(gross_eu_genes)

prioritised_list_genes <- list(BPD = human_bpd_priority_genes, ARDS = clinical_maic_priority_genes)

prioritised_eu_genes <- eulerr::euler(prioritised_list_genes)

euler_plot(prioritised_eu_genes)

#   ____________________________________________________________________________
#   write out overlap genes                                                 ####


publish_gosttable(func_enrichment, highlight_terms = func_enrichment$result[c(1:2,10,120),],
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size"),
                        filename = NULL)

write.csv(func_enrichment,file='human_func_enrichment.csv',quote=FALSE)

write(func_enrichment, file = "human_func_enrichment.csv",
      ncolumns = if(is.character(func_enrichment)) 1 else 5,
      append = FALSE, sep = ",")


writeLines(overlap_prioritised, "prioritised_overlap.txt")