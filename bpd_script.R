library('clusterProfiler')
library(ARDSMAICR)
library(data.table) # 1.15.4
library(dplyr) # 1.1.4
library(eulerr) # 7.0.0
library(readr) # 2.1.5

====
#Read in and prioritise
====
human_maic_data <- data.table::fread('/Users/sclohise/Dropbox/03_BPD/BPD_maic_processed/05_maic_output/human_bpd_250210_v4/MAIC_input_2025-02-10-2025-02-10--11-08/maic_raw.txt')
human_priority <- ARDSMAICR::inflection_point(human_maic_data)
human_priority_plot <- ARDSMAICR::inflection_point_plot(human_maic_data, first_break = 1500, increment = 2000)
#image
tiff(file="human_priority.tiff",
     width=12, height=8, units="in", res=100)
human_priority_plot
dev.off()

====
rodent_maic_data <- data.table::fread('/Users/sclohise/Dropbox/03_BPD/BPD_maic_processed/05_maic_output/rodent_bpd_250409_v5/maic_raw.txt')
rodent_priority <- ARDSMAICR::inflection_point(rodent_maic_data)
rodent_priority_plot <- ARDSMAICR::inflection_point_plot(rodent_maic_data, first_break = 1000, increment = 2000)
#image
tiff(file="rodent_priority.tiff",
     width=12, height=8, units="in", res=100)
rodent_priority_plot
dev.off()


====
#Overlap
====

  
#Extract genes  
gene_list_human <- human_maic_data$gene  
gene_list_rodent <- rodent_maic_data$gene  

#Eulerplot function (JM)

euler_plot <- function(eu_obj) {
  plot <- plot(eu_obj,
               quantities = list(type = "counts"),
               shape = "ellipse")
  return(plot)
}

#Overlap 

overlap <- gene_list_rodent[gene_list_rodent %in% gene_list_human]
gross_list_genes <- list(Rodent = gene_list_rodent, Human = gene_list_human)
gross_eu_genes <- eulerr::euler(gross_list_genes)
euler_plot(gross_eu_genes)
#image
tiff(file="gross_euler_overlap.tiff",
     width=6, height=4, units="in", res=100)
euler_plot(gross_eu_genes)
dev.off()

# Get prioritised gene list
  
human_bpd_maic_slice <- human_priority$gene_number
human_priority_genes <- human_maic_data |>
     dplyr::slice(1:human_bpd_maic_slice) |>
     dplyr::pull(gene)

rodent_bpd_maic_slice <- rodent_priority$gene_number
rodent_priority_genes <- rodent_maic_data |>
  dplyr::slice(1:rodent_bpd_maic_slice) |>
  dplyr::pull(gene)

#Overlap prioritised gene lists

prioritised_list_genes <- list(Rodent = rodent_priority_genes, Human = human_priority_genes)
prioritised_eu_genes <- eulerr::euler(prioritised_list_genes)
euler_plot(prioritised_eu_genes)
#image
tiff(file="prioritised_euler_overlap_ards_bod.tiff",
     width=6, height=4, units="in", res=100)
euler_plot(prioritised_eu_genes)
dev.off()

====
#Semantic similarity
====
  
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="BP")


======
  
library("org.Hs.eg.db")
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(DOSE)
library(ReactomePA)

  
h_ego <- enrichGO(gene= human_priority_genes,
                    OrgDb= org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont= "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE
                    
)
#simplify based on redundant GO terms ()
s_h_ego<-clusterProfiler::simplify(h_ego)


human_go_ora_dotplot <- dotplot(s_h_ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO ORA Enrichment of Human BPD genes")
#image
tiff(file="human_GO_ORA.tiff",
     width=12, height=8, units="in", res=100)
human_go_ora_dotplot 
dev.off()

#Functional Enrichment
##Trying to rearrange data to be suitable for cluster profilr
#genes with gene score (maic score)

human_maic_scores <- data.frame('genes' = human_maic_data$gene, 
                     'maic_score' = human_maic_data$maic_score 
                     )
human_priority_entrez = bitr(human_priority_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
human_maic_scores_fe = human_maic_scores[,2]
names(human_maic_scores_fe) = as.character(human_priority_entrez$ENTREZID)
#order appropriately
human_maic_scores_fe = sort(human_maic_scores_fe, decreasing = TRUE)
#convert to numerical

#GSEA

h_ego2 <- enrichGO(gene         = gene_list_human,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

s_h_ego2<-clusterProfiler::simplify(h_ego2)

human_go_gsea_dotplot <- dotplot(s_h_ego2,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO GSEA Enrichment of Human BPD genes")
#image
tiff(file="human_GO_GSEA.tiff",
     width=12, height=8, units="in", res=100)
human_go_gsea_dotplot 
dev.off()






====

  #KEGG
  
====
  kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)


mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.5)


====
  
  #Reactome

====
  
#ORA
  
x <- enrichPathway(gene=entrez_priority_h, pvalueCutoff = 0.05, readable = TRUE)
human_reactome_ora_dotplot <- dotplot(x,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("Reactome ORA Enrichment of Human BPD genes")
#image
tiff(file="human_Reactome_ORA.tiff",
     width=12, height=8, units="in", res=100)
human_reactome_ora_dotplot
dev.off()

#GSEA

y <- gsePathway(gene_list_human,
                keyType = "SYMBOL",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE)


=====
  #ARDS
  
====

ards_maic <- ARDSMAICR::data_genes
ards_priority <- ARDSMAICR::inflection_point(ards_maic)
ARDSMAICR::inflection_point_plot(ards_maic, first_break = 1500, increment = 2000)

ards_maic_slice <- ards_priority$gene_number
ards_priority_genes <- ards_maic |>
  dplyr::slice(1:ards_maic_slice) |>
  dplyr::pull(gene)


human_prioritised_list_genes <- list(ARDS = ards_priority_genes, BPD = human_priority_genes)
> prioritised_eu_genes <- eulerr::euler(human_prioritised_list_genes)
> euler_plot(prioritised_eu_genes)
> tiff(file="prioritised_euler_overlap_ards_bpd.tiff",
       +      width=6, height=4, units="in", res=100)
> euler_plot(prioritised_eu_genes)
> dev.off()

a_ego <- enrichGO(gene= ards_priority_genes,
                  OrgDb= org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont= "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE
                  
)
#simplify based on redundant GO terms ()
s_a_ego<-clusterProfiler::simplify(a_ego)


human_go_ards_dotplot <- dotplot(s_a_ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO ORA Enrichment of Human ARDS genes")
#image
tiff(file="ards_GO_ORA.tiff",
     width=12, height=8, units="in", res=100)
human_go_ards_dotplot 
dev.off()


mlist<-list(s_h_ego,s_a_ego)
names(mlist)<-c("BPD","ARDS")
mresult<-merge_result(mlist)
dotplot(mresult,showCategory=10) 

ards_bpd_compare <- dotplot(mresult,font.size = 12,showCategory=20, title ="GO ORA Enrichment of BPD and ARDS genes")+ scale_y_discrete(labels=function(x) str_wrap(x, width=80))

#image
tiff(file="ards_bpd_GO_ORA_20.tiff",
     width=12, height=8, units="in", res=100)
ards_bpd_compare 
dev.off()

overlap <- ards_priority_genes[ards_priority_genes %in% human_priority_genes]

o_ego <- enrichGO(gene= overlap,
                  OrgDb= org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont= "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE
                  
)
s_a_ego<-clusterProfiler::simplify(a_ego)


human_go_priority_overlap_ards_dotplot <- dotplot(o_ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO ORA Enrichment of Human ARDS/BPD genes")
#image
tiff(file="ards_BPD_GO_ORA.tiff",
     width=12, height=8, units="in", res=100)
human_go_priority_overlap_ards_dotplot
dev.off()




#==========
  #Hypergeometric overlap
#==========
#ards bpd overlap
n_A = 945;n_B = 1306; n_C = 20000; n_A_B = 112
1-phyper(n_A_B, n_B, n_C-n_B, n_A)

#human v rodent bpd overlap
n_A = 945;n_B = 1783; n_C = 19800; n_A_B = 163
1-phyper(n_A_B, n_B, n_C-n_B, n_A)

#gross human v rodent overlap
n_A = (4889 + 2747) ;n_B = (5270+2747); n_C = 20000; n_A_B = 2747
1-phyper(n_A_B, n_B, n_C-n_B, n_A)


========
  #rodent
========

  
r_ego <- enrichGO(gene= rodent_priority_genes,
                    OrgDb= org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont= "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE
                    
  )
#simplify based on redundant GO terms ()
r_h_ego<-clusterProfiler::simplify(r_ego)


rodent_go_ora_dotplot <- dotplot(r_h_ego,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO ORA Enrichment of Rodent BPD genes")
#image
tiff(file="rodent_GO_ORA.tiff",
     width=12, height=8, units="in", res=100)
rodent_go_ora_dotplot 
dev.off()


mlist<-list(s_h_ego,r_h_ego)
names(mlist)<-c("Human","Rodent")
mresult<-merge_result(mlist)
dotplot(mresult,showCategory=10) 

rodent_human_bpd_compare <- dotplot(mresult,font.size = 12,showCategory=20, label_format=70, title ="GO ORA Enrichment of human and rodent BPD genes")+ scale_y_discrete(labels=function(x) str_wrap(x, width=80))

#image
tiff(file="rodent_human_BPD_GO_ORA.tiff",
     width=12, height=8, units="in", res=100)
rodent_human_bpd_compare
dev.off()









================================
  #FGSEA
===============================  
https://biostatsquid.com/fgsea-tutorial-gsea/  
gene_sets_df <- msigdbr(species = 'Homo sapiens', category = 'H')
================================
#G profiler
===============================  
library(gprofiler2)

gostres <- gost(
    gene_list_human,
    organism = "hsapiens",
    ordered_query = TRUE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = TRUE,
    measure_underrepresentation = TRUE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = c("fdr"),
    domain_scope = c("annotated"),
    custom_bg = NULL,
    numeric_ns = "",
    sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP"),
    as_short_link = FALSE,
    highlight = FALSE
  )

gostplot(
gostres,
capped = FALSE,
interactive = FALSE,
pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG ="#dd4477", REAC = "#3366cc", WP = "#0099c6")
)



===================
  GSEA() Human

====================
library("org.Hs.eg.db")
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")


d = read.csv("human_genelist.csv")
genes = d$gene
entrez_human <- mapIds(org.Hs.eg.db, keys = genes, keytype="SYMBOL", column = "ENTREZID")
geneList_n = d[,2]
names(geneList_n) = entrez_human
human_geneList = sort(geneList_n, decreasing = TRUE)

#GO GSEA

ego3 <- gseGO(geneList     = human_geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

human_go_gsea_dotplot <- dotplot(ego3,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO:BP GSEA of Human BPD genes")
#image
tiff(file="human_GO_BP_gsea.tiff",
     width=12, height=8, units="in", res=100)
human_go_gsea_dotplot
dev.off()
human_go_gsea_dotplot


#Reactome GSEA
human_reacty <- gsePathway(human_geneList, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE)

human_react_gsea_dotplot <- dotplot(human_reacty,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("Reactome GSEA of Human BPD genes")
#image
tiff(file="human_Reactome_gsea.tiff",
     width=12, height=8, units="in", res=100)
human_react_gsea_dotplot
dev.off()
human_react_gsea_dotplot



#image individual pathways, where 1 is the first pathway in the list

react1 <- gseaplot(human_reacty, geneSetID = 1, by = "runningScore", title = human_reacty$Description[1])
tiff(file="human_Reactome_gsea_1.tiff",
     width=12, height=8, units="in", res=100)
react1
dev.off()
react1


react1_7 <- gseaplot2(human_reacty, geneSetID = 1:7)
tiff(file="human_Reactome_gsea_1_7.tiff",
     width=12, height=8, units="in", res=100)
react1_7
dev.off()
react1_7

gseaplot2(human_reacty, geneSetID = 1, title = human_reacty$Description[1])

gseaplot2(human_reacty, geneSetID = 1:7)
#KEGG GSEA
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

human_kegg_gsea_dotplot <- dotplot(kk2,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("KEGG GSEA of Human BPD genes")
#image
tiff(file="human_kegg_gsea.tiff",
     width=12, height=8, units="in", res=100)
human_kegg_gsea_dotplot
dev.off()
human_kegg_gsea_dotplot



===============
GSEA() rodent
===============
library("org.Hs.eg.db")
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")


d = read.csv("rodent_genelist.csv")
genes = d$gene
entrez_rodent <- mapIds(org.Hs.eg.db, keys = genes, keytype="SYMBOL", column = "ENTREZID")
geneList_n = d[,2]
names(geneList_n) = entrez_rodent
rodent_geneList = sort(geneList_n, decreasing = TRUE)

#GO GSEA

ego3 <- gseGO(geneList     = rodent_geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

rodent_go_gsea_dotplot <- dotplot(ego3,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO:BP GSEA of Rodent BPD genes")
#image
tiff(file="rodent_GO_BP_gsea.tiff",
     width=12, height=8, units="in", res=100)
rodent_go_gsea_dotplot
dev.off()
rodent_go_gsea_dotplot


#Reactome GSEA
rodent_reacty <- gsePathway(rodent_geneList, 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH", 
                            verbose = FALSE)

rodent_react_gsea_dotplot <- dotplot(rodent_reacty,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("Reactome GSEA of Rodent BPD genes")
#image
tiff(file="rodent_Reactome_gsea.tiff",
     width=12, height=8, units="in", res=100)
rodent_react_gsea_dotplot
dev.off()
rodent_react_gsea_dotplot

rreact1_10 <- gseaplot2(rodent_reacty, geneSetID = 1:10)
tiff(file="rodent_Reactome_gsea_1_10.tiff",
     width=12, height=8, units="in", res=100)
rreact1_10
dev.off()
rreact1_10




#KEGG GSEA
kk2 <- gseKEGG(geneList     = rodent_geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

rodent_kegg_gsea_dotplot <- dotplot(kk2,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("KEGG GSEA of rodent BPD genes")
#image
tiff(file="rodent_kegg_gsea.tiff",
     width=12, height=8, units="in", res=100)
rodent_kegg_gsea_dotplot
dev.off()
rodent_kegg_gsea_dotplot



==============
  #PPI overlap
=============
library(ggvenn)

hdegree <- (read.csv("human_bpd_ppi/human_bpd_100_degree.csv", skip = 1, header=T))$Name 
hdmnc <- (read.csv("human_bpd_ppi/human_bpd_100_dmnc.csv", skip = 1, header=T))$Name 
hepc <- (read.csv("human_bpd_ppi/human_bpd_100_epc.csv", skip = 1, header=T))$Name 
hmcc <- (read.csv("human_bpd_ppi/human_bpd_100_mcc.csv", skip = 1, header=T))$Name
hmnc <- (read.csv("human_bpd_ppi/human_bpd_100_mnc.csv", skip = 1, header=T))$Name

venndata <- list(hdegree, hepc, hmcc, hmnc, hdmnc)
hppi_upset <- ggVennDiagram(venndata, force_upset = TRUE, order.set.by = "name", category.names = c("Degree","EPC","MCC", "MNC", "DMNC"))

tiff(file="hppi_upset.tiff",
     width=12, height=8, units="in", res=100)
hppi_upset 
dev.off()
hppi_upset 

write.table(
  overlap_all,
  file = "overlapping_genes_human_ppi_hub.txt",
  quote = FALSE,     # don't add quotes around values
  row.names = FALSE, # don't write row numbers
  col.names = FALSE  # don't write a header
)


rdegree <- (read.csv("rodent_bpd_ppi/rodent_bpd_100_degree.csv", skip = 1, header=T))$Name 
rdmnc <- (read.csv("rodent_bpd_ppi/rodent_bpd_100_dmnc.csv", skip = 1, header=T))$Name 
repc <- (read.csv("rodent_bpd_ppi/rodent_bpd_100_epc.csv", skip = 1, header=T))$Name 
rmcc <- (read.csv("rodent_bpd_ppi/rodent_bpd_100_mcc.csv", skip = 1, header=T))$Name
rmnc <- (read.csv("rodent_bpd_ppi/rodent_bpd_100_mnc.csv", skip = 1, header=T))$Name

venndata <- list(rdegree, repc, rmcc, rmnc, rdmnc)
hppi_upset_r <- ggVennDiagram(venndata, force_upset = TRUE, order.set.by = "name", category.names = c("Degree","EPC","MCC", "MNC", "DMNC"))
overlap_all <- Reduce(intersect, venndata)
tiff(file="rppi_upset.tiff",
     width=12, height=8, units="in", res=100)
hppi_upset_r
dev.off()
hppi_upset_r

write.table(
  overlap_all,
  file = "overlapping_genes_rodent_ppi_hub.txt",
  quote = FALSE,     # don't add quotes around values
  row.names = FALSE, # don't write row numbers
  col.names = FALSE  # don't write a header
)


==========
  #STRING
=========  
human_priority_p_mapped <- string_db$map( human_priority_genes, "gene", removeUnmappedRows = TRUE )


