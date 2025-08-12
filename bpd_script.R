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
