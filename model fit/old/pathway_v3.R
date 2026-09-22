library(tidyverse)
library(ggpubr)
library(data.table)
library(readr)
library(readxl)
library(formattable)
library(OlinkAnalyze)
library(gtsummary)
library(flextable)
library(pheatmap)
library(plotly)
library(gapminder)
library(manipulateWidget)
library(ggfortify)
library(reactable)
library(mice)
library(naniar)
library(ggrepel)
library(writexl)
library(plotly)
library(OlinkAnalyze)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(fgsea)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(mixOmics)
library(ComplexHeatmap)
library(UpSetR)

setwd("/Users/meizhenyao/Documents/projects/proteomics_ht/data/model/paper1")

# List all CSV files ending with "_noweight_sva.csv"
files <- list.files(pattern = "pth_v3_ad_noweight_robust\\.csv$")

# restricted to proteins with good detection rate
phthalates <- do.call(rbind, lapply(files, read.csv))%>% 
  filter(!Outcome %in% OlinkID_lod_remove75)


phthalates$estimate<- phthalates$Value
phthalates$OlinkID<- phthalates$Outcome


phthalates$Exposure <- factor(phthalates$Exposure,
                              levels = (phthalates %>% distinct(Exposure))$Exposure,
                              labels = toupper(gsub("_SG_3_l","",gsub("pth_", "", (phthalates %>% distinct(Exposure))$Exposure))))

assay_list<- read_excel("/Users/meizhenyao/Documents/projects/proteomics_ht/data/raw/proteins/olink-explore-ht-assay-list.xlsx")
names(assay_list)[1]<- "OlinkID"
names(assay_list)[2]<- "UniProt"
names(assay_list)[3]<- "Protein_name"
names(assay_list)[4]<- "Assay"


protein_lod_remove75<- fread("/Users/meizhenyao/Documents/projects/proteomics_ht/data/intermediate/protein_lod_remove75.txt")
names(protein_lod_remove75)[1]<- "Assay"
names(protein_lod_remove75)[2]<- "OlinkID"
names(protein_lod_remove75)[3]<- "UniProt"
names(protein_lod_remove75)[4]<- "Protein_name"

UniProt_list<- setdiff(assay_list$UniProt, protein_lod_remove75$UniProt)



data<- phthalates %>% 
  left_join(assay_list, by = "OlinkID") 

###### function ######

pos_pathway_info<- function(universe_prot, expo_names, data){
  
  out_list <- vector("list", length(expo_names))
  names(out_list) <- expo_names
  
  entrez_universe_ids <- mapIds(org.Hs.eg.db, keys = universe_prot, keytype = "UNIPROT", column = "ENTREZID")
  entrez_universe_ids<- na.omit(entrez_universe_ids)
  
  for (i in 1:length(expo_names)) {
    
    input_data<- data %>% 
      filter(Exposure == expo_names[i]) %>% 
      filter(p.value < 0.05 & Value > 0)
    
    
    uniprot_sig_name<- input_data$UniProt
    entrez_ids <- mapIds(org.Hs.eg.db, keys = uniprot_sig_name, keytype = "UNIPROT", column = "ENTREZID")
    entrez_ids<- na.omit(entrez_ids)
    
    
    cluster_ora_results<- enrichPathway(gene=entrez_ids, 
                                        universe=entrez_universe_ids,
                                        organism = "human",
                                        pvalueCutoff = 0.99, qvalueCutoff = 0.99,
                                        readable=TRUE, 
                                        pAdjustMethod = "BY")
    
    
    res <- cluster_ora_results@result
    if (is.null(res) || nrow(res) == 0) {
      out_list[[i]] <- tibble::tibble()
    } else {
      out_list[[i]] <- dplyr::mutate(res, exposure = expo_names[i])
    }
  }
  # bind all exposures, keep exposure column
  final <- dplyr::bind_rows(out_list, .id = "ExposureName")
  return(final)
}


neg_pathway_info<- function(universe_prot, expo_names, data){
  
  out_list <- vector("list", length(expo_names))
  names(out_list) <- expo_names
  
  entrez_universe_ids <- mapIds(org.Hs.eg.db, keys = universe_prot, keytype = "UNIPROT", column = "ENTREZID")
  entrez_universe_ids<- na.omit(entrez_universe_ids)
  
  for (i in 1:length(expo_names)) {
    
    input_data<- data %>% 
      filter(Exposure == expo_names[i]) %>% 
      filter(p.value < 0.05 & Value < 0)
    
    
    uniprot_sig_name<- input_data$UniProt
    entrez_ids <- mapIds(org.Hs.eg.db, keys = uniprot_sig_name, keytype = "UNIPROT", column = "ENTREZID")
    entrez_ids<- na.omit(entrez_ids)
    
    
    cluster_ora_results<- enrichPathway(gene=entrez_ids, 
                                        universe=entrez_universe_ids,
                                        organism = "human", 
                                        pvalueCutoff = 0.99, qvalueCutoff = 0.99,
                                        readable=TRUE, 
                                        pAdjustMethod = "BY")
    
    
    res <- cluster_ora_results@result
    if (is.null(res) || nrow(res) == 0) {
      out_list[[i]] <- tibble::tibble()
    } else {
      out_list[[i]] <- dplyr::mutate(res, exposure = expo_names[i])
    }
  }
  # bind all exposures, keep exposure column
  final <- dplyr::bind_rows(out_list, .id = "ExposureName")
  return(final)
}


###############

###### function ######

pos_go_pathway_info<- function(universe_prot, expo_names, data){
  
  out_list <- vector("list", length(expo_names))
  names(out_list) <- expo_names
  
  entrez_universe_ids <- mapIds(org.Hs.eg.db, keys = universe_prot, keytype = "UNIPROT", column = "ENTREZID")
  entrez_universe_ids<- na.omit(entrez_universe_ids)
  
  for (i in 1:length(expo_names)) {
    
    input_data<- data %>% 
      filter(Exposure == expo_names[i]) %>% 
      filter(p.value < 0.05 & Value > 0)
    
    
    uniprot_sig_name<- input_data$UniProt
    entrez_ids <- mapIds(org.Hs.eg.db, keys = uniprot_sig_name, keytype = "UNIPROT", column = "ENTREZID")
    entrez_ids<- na.omit(entrez_ids)
    
    
    
    cluster_ora_results<- enrichGO(gene=entrez_ids, 
                                   universe=entrez_universe_ids,
                                   OrgDb = org.Hs.eg.db,
                                   pvalueCutoff = 0.99, qvalueCutoff = 0.99,
                                   ont           = "BP",
                                   pAdjustMethod = "BY")
    
    cluster_ora_results <- simplify(
      cluster_ora_results,
      cutoff = 0.7,                # similarity threshold (0.5–0.7 is common)
      by = "p.adjust",             # keep terms with smallest p.adjust
      select_fun = min,
      measure = "Wang",            # semantic similarity measure: "Wang", "Resnik", etc.
      semData = NULL
    )
    
    
    res <- cluster_ora_results@result
    if (is.null(res) || nrow(res) == 0) {
      out_list[[i]] <- tibble::tibble()
    } else {
      out_list[[i]] <- dplyr::mutate(res, exposure = expo_names[i])
    }
  }
  # bind all exposures, keep exposure column
  final <- dplyr::bind_rows(out_list, .id = "ExposureName")
  return(final)
}


neg_go_pathway_info<- function(universe_prot, expo_names, data){
  
  out_list <- vector("list", length(expo_names))
  names(out_list) <- expo_names
  
  entrez_universe_ids <- mapIds(org.Hs.eg.db, keys = universe_prot, keytype = "UNIPROT", column = "ENTREZID")
  entrez_universe_ids<- na.omit(entrez_universe_ids)
  
  for (i in 1:length(expo_names)) {
    
    input_data<- data %>% 
      filter(Exposure == expo_names[i]) %>% 
      filter(p.value < 0.05 & Value < 0)
    
    
    uniprot_sig_name<- input_data$UniProt
    entrez_ids <- mapIds(org.Hs.eg.db, keys = uniprot_sig_name, keytype = "UNIPROT", column = "ENTREZID")
    entrez_ids<- na.omit(entrez_ids)
    
    
    cluster_ora_results<- enrichGO(gene=entrez_ids, 
                                   universe=entrez_universe_ids,
                                   OrgDb = org.Hs.eg.db,
                                   pvalueCutoff = 0.99, qvalueCutoff = 0.99,
                                   ont           = "BP",
                                   pAdjustMethod = "BY")
    
    cluster_ora_results <- simplify(
      cluster_ora_results,
      cutoff = 0.7,                # similarity threshold (0.5–0.7 is common)
      by = "p.adjust",             # keep terms with smallest p.adjust
      select_fun = min,
      measure = "Wang",            # semantic similarity measure: "Wang", "Resnik", etc.
      semData = NULL
    )
    
    
    res <- cluster_ora_results@result
    if (is.null(res) || nrow(res) == 0) {
      out_list[[i]] <- tibble::tibble()
    } else {
      out_list[[i]] <- dplyr::mutate(res, exposure = expo_names[i])
    }
  }
  # bind all exposures, keep exposure column
  final <- dplyr::bind_rows(out_list, .id = "ExposureName")
  return(final)
}


###############
UniProt_list<- (data %>% distinct(UniProt))$UniProt

Exposure <- (data %>% distinct(Exposure))$Exposure

############ reactome pos #############

pathway<- pos_pathway_info(UniProt_list, Exposure, data)

write_xlsx(pathway, "/Users/meizhenyao/Documents/projects/proteomics_ht/data/model/paper1/reactome_pos_v3.xlsx")

############ reactome neg #############


pathway<- neg_pathway_info(UniProt_list, Exposure, data)

write_xlsx(pathway, "/Users/meizhenyao/Documents/projects/proteomics_ht/data/model/paper1/reactome_neg_v3.xlsx")


############ go pos #############


pathway<- pos_go_pathway_info(UniProt_list, Exposure, data)

write_xlsx(pathway, "/Users/meizhenyao/Documents/projects/proteomics_ht/data/model/paper1/go_pos_v3.xlsx")

############ go neg #############


pathway<- neg_go_pathway_info(UniProt_list, Exposure, data)

write_xlsx(pathway, "/Users/meizhenyao/Documents/projects/proteomics_ht/data/model/paper1/go_neg_v3.xlsx")
