library(readr)
library(dplyr)
library(ggplot2)
library(data.table)
library(org.Hs.eg.db)
library (EDASeq)
library(biomaRt)
library(RColorBrewer)
library(here)
library(tidyr)

#NOTEPAD:
  #psychiatric ones, ones that I've already presented 
  #both tests for one ancestry for one slide (5 slides per disease)



#Information about ancestry names, population sizes, color palette
Ancestry <- c("AF_afr", "AF_amr", "AF_asj", "AF_eas", "AF_eas_jpn", "AF_eas_kor", "AF_eas_oea", "AF_fin", "AF_nfe", "AF_nfe_bgr", "AF_nfe_est", "AF_nfe_nwe", "AF_nfe_onf", "AF_nfe_seu", "AF_nfe_swe", "AF_oth", "AF_sas")
Ancestry_enhanced <- c("afr", "amr", "asj", "eas", "eas_jpn", "eas_kor", "eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est", "nfe_nwe", "nfe_onf", "nfe_seu", "nfe_swe", "oth", "sas")                                                                                                                                                   
num_exomes <- c(8128, 17296, 5040, 9197, 76, 1909, 7212, 10824, 56885, 1335, 121, 21111, 15499, 5752, 13067, 3070, 15308)
population_sizes <- data.frame(Ancestry , num_exomes)
population_sizes_enhanced <- data.frame(Ancestry_enhanced, num_exomes)
Ancestry1 <- factor(Ancestry, levels = c("afr", "amr", "asj", "eas", "eas_jpn", "eas_kor", "eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est", "nfe_nwe", "nfe_onf", "nfe_seu", "nfe_swe", "oth", "sas"))                                                                                                                                                   
colors <- c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256", "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F")
color_mapping <- data.frame(ancestry = Ancestry_enhanced, color = colors)







#Read in Cosmic and OMIM files
Cosmic.tissue <- readr::read_csv("/blue/kgraim/leslie.smith1/phyloFrame/diseases/census_genes_with_cancer_type_occurence.csv") 
Omim <- readr::read_tsv("/blue/kgraim/leslie.smith1/OMIM/genemap2.txt", skip = 3) 


head(Omim$`Gene/Locus And Other Related Symbols`, 10)
#Clean OMIM file
#Make a separate row for each disease in the phenotype colum and clean the disease names 
omim_data <- Omim %>%
  separate_rows(Phenotypes, sep = ";") %>%
  mutate(Gene = gsub(" ", "", `Gene/Locus And Other Related Symbols`))
omim_data$Phenotypes <- sub("^[^a-zA-Z]*", "", omim_data$Phenotypes)
omim_data$Phenotypes <- sub(",.*", "", omim_data$Phenotypes)
omim_data$Phenotypes <- sub("[^a-zA-Z]+$", "", omim_data$Phenotypes)
omim_data$Phenotypes <- gsub("\\(([^)]*)$", "(\\1)", omim_data$Phenotypes)







#Read in and clean up exome and enhanced exome files
exome_with_germline_mapping <- readr::read_tsv("/blue/kgraim/lucaspereira/ancestry_project/ClinvarMappingResults/exomeAF_with_germline_clinsig.tsv")
exome_with_somatic_mapping <- readr::read_tsv("/blue/kgraim/lucaspereira/ancestry_project/ClinvarMappingResults/clinvar_somatic_with_af.tsv")
exome_with_germline_mapping1 <- exome_with_germline_mapping[-c(8)]
exome_with_somatic_mapping1 <- exome_with_somatic_mapping[-c(10)]

clinvar_germline_with_enhanced_exome <- readr::read_tsv(here("processed-data","clinvar_germline_with_enhanced_exome.tsv"))
clinvar_somatic_with_enhanced_exome <- readr::read_tsv(here("processed-data","clinvar_somatic_with_enhanced_exome.tsv"))
clinvar_germline_with_enhanced_exome <- clinvar_germline_with_enhanced_exome[-c(8)]
clinvar_somatic_with_enhanced_exome <- clinvar_somatic_with_enhanced_exome[-c(10)]

### !!!
clinvar_with_enhanced_exome <- rbind(clinvar_somatic_with_enhanced_exome, clinvar_germline_with_enhanced_exome)
clinvar_with_enhanced_exome$CLNSIG <- toupper(clinvar_with_enhanced_exome$CLNSIG)








#Get information about gene size using biomart
Cosmic.tissue <- data.frame("ENSEMBL" = mapIds(org.Hs.eg.db, keys = Cosmic.tissue$`Gene Symbol`, column = "ENSEMBL", keytype = "SYMBOL"), Cosmic.tissue)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=unique(Cosmic.tissue$ENSEMBL), mart=human)
gene_coords$size=gene_coords$end_position - gene_coords$start_position
colnames(gene_coords)[1] <- "gene"







#get list of associated genes from Cosmic file
get_associated_genes_cosmic <- function(cancer){
  associated_genes <- Cosmic.tissue %>%
    filter(grepl(cancer, Tumour.Types.Germline.) | grepl(cancer, Tumour.Types.Somatic.)) %>%
    distinct(Gene.Symbol) %>%
    pull(Gene.Symbol)
  return(associated_genes)
}

a <- get_associated_genes_omim("Breast Cancer")
b <- get_associated_genes_cosmic("breast")
intersect(a,b)

#get list of associated genes from OMIM file
get_associated_genes_omim <- function(disease){
  disease_data <- omim_data[grep(tolower(disease), tolower(omim_data$Phenotypes)), ]
  genelist <- unlist(strsplit(disease_data$`Gene/Locus And Other Related Symbols`, ", "))
  return(genelist)
}







#Gets all snps with the desired pathogenicity that're associated with the list of genes "genes."
get_associated_snps <- function(genes){
  associated_snps <- snps_with_filtered_clinsig %>%
    filter(gene %in% genes)
  return(associated_snps)
}

get_non_associated_snps <- function(genes){
  non_associated_snps <- snps_with_filtered_clinsig %>%
    filter(!gene %in% genes)
  return(non_associated_snps)
}







#filters pathogenicity in clinvar file
filter_path <- function(df, path){
  return(df[df$CLNSIG %like% toupper(path),])
}








#FOR MUTATION BURDEN ANALYSIS
get_mb <- function(df){
  df.mb <- dplyr::select(df, "AF_nfe_seu":"gene",) %>% group_by(gene) %>% summarise_if(is.numeric, ~length(.[.>0]))
  df.mb.melt <- melt(setDT(df.mb), id.vars = c("gene"), variable.name = "Ancestry")
  print(df.mb.melt)
  return(df.mb.melt)
  
}

scale_data <- function(df){
  df <- merge(df, population_sizes, by="Ancestry")
  df <- merge(df, gene_coords, by="gene")
  df$standardized <- 0
  df$standardized <- df$value/df$num_exomes/df$size*100000
  return(df)
}

graph <- function(data, title1){
  x <- ggplot(data, aes(fill=Ancestry, y=standardized, x=gene)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256",
                               "#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283"))+
    scale_color_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256",
                              "#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283"))+                                                                                                                           
labs(title=title1, x ="Gene", y = "Mutation Burden")
  print(x)
  #ggsave(here("plots", paste(title1, ".png")), width = 14, height = 8, dpi = 300)

}

generate_plot <- function(df, mut_type, cancer, clinsig){
  tit <- paste(mut_type, toupper(cancer), clinsig)
  graph(df, tit)
}

generate_df <- function(mut_type, cancer, clinsig){
  cosmic_somatic <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Somatic.),]
  cosmic_germline <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Germline.),]
  
  common.genes.names <- intersect(cosmic_somatic$Gene.Symbol, cosmic_germline$Gene.Symbol)
  somatic.genes.names <- setdiff(cosmic_somatic$Gene.Symbol , common.genes.names)
  germline.genes.names <- setdiff(cosmic_germline$Gene.Symbol, common.genes.names)
  
  somatic.genes.names
  germline.genes.names
  
  df <- data.frame()
  if(mut_type=="SOMATIC"){
    df <- exome_with_somatic_mapping1[exome_with_somatic_mapping1$gene %in% somatic.genes.names,]
  }
  if(mut_type=="GERMLINE"){
    df <- exome_with_germline_mapping1[exome_with_germline_mapping1$gene %in% germline.genes.names,]
  }
  if(mut_type == "COMMON"){
    df1 <- exome_with_somatic_mapping1[exome_with_somatic_mapping1$gene %in% common.genes.names,]
    df2 <- exome_with_germline_mapping1[exome_with_germline_mapping1$gene %in% common.genes.names,]
    df <- rbind(df1, df2)
  }
  
  head(df)
  df$CLNSIG <- toupper(df$CLNSIG)
  df_with_filtered_clinsig <- filter_path(df, clinsig)
  df_with_filtered_clinsig.mb <- get_mb(df_with_filtered_clinsig)
  df_with_filtered_clinsig.mb <- scale_data(df_with_filtered_clinsig.mb)
  
  return(df_with_filtered_clinsig.mb)
  
}









# FOR ENHANCED MUTATION BURDEN ANALYSIS
get_mb_enhanced <- function(df){
  print("got mb")
  df.mb <- dplyr::select(df, "nfe_seu":"gene",) %>% group_by(gene) %>% summarise_if(is.numeric, ~length(.[.>0]))
  df.mb.melt <- melt(setDT(df.mb), id.vars = c("gene"), variable.name = "Ancestry_enhanced")
  print(df.mb.melt)
  return(df.mb.melt)
}

scale_data_enhanced <- function(df){
  df <- merge(df, population_sizes_enhanced, by="Ancestry_enhanced")
  df <- merge(df, gene_coords, by="gene")
  df$standardized <- 0
  df$standardized <- df$value/df$num_exomes/df$size*100000
  return(df)
}

graph_enhanced <- function(data, title1){
  x <- ggplot(data, aes(fill=Ancestry_enhanced, y=standardized, x=gene)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256",
                               "#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283"))+
    scale_color_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#238B45", "#41AB5D","#74C476","#CE1256",
                                "#CA4136","#67000D", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#825283"))+                                                                                                                           
    labs(title=title1, x ="Gene", y = "Mutation Burden")
  print(x)
  #ggsave(here("plots", paste(title1, ".png")), width = 14, height = 8, dpi = 300)
}

generate_plot_enhanced <- function(df, cancer, clinsig){
  tit <- paste(toupper(cancer), clinsig)
  graph_enhanced(df, tit)
}

generate_df_enhanced <- function(cancer, clinsig){
  somatic <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Somatic),]
  germline <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Germline),]
  
  combined <- rbind(somatic, germline)
  unique_combined <- unique(combined)
  
  df1 <- clinvar_germline_with_enhanced_exome[clinvar_germline_with_enhanced_exome$gene %in% unique_combined$`Gene.Symbol`,]
  df2 <- clinvar_somatic_with_enhanced_exome[clinvar_somatic_with_enhanced_exome$gene %in% unique_combined$`Gene.Symbol`,]
  
  df <- rbind(df1, df2)
  
  
  head(df)
  df$CLNSIG <- toupper(df$CLNSIG)
  df_with_filtered_clinsig <- filter_path(df, clinsig)
  df_with_filtered_clinsig.mb <- get_mb_enhanced(df_with_filtered_clinsig)
  df_with_filtered_clinsig.mb <- scale_data_enhanced(df_with_filtered_clinsig.mb)
  
  return(df_with_filtered_clinsig.mb)
}







#MAKE DENSITY PLOT!!!
#associated snps for one ancestry vs non associated snps for that same ancestry
t_test_1 <- function(cancer, pathogenicity, database){
  if(toupper(database)=="COSMIC"){
    associated_genes <- get_associated_genes_cosmic(cancer)
  } else if(toupper(database)=="OMIM"){
    associated_genes <- get_associated_genes_omim(cancer)
  } else {
    print("invalid database name")
  }
  
  associated_snps <- get_associated_snps(associated_genes, pathogenicity)
  non_associated_snps <- get_non_associated_snps(associated_genes, pathogenicity)
  

  #(nfe, afr, eas, sas, amr)
  nfe_result <- t.test(associated_snps$nfe, non_associated_snps$nfe)
  afr_result <- t.test(associated_snps$afr, non_associated_snps$afr)
  eas_result <- t.test(associated_snps$eas, non_associated_snps$eas)
  sas_result <- t.test(associated_snps$sas, non_associated_snps$sas)
  amr_result <- t.test(associated_snps$amr, non_associated_snps$amr)

  dat <- data.frame("Ancestry"= c("afr", "nfe", "eas", "sas", "amr"), 
                    "P value"= c(nfe_result$p.value, afr_result$p.value, eas_result$p.value, sas_result$p.value, amr_result$p.value), 
                    "df"= c(nfe_result$parameter, afr_result$parameter, eas_result$parameter, sas_result$parameter, amr_result$parameter),
                    "CI lower"= c(nfe_result$conf.int[1], afr_result$conf.int[1], eas_result$conf.int[1], sas_result$conf.int[1], amr_result$conf.int[1]),                    
                    "CI upper"= c(nfe_result$conf.int[2], afr_result$conf.int[2], eas_result$conf.int[2], sas_result$conf.int[2], amr_result$conf.int[2])
                    )
  return(dat)
}

#MAKE DENSITY PLOT (the one i already have)
# associated snps for one ancestry vs associated snps for the other 4 ancestries
t_test_2 <- function(cancer, pathogenicity, database){
  cat("disease: ", cancer, "\n")
  
  if(toupper(database)=="COSMIC"){
    associated_genes <- get_associated_genes_cosmic(cancer)
  } else if(toupper(database)=="OMIM"){
    associated_genes <- get_associated_genes_omim(cancer)
  } else {
    print("invalid database name")
  }
  
  cat("number of associated genes: ", length(associated_genes), "\n")
  
  associated_snps <- get_associated_snps(associated_genes)
  non_associated_snps <- get_non_associated_snps(associated_genes)
  
  non_afr_associated <- dplyr::select(associated_snps, c("ID", "gene", "nfe", "eas", "sas", "amr"),) %>% group_by(ID) 
  non_afr_associated.melt <- melt(setDT(non_afr_associated), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  
  if (length(associated_snps$afr)<2 | nrow(non_afr_associated.melt)<2){
    print("failed check")
    dat <- data.frame("Ancestry" = "NA", "Pvalue" = "NA", "df" = "NA", "CI lower" = "NA", "CI upper" = "NA", "Disease" = cancer)
    return(dat)
    #length(associated_snps$nfe)<2| length(associated_snps$eas)<2| length(associated_snps$sas)<2| length(associated_snps$amr)<2 
    # nrow(non_nfe_associated.melt)<2 | nrow(non_eas_associated.melt)<2 | nrow(non_sas_associated.melt)<2 | nrow(non_amr_associated.melt)<2
    
  }
  
  non_nfe_associated <- dplyr::select(associated_snps, c("ID", "gene", "afr", "eas", "sas", "amr"),) %>% group_by(ID) 
  non_nfe_associated.melt <- melt(setDT(non_nfe_associated), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  
  non_eas_associated <- dplyr::select(associated_snps, c("ID", "gene", "nfe", "afr", "sas", "amr"),) %>% group_by(ID) 
  non_eas_associated.melt <- melt(setDT(non_eas_associated), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  
  non_sas_associated <- dplyr::select(associated_snps, c("ID", "gene", "nfe", "eas", "afr", "amr"),) %>% group_by(ID) 
  non_sas_associated.melt <- melt(setDT(non_sas_associated), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  
  non_amr_associated <- dplyr::select(associated_snps, c("ID", "gene", "nfe", "eas", "sas", "afr"),) %>% group_by(ID) 
  non_amr_associated.melt <- melt(setDT(non_amr_associated), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  

  
  dat <- tryCatch({
    t_test_result_non_afr <- t.test(associated_snps$afr, non_afr_associated.melt$value)
    t_test_result_non_nfe <- t.test(associated_snps$nfe, non_nfe_associated.melt$value)
    t_test_result_non_eas <- t.test(associated_snps$eas, non_eas_associated.melt$value)
    t_test_result_non_sas <- t.test(associated_snps$sas, non_sas_associated.melt$value)
    t_test_result_non_amr <- t.test(associated_snps$amr, non_amr_associated.melt$value)
    dat <- data.frame("Ancestry" = c("afr", "nfe", "eas", "sas", "amr"), 
                      "Pvalue" = c(t_test_result_non_afr$p.value, t_test_result_non_nfe$p.value, t_test_result_non_eas$p.value, t_test_result_non_sas$p.value, t_test_result_non_amr$p.value), 
                      "df" = c(t_test_result_non_afr$parameter, t_test_result_non_nfe$parameter, t_test_result_non_eas$parameter, t_test_result_non_sas$parameter, t_test_result_non_amr$parameter),
                      "CI lower" = c(t_test_result_non_afr$conf.int[1], t_test_result_non_nfe$conf.int[1], t_test_result_non_eas$conf.int[1], t_test_result_non_sas$conf.int[1], t_test_result_non_amr$conf.int[1]), 
                      "CI upper" = c(t_test_result_non_afr$conf.int[2], t_test_result_non_nfe$conf.int[2], t_test_result_non_eas$conf.int[2], t_test_result_non_sas$conf.int[2], t_test_result_non_amr$conf.int[2]),
                      "Disease" = cancer
    )
    return(dat)
  }, error = function(err) {
    dat1 <- data.frame("Ancestry" = "NA", "Pvalue" = "NA", "df" = "NA", "CI lower" = "NA", "CI upper" = "NA", "Disease" = cancer)
    return(dat1)
  })
  
  
  #t_test_result_non_afr <- t.test(associated_snps$afr, non_afr_associated.melt$value)
  #t_test_result_non_nfe <- t.test(associated_snps$nfe, non_nfe_associated.melt$value)
  #t_test_result_non_eas <- t.test(associated_snps$eas, non_eas_associated.melt$value)
  #t_test_result_non_sas <- t.test(associated_snps$sas, non_sas_associated.melt$value)
  #t_test_result_non_amr <- t.test(associated_snps$amr, non_amr_associated.melt$value)

  
  #dat <- data.frame("Ancestry" = c("afr", "nfe", "eas", "sas", "amr"), 
                    #"Pvalue" = c(t_test_result_non_afr$p.value, t_test_result_non_nfe$p.value, t_test_result_non_eas$p.value, t_test_result_non_sas$p.value, t_test_result_non_amr$p.value), 
                    #"df" = c(t_test_result_non_afr$parameter, t_test_result_non_nfe$parameter, t_test_result_non_eas$parameter, t_test_result_non_sas$parameter, t_test_result_non_amr$parameter),
                    ##"CI lower" = c(t_test_result_non_afr$conf.int[1], t_test_result_non_nfe$conf.int[1], t_test_result_non_eas$conf.int[1], t_test_result_non_sas$conf.int[1], t_test_result_non_amr$conf.int[1]), 
                    #"CI upper" = c(t_test_result_non_afr$conf.int[2], t_test_result_non_nfe$conf.int[2], t_test_result_non_eas$conf.int[2], t_test_result_non_sas$conf.int[2], t_test_result_non_amr$conf.int[2]),
                    #"Disease" = cancer
                    #)
  return(dat)
  
}

# all snps associated with a cancer with all ancestries vs all snps not associated with a cancer with all ancestries 
t_test_3 <- function(cancer, pathogenicity){
  associated_genes <- get_associated_genes_cosmic(cancer)
  snps_with_filtered_clinsig <- filter_path(clinvar_with_enhanced_exome, pathogenicity)
  
  # Extract SNPs associated with the selected genes
  associated_snps <- snps_with_filtered_clinsig %>%
    filter(gene %in% associated_genes)
  non_associated_snps <- snps_with_filtered_clinsig %>%
    filter(!gene %in% associated_genes)
  
  associated_select <- dplyr::select(associated_snps, c("ID", "gene", "nfe_seu":"oth"),) %>% group_by(ID)
  associated.melt <- melt(setDT(associated_select), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  
  nonassociated_select <- dplyr::select(non_associated_snps, c("ID", "gene", "nfe_seu":"oth"),) %>% group_by(ID)
  non_associated.melt <- melt(setDT(nonassociated_select), id.vars = c("ID", "gene"), variable.name = "Ancestry")
  
  t_test_result <- t.test(associated.melt$value, non_associated.melt$value)
  
  print(t_test_result)
  
}

#random selection associated snps for one ancestry vs non associated snps for that same ancestry
t_test_4 <- function(tumor_type, pathogenicity, num_permutations){
  associated_genes <- get_associated_genes_cosmic(tumor_type)
  snps_with_filtered_clinsig <- filter_path(clinvar_with_enhanced_exome, pathogenicity)
  
  # Extract SNPs associated with the selected genes (snip set A)
  associated_snps <- snps_with_filtered_clinsig %>%
    filter(gene %in% associated_genes)
  
  non_associated_snps <- snps_with_filtered_clinsig %>%
    filter(!gene %in% associated_genes)
  
  permuted_p_values <- numeric(num_permutations)
  
  set.seed(42)  # For reproducibility
  for (i in 1:num_permutations) {
    # Randomly select the same number of snips as in snip set A
    random_snip_set_B <- snps_with_filtered_clinsig %>%
      sample_n(nrow(associated_snps))
    #make sure those are actually a random set of snips###########################################
    
    # Perform a t-test between snip set A and the random snip set B
    t_test_result <- t.test(associated_snps$afr, random_snip_set_B$afr)
    
    # Store the p-value from this permutation
    permuted_p_values[i] <- t_test_result$p.value
  }
  
  # Calculate the empirical p-value
  observed_t_test_result <- t.test(associated_snps$afr, non_associated_snps$afr)
  empirical_p_value <- mean(permuted_p_values <= observed_t_test_result$p.value)
  print(permuted_p_values)
  # Print the observed t-test result and empirical p-value
  print("Observed T-Test Result:")
  print(observed_t_test_result)
  print("Empirical P-Value:")
  print(empirical_p_value)
  
}








#density plot for snps associated with a disease for this ancestry vs other ancestries
density_plot_2 <- function(cancer, pathogenicity, a, low, high, database){
  if(toupper(database)=="COSMIC"){
    associated_genes <- get_associated_genes_cosmic(cancer)
  } else if(toupper(database)=="OMIM"){
    associated_genes <- get_associated_genes_omim(cancer)
  } else {
    print("invalid database name")
  }
  
  associated_snps <- snps_with_filtered_clinsig %>%
    filter(gene %in% associated_genes)
  print(associated_snps)
  
  data_long_associated <- gather(associated_snps, ancestry, EAF, nfe_seu:oth, factor_key = TRUE)
  
  data_long_associated$ancestry <- factor(data_long_associated$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn",
                                                              "eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))
  
  data_long_associated$group <- ifelse(data_long_associated$ancestry == a, "this ancestry", "other ancestries")
  
  print(data_long_associated)

  e <- ggplot(data = data_long_associated, aes(x = EAF, color = group, fill = group)) + 
    ggtitle(paste("Density Plot")) +
    geom_density(alpha = 0.7) +
    theme_minimal() + 
    scale_fill_manual(values = c("this ancestry" = color_mapping$color[color_mapping$ancestry == a], "other ancestries" = "#FFA92F")) +
    scale_color_manual(values = c("this ancestry" = color_mapping$color[color_mapping$ancestry == a], "other ancestries" = "#FFA92F")) +
    scale_x_continuous(limits = c(low, high)) + scale_y_continuous(limits = c(0, 1000))
  
  print(e)
  
}


#density plot for snps associated with a disease vs snps that are not assocaited with a disease for one ancestry
density_plot_1 <- function(cancer, pathogenicity, a, low, high){
  associated_genes <- get_associated_genes_cosmic(cancer)
  
  snps <- clinvar_with_enhanced_exome %>%
    mutate(Association = ifelse(gene %in% associated_genes, "associated", "non-associated"))
  
  snps <- filter_path(snps, pathogenicity)
  
  data_long <- gather(snps, ancestry, EAF, nfe_seu:oth, factor_key = TRUE)
  data_long$ancestry <- factor(data_long$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn",
                                                              "eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))
  
  data_long_filtered <- data_long[data_long$ancestry == a, ]
  

  e <- ggplot(data = data_long_filtered, aes(x = EAF, color = Association, fill = Association)) + 
    ggtitle(paste("Density of snps from this ancestry that are associated vs not associated with this cancer")) +
    geom_density(alpha = 0.7) + 
    theme_minimal() + 
    scale_fill_manual(values = c("associated" = color_mapping$color[color_mapping$ancestry == a], "non-associated" = "#FFA92F")) +
    scale_color_manual(values = c("associated" = color_mapping$color[color_mapping$ancestry == a], "associated" = "#FFA92F")) +
    scale_x_continuous(limits = c(low, high)) + scale_y_continuous(limits = c(0, 1000))
  
  print(e)

}







#General density plot for all snps 
generate_density_plot <- function(cancer, low, high){
  cosmic_this_cancer <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Germline.),]
  
  print(unique(cosmic_this_cancer$`Gene.Symbol`))
  #snps <- clinvar_germline_with_enhanced_exome[clinvar_germline_with_enhanced_exome$gene %in% cosmic_this_cancer$`Gene.Symbol`,]
  snps <- clinvar_somatic_with_enhanced_exome[clinvar_somatic_with_enhanced_exome$gene %in% cosmic_this_cancer$`Gene.Symbol`,]
  
  print(nrow(snps))
  
  data_long <- gather(snps, ancestry, EAF, nfe_seu:oth, factor_key = TRUE)
  
  print(nrow(data_long))
  
  data_long$ancestry <- factor(data_long$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn",
                                                              "eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))
  e <- ggplot(data = data_long, aes(x=EAF, fill = ancestry, color = ancestry)) + 
    ggtitle(paste("Density Plot for", cancer, "Cancer")) +
    geom_density(alpha = 0.7) + 
    theme_minimal() + 
    scale_fill_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                                    "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F")) + 
    scale_color_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                                   "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F")) + 
    scale_x_continuous(limits = c(low, high)) + scale_y_continuous(limits = c(0, 1000))
  #ggsave(here("plots", paste("Density Plot for", cancer, "Cancer.png")), width = 14, height = 8, dpi = 300)
  return(e)
}


disease_list <- unique(omim_data$Phenotypes)
disease_list <- sort(disease_list)

snps_with_filtered_clinsig <- filter_path(clinvar_with_enhanced_exome, "PATHOGENIC")

ttest2_omim <- data.frame(Pvalue = numeric(), df = numeric(), CI_lower = numeric(), CI = numeric(), disease = character(), stringsAsFactors = FALSE)

ttest2_omim$Pvalue_adj_BH <- p.adjust(ttest2_omim$Pvalue, "BH")
ttest2_omim$Pvalue_adj_bon <- p.adjust(ttest2_omim$Pvalue, "bonferroni")
ttest2_omim$Pvalue_adj_holm <- p.adjust(ttest2_omim$Pvalue, "holm")
ttest2_omim$Pvalue_adj_hochberg <- p.adjust(ttest2_omim$Pvalue, "hochberg")
ttest2_omim$Pvalue_adj_hommel <- p.adjust(ttest2_omim$Pvalue, "hommel")
ttest2_omim$Pvalue_adj_BY <- p.adjust(ttest2_omim$Pvalue, "BY")
ttest2_omim$Pvalue_adj_fdr <- p.adjust(ttest2_omim$Pvalue, "fdr")


sorted_pvalue <- ttest2_omim[order(ttest2_omim$Pvalue_adj_holm), ]

get_associated_snps(get_associated_genes_omim("Macrocephaly/megalencephaly syndrome"))

sorted_pvalue


n <- length(disease_list)
for (i in 1:n) {
  print(i)
  result_df <- t_test_2(disease_list[i], "pathogenic", "omim")
  ttest2_omim <- rbind(ttest2_omim, result_df)
}


for(i in 1:n){
  disease <- disease_list[i]
  if(length(get_associated_genes_omim(disease))==0){
    disease_data <- omim_data[grep(tolower(disease), tolower(omim_data$Phenotypes)), ]
    print(tolower(disease))
    print(disease_data$`Gene/Locus And Other Related Symbols`)
  }
}

print(ttest_2_omim)
write.table(ttest2_omim,"/blue/kgraim/rachel.fenton/ttest", col.names = TRUE, row.names = FALSE, sep = "\t")

density_plot_2("Lung Cancer", "PATHOGENIC", "afr", 0, 0.1, "OMIM")
t_test_2("Lung Cancer", "pathogenic", "omim")



x <- omim_data[grep(tolower("Ectodermal dysplasia 16"), tolower(omim_data$Phenotypes)), ]
x$Phenotypes

get_associated_genes_cosmic("Breast")
Cosmic.tissue

sorted_pvalue <- ttest2_omim[order(ttest2_omim$Pvalue), ]
sorted_pvalue








print(typeof(associated_snps$afr))
print(typeof(non_afr_associated.melt))
print(disease)
print(length(associated_snps))
print(length(non_associated_snps))
print(length(associated_snps$afr))
print(length(associated_snps$nfe))
print(length(associated_snps$eas))
print(length(associated_snps$sas))
print(length(associated_snps$amr))
print(nrow(non_afr_associated.melt))
print(nrow(non_nfe_associated.melt))
print(nrow(non_eas_associated.melt))
print(nrow(non_sas_associated.melt))
print(nrow(non_amr_associated.melt))
