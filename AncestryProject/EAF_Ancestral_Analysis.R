library(readr)
library(dplyr)
library(ggplot2)
library(data.table)
library(org.Hs.eg.db)
library (EDASeq)
library(biomaRt)
library(RColorBrewer)
library(tidyr)
library(here)


Cosmic.tissue <- readr::read_csv("/blue/kgraim/leslie.smith1/phyloFrame/diseases/census_genes_with_cancer_type_occurence.csv") 
clinvar_germline_with_enhanced_exome <- readr::read_tsv(here("processed-data","clinvar_germline_with_enhanced_exome.tsv"))
clinvar_somatic_with_enhanced_exome <- readr::read_tsv(here("processed-data","clinvar_somatic_with_enhanced_exome.tsv"))




Ancestry <- c("nfe_seu", "nfe_bgr", "afr", "sas", "nfe_onf", "amr", "eas", "nfe_swe", "nfe_nwe", "eas_jpn", "eas_kor", "eas_oea", "nfe_est", "nfe", "fin", "asj", "oth")
num_exomes <- c(5752, 1335, 8128, 15308, 15499, 17296, 9197, 13067, 21111, 76, 1909, 7212, 121, 56885, 10824, 5040, 3070)
population_sizes <- data.frame(Ancestry , num_exomes)

#Create gene size lookup table
Cosmic.tissue <- data.frame("ENSEMBL" = mapIds(org.Hs.eg.db, keys = Cosmic.tissue$`Gene Symbol`, column = "ENSEMBL", keytype = "SYMBOL"), Cosmic.tissue)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=unique(Cosmic.tissue$ENSEMBL), mart=human)
gene_coords$size=gene_coords$end_position - gene_coords$start_position
colnames(gene_coords)[1] <- "gene"


filter_path <- function(df, path){
  return(df[df$CLNSIG %like% path,])
}

get_mb <- function(df){
  print("got mb")
  df.mb <- dplyr::select(df, "nfe_seu":"gene",) %>% group_by(gene) %>% summarise_if(is.numeric, ~length(.[.>0]))
  df.mb.melt <- melt(setDT(df.mb), id.vars = c("gene"), variable.name = "Ancestry")
  print(df.mb.melt)
  return(df.mb.melt)

}

head(clinvar_germline_with_enhanced_exome)

scale_data <- function(df){
  print("data scaled")
  df <- merge(df, population_sizes, by="Ancestry")
  df <- merge(df, gene_coords, by="gene")
  df$standardized <- 0
  df$standardized <- df$value/df$num_exomes/df$size*100000
  return(df)
}

graph <- function(data, title){
  ggplot(data, aes(fill=Ancestry, y=standardized, x=gene)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(values=c("#1B9E77" ,"#D95F02", "#7570B3" ,"#E7298A" ,"#66A61E" ,"#E6AB02" ,"#A6761D", "#666666", "#8DD3C7" ,"#FFFFB3", "#BEBADA" ,"#FB8072", "#80B1D3", "#FF0000", "#B3DE69", "#FCCDE5","#641E16"))+
    labs(title=title, x ="Gene", y = "Mutation Burden")
}

generate_plot <- function(mut_type, cancer, clinsig){
  cosmic_germline_for_this_cancer <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Germline.),]
  cosmic_somatic_for_this_cancer <- Cosmic.tissue[grep(cancer, Cosmic.tissue$Tumour.Types.Somatic.),]

  df <- data.frame()
  if(mut_type == "GERMLINE"){
    df <- clinvar_germline_with_enhanced_exome[clinvar_germline_with_enhanced_exome$gene %in% cosmic_germline_for_this_cancer$`Gene.Symbol`,]
    
  } else if(mut_type == "SOMATIC"){
    df <- clinvar_somatic_with_enhanced_exome[clinvar_somatic_with_enhanced_exome$gene %in% cosmic_somatic_for_this_cancer$`Gene.Symbol`,]
  }
  
  print(df)
  
  head(df)
  df$CLNSIG <- toupper(df$CLNSIG)
  df_with_filtered_clinsig <- filter_path(df, clinsig)
  df_with_filtered_clinsig.mb <- get_mb(df_with_filtered_clinsig)
  df_with_filtered_clinsig.mb <- scale_data(df_with_filtered_clinsig.mb)
  tit <- paste(toupper(cancer), mut_type, clinsig)
  graph(df_with_filtered_clinsig.mb, tit)
}

generate_plot("SOMATIC", "breast", "PATHOGENIC")
generate_plot("GERMLINE", "breast", "PATHOGENIC")
generate_plot("GERMLINE", "ovarian", "PATHOGENIC")
generate_plot("GERMLINE", "endomet", "PATHOGENIC")
generate_plot("GERMLINE", "prostate", "PATHOGENIC")
generate_plot("GERMLINE", "thyroid", "PATHOGENIC")

generate_plot("GERMLINE", "breast", "UNCERTAIN")
generate_plot("GERMLINE", "ovarian", "UNCERTAIN")
generate_plot("GERMLINE", "endomet", "UNCERTAIN")
generate_plot("GERMLINE", "prostate", "UNCERTAIN")
generate_plot("GERMLINE", "thyroid", "UNCERTAIN")


head(clinvar_germline_with_enhanced_exome)

cosmic_germline_breast <- Cosmic.tissue[grep("breast", Cosmic.tissue$Tumour.Types.Germline.),]
breast_snps <- clinvar_germline_with_enhanced_exome[clinvar_germline_with_enhanced_exome$gene %in% cosmic_germline_breast$`Gene.Symbol`,]
head(breast_snps)
dim(breast_snps)

data_long <- gather(breast_snps, ancestry, EAF, nfe_seu:oth, factor_key = TRUE)

head(data_long)


data_long$ancestry <- factor(data_long$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn",
                                                            "eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))

data_long


e <- ggplot(data = data_long, aes(x=EAF, fill = ancestry, color = ancestry)) + geom_density(alpha = 0.7) + theme_minimal()
e+scale_fill_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                             "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F"
)) + scale_color_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                                 "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F"
)) + scale_x_continuous(limits = c(0.45, .5)) + scale_y_continuous(limits = c(0, 1000)) 

e


cosmic_somatic_breast <- Cosmic.tissue[grep("breast", Cosmic.tissue$Tumour.Types.Somatic.),]
breast_snps_2 <- clinvar_somatic_with_enhanced_exome[clinvar_somatic_with_enhanced_exome$gene %in% cosmic_somatic_breast$`Gene.Symbol`,]
head(breast_snps_2)
dim(breast_snps_2)

data_long_2 <- gather(breast_snps_2, ancestry, EAF, nfe_seu:oth, factor_key = TRUE)

head(data_long_2)


data_long_2$ancestry <- factor(data_long_2$ancestry, levels = c("afr", "amr", "asj","eas", "eas_jpn",
                                                            "eas_kor" ,"eas_oea", "fin", "nfe", "nfe_bgr", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "nfe_swe","oth", "sas"))

f <- ggplot(data = data_long_2, aes(x=EAF, fill = ancestry, color = ancestry)) + geom_density(alpha = 0.7) + theme_minimal()
f+scale_fill_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                             "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F"
)) + scale_color_manual(values=c("#80B1D3","#F781BF","#B2ABD2","#005A32","#238B45", "#41AB5D","#74C476","#CE1256",
                                 "#67000D", "#A50F15","#CB181D", "#EF3B2C","#FB6A4A","#FC9272","#FCBBA1" ,"#BF812D","#FFD92F"
)) + scale_x_continuous(limits = c(.1, .4)) + scale_y_continuous(limits = c(0, 1000)) 

f



