library(R6)
library(data.table)

if(!exists("disableBiom")){
  library(biom)
}

source("../common/mbiome/mbiome-data.R")
source("../common/mbiome/mbiome-utils.R")
source("../common/mbiome/otu-table.R")
source("../common/mbiome/sample-table.R")


import.biom <- function(biom_path, metadata_details=NULL, aggregate_by="Phylum", use_relative_abundance=T){

  process_column <- function(smpl_id, vector){
    line<-vector[vector>0]
    otus<-names(line)
    new_row<-data.table(OtuId=otus, AbsoluteAbundance=line, Sample=smpl_id)
    new_row
  }
  
  biom_obj<-read_biom(biom_path)
  abundance_data<-biom_data(biom_obj)
  sample_ids <- colnames(abundance_data)
  main_dt<-data.table()
  for(i in 1:ncol(abundance_data)){
    smpl<-sample_ids[i]
    main_dt<-rbindlist(list(main_dt, process_column(smpl, abundance_data[,i]) ) )
  }
  otus<-observation_metadata(biom_obj)
  otu_dt<-as.data.table(otus)
  colnames(otu_dt)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # fix the names of the species here
  otu_dt[,OtuId:=rownames(otus)]
  
  columns_starters<-c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  gg_starters<-c("s__", "g__", "f__", "o__", "c__", "p__", "k__")
  
  for(j in 1:7){
    set(otu_dt, j=columns_starters[j], value=gsub(sprintf("^%s(.*)$", gg_starters[j]), "\\1", otu_dt[[columns_starters[j]]]))
    set(otu_dt, i=which(otu_dt[[columns_starters[j]]]==""), j=columns_starters[j], value="N/A")
  }
  
  complete_otus<-merge(main_dt, otu_dt)
  otu_object <- OTUClass$new(complete_otus, aggregate_by= "Phylum", use_relative_abundance=F)
  if(use_relative_abundance){
    otu_object$calculate_relative_abundance() 
  }
  
  sample_df<-sample_metadata(biom_obj)
  sample_dt<-as.data.table(sample_df)
  sample_dt[,SampleName:=rownames(sample_df)]
  
  if(!is.null(metadata_details)){
    metadata_details<-fread(metadata_details)
  }
  sample_object <- SampleClass$new(biom_dt=sample_dt, details_dt=metadata_details)
  
  MicrobiomeData$new(otu_object, sample_object)
  
}

#' Function to read input files for MicrobiomeDB and returns a MicrobiomeData
#'   object.
#' @param taxa_abundance_path The path to the taxa abundance file.
#' @param sample_path The path to the sample characteristics file.
#' @param aggregate_by Taxonomic level to aggregate your data (Phylum, Class, Order,
#' Family, Genus, Species). Defaults to Phylum.
#' @param use_relative_abundance If it's true the relative abundance will be used to manipulate
#'   the data, otherwise the absolute abundance will be used. Defaults to TRUE.
#' @return Returns a MicrobiomeData object.
#' @keywords microbiome microbiomedb
#' @export
#' @example
#' import.eupath("~/Projects/microbiomedb/TaxaRelativeAbundance.txt",
#'     "~/Projects/microbiomedb/Characteristics.txt")
import.eupath <- function(taxa_abundance_path, sample_path, aggregate_by="Phylum", use_relative_abundance=T){
  df_abundance <-
    fread(
      taxa_abundance_path,
      col.names = c("Sample","Taxon", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "RelativeAbundance", "AbsoluteAbundance", "EmptyColumn"),
      colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character", "character", "numeric", "integer", "character")
    )

  df_sample <-
    fread(
      sample_path,
      col.names = c("SampleName", "Source", "Property", "Value", "Type", "Filter", "EmptyColumn"),
      colClasses = c("character", "character", "character", "character", "character", "character", "character")
    )
  # Removing the last unuseless column
  
  df_abundance[,EmptyColumn:=NULL]
  
  
  otu_object <- OTUClass$new(df_abundance, aggregate_by = aggregate_by, use_relative_abundance=use_relative_abundance)
  
  sample_object <- SampleClass$new(df_sample)
  MicrobiomeData$new(otu_object, sample_object)
}
