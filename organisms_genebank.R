################################################################
##        Downloading Sequences and MetaData for Genomes 
##                  of certain ORGANISM and SIZE
##
##  date: 21 October 2020
###############################################################
  

load.lib<-c("ggplot2","lubridate","pander","dplyr","tidyverse","stringr",
            "shiny", "leaflet", "tidyr", "reutils", "BiocManager","ape",
            "rentrez")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

BiocManager::install()
BiocManager::install(c("IRanges", "Biostrings"))
devtools::install_github("ropensci/rentrez")

## Installing BiocManager and BiocGenerics for the R version 3.6

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

devtools::install_github("gschofl/biofiles")
BiocManager::install()
BiocManager::install("IRanges")
BiocManager::install("Biostrings")

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(reutils)
library(ape)
library(rentrez)
library(reutils)
library(biofiles)
library(Biostrings)


## Preparing environment
rm(list=ls())

## Some invariant parameters-------------------------------
# Feature terms to be retrieved
ft_terms <- c("accession", "isolate", "host", "db_xref", "country", "collection_date", "note")

# Function to retrieve the features
features_func <- function(features_list, term) {
  
  out <- rep(NA, length(features_list))
  
  for(i in 1:length(features_list)){
    search <- paste(term," = ", sep = "")
    start = regexpr(search, features_list[[i]])[1] + nchar(search)
    
    if((start - nchar(search)) > 0){
      str <- NA
      str <- substr(features_list[[i]], start, nchar(features_list[[i]]))
      str <- gsub("))", ",", str)
      
      end <- regexpr(",", str)[1] - 1
      out[i] <- substr(str, 1, end)
    }
  }
  
  return(out)
}


## Retrieve GenBank data-----------------------------------

#Define search
#rab_search = esearch(term = "Rabies lyssavirus[ORGN] AND (Algeria OR Morocco OR Tunisia OR Ceuta OR Melilla 
#              OR Western Sahara OR Egypt OR Sudan OR Libya)", db = "nucleotide", retmax = 1000, usehistory = TRUE)

rab_search11920 = esearch(term = "Rabies lyssavirus[ORGN] AND 11920:12000[Sequence Length]", db = "nucleotide", retmax = 3000, usehistory = TRUE)
rab_search11500 = esearch(term = "Rabies lyssavirus[ORGN] AND 11500:11919[Sequence Length]", db = "nucleotide", retmax = 3000, usehistory = TRUE)

#retrieve the search from GenBank
efetch(uid = rab_search11920, db="nucleotide", rettype = "gb", outfile = "rab11920.gb")
efetch(uid = rab_search11500, db="nucleotide", rettype = "gb", outfile = "rab11500.gb")

# Convert GenBank files to gbRecord
rab11920_gb <- biofiles::gbRecord("rab11920.gb")
rab11500_gb <- biofiles::gbRecord("rab11500.gb")


## Create the features table-------------------------------
# Create the features list
ft_list_raw <- getFeatures(rab11920_gb)
ft_list <- vector(mode = "list", length = length(rab11920_gb))

for (i in 1:length(rab11920_gb)) {
  ft <- ft_list_raw[[i]][1]
  ft_list[[i]] <- gsub("\"", "", ft)
}
# ft_list contains the features of each sequence


# Create the frame of the features table  
ft_table <- sapply(ft_terms, function(x,y) features_func(x, y), x=ft_list)
ft_table <- as.data.frame(ft_table)

# Suppress "taxon:" in the db_xref column
ft_table$db_xref <- as.numeric(sub("taxon:", "", ft_table$db_xref))
summary(ft_table$db_xref) # We have only Rabies lyssaviruses, no mistake was made

# Complete the data frame with definition, sequence length and sequence version
ft_table$length <- biofiles::getLength(rab11920_gb)
ft_table$definition <- getDefinition(rab11920_gb)
ft_table$version <- getVersion(rab11920_gb)

# Get the sequences from the data file
seq <- biofiles::getSequence(rab11920_gb)

# Write the sequences to a Fasta file
writeXStringSet(seq, file = 'rabies11920.fasta')

# Write the features meta data table to a csv
write.csv(ft_table, file = "metadata11920.csv")


#FASTA Files were then combined together for both 11920 and 11500 datasets
#Metadata tables were also combined together and saved as a whole