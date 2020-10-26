################################################################
##        REMOVE SEQUENCES THAT ARE NOT IN METADATA 
##
##  date: 26 October 2020
###############################################################




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
library(seqinr)
library("argparse")


## Preparing environment
rm(list=ls())

setwd('~/Documents/RabiesPractice/data')

parser <- ArgumentParser()

parser$add_argument("-f", "--fasta", type="character",
                    help="sequences in fasta format")
parser$add_argument("-m", "--meta", type="character",
                    help="meta data, in tab separated format")
parser$add_argument("-o", "--output_fasta", type="character",
                    help="path to output fasta with sequences removed")

args <- parser$parse_args()

if (is.null(args$fasta) || is.null(args$meta)) {
  arg_parser$print_help()
  stop("You must specify the path for all three files: fasta, meta and output_fasta", call.=FALSE)
}

##WITH PARSER
#meta = read.table(args$meta, header = TRUE)
#seq = read.fasta(args$fasta)

##HARD CODE
meta <- read.delim("/Volumes/@home/rabies/data/ncbi_cleaned.tab")
seq <- read.fasta('/Volumes/@home/rabies/data/ncbi_RABV_11500_12000.fasta')

names(seq) <- gsub("\\..*","",names(seq))

seq <- seq[names(seq) %in% meta$Accession]

##HARD CODE
#write.fasta(seq, names = names(seq), file = 'rabies2.fasta')

#WITH PARSER
write.fasta(seq, names = names(seq), file = args$output_fasta)








