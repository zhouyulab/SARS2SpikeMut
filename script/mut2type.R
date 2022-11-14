library(readxl)
library(argparse)
library(dplyr)
library(readr)
library(Biostrings)

parser <- ArgumentParser()
parser$add_argument("--mut-pos-info", required=TRUE, type="character", dest = "mut_pos_info", metavar="MutPosInfo.xlsx")
parser$add_argument("--ref-fa", required=TRUE, type="character", dest = "ref_fa", metavar="ref_fa.fa")
parser$add_argument("--prot-mut", required=TRUE, type="character", dest = "prot_mut", metavar="prot_mut.tsv")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output.tsv")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)


# args <- c(
#   "--mut-pos-info", "~/mu01/project/cov2mut/data/MutPosInfo.xlsx",
#   "--ref-fa", "~/mu01/project/cov2mut/data/prot/spike.ref.fasta",
#   "--prot-mut", "~/mu01/project/cov2mut/analysis/seq_info/uniq_prot.mut.tsv",
#   "-o", "~/mu01/project/cov2mut/analysis/seq_info/uniq_prot.mut_type.tsv"
# )
# args <- parser$parse_args(args)

MutPosInfo <- read_excel(args$mut_pos_info)
MutPosInfo <- MutPosInfo[!duplicated(MutPosInfo),]
MutNumInfo <- MutPosInfo %>% group_by(MutName) %>% summarise(MutNum=n())
ref_fa <- readBStringSet(args$ref_fa)
ref_prot_seq <- strsplit(as.vector(ref_fa), "")[[1]]

ref_prot_seq_df <- data.frame(Position=1:length(ref_prot_seq), OriAa=ref_prot_seq, IsOriAaMatch=TRUE)
MutPosInfo <- left_join(MutPosInfo, ref_prot_seq_df)
MutPosInfo$IsOriAaMatch[is.na(MutPosInfo$IsOriAaMatch)] <- FALSE
assertthat::assert_that(all(MutPosInfo$IsOriAaMatch))
MutPosInfo$IsOriAaMatch <- NULL

is_contain <- function(A, B, MutPosInfo){
  A_df <- MutPosInfo[as.character(MutPosInfo$MutName)==A,]
  B_df <- MutPosInfo[as.character(MutPosInfo$MutName)==B,]
  A_num <- nrow(A_df)
  B_num <- nrow(B_df)
  if(A_num < B_num){
    return(FALSE)
  }
  both_df <- inner_join(A_df[,c("Position", "MutAa")], B_df[,c("Position", "MutAa")])
  both_num <- nrow(both_df)
  return(both_num == B_num)
}

find_mut_type_relationship <- function(MutPosInfo){
  uniq_mut_type <- unique(MutPosInfo$MutName)
  big_set <- c()
  small_set <- c()
  for(mut1 in uniq_mut_type){
    for(mut2 in uniq_mut_type){
      if(mut1 == mut2) next()
      if(is_contain(mut1, mut2, MutPosInfo)){
        big_set <- c(big_set, mut1)
        small_set <- c(small_set, mut2)
      }
    }
  }
  df <- data.frame(BigSet=big_set, SmallSet=small_set)
  return(df)
}

mut_type_relationship <- find_mut_type_relationship(MutPosInfo)

prot_mut_df <- read_delim(args$prot_mut, "\t", escape_double = FALSE, trim_ws = TRUE)

ann_prot_mut_df <- inner_join(prot_mut_df, MutPosInfo)
ann_prot_mut_num_df <- ann_prot_mut_df %>% group_by(UniqProtID, MutName) %>% summarise(SeqMutNum=n())
ann_prot_mut_num_df <- inner_join(ann_prot_mut_num_df, MutNumInfo)
protential_prot_mut_label_df <- ann_prot_mut_num_df[ann_prot_mut_num_df$SeqMutNum==ann_prot_mut_num_df$MutNum,]

prot_mut_label_df <- plyr::ddply(protential_prot_mut_label_df, "UniqProtID", function(tmp_mut_num_df, mut_type_relationship){
  tmp_mut_label <- tmp_mut_num_df$MutName
  tmp_big_set <- intersect(tmp_mut_label, mut_type_relationship$BigSet)
  for(big_label in tmp_big_set){
    small_label_set <- mut_type_relationship$SmallSet[mut_type_relationship$BigSet==big_label]
    tmp_mut_label <- tmp_mut_label[!tmp_mut_label %in% small_label_set]
  }
  res <- data.frame(UniqProtID=tmp_mut_num_df$UniqProtID[1], MutName=tmp_mut_label)
  return(res)
}, mut_type_relationship=mut_type_relationship)
write_tsv(prot_mut_label_df, args$output)