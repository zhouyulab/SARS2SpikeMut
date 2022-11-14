library(readxl)
library(argparse)
library(dplyr)
library(readr)
library(Biostrings)
library(BiocParallel)
library(ComplexHeatmap)
library(ggplot2)

parser <- ArgumentParser()
parser$add_argument("--variant-info", required=TRUE, type="character", dest = "mut_pos_info", metavar="VariantInfo.tsv")
parser$add_argument("--ref-fa", required=TRUE, type="character", dest = "ref_fa", metavar="ref_fa.fa")
parser$add_argument("--prot-mut", required=TRUE, type="character", dest = "prot_mut", metavar="prot_mut.tsv")
parser$add_argument("--seq-info", required=TRUE, type="character", dest = "seq_info", metavar="seq_info.tsv")
parser$add_argument("--seq2uniq", required=TRUE, type="character", dest = "seq2uniq", metavar="seq2uniq.tsv")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output.tsv")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)

MutPosInfo <- MutPosInfo[!duplicated(MutPosInfo),]
MutNumInfo <- MutPosInfo %>% group_by(MutName) %>% summarise(MutNum=n())
ref_fa <- readBStringSet(args$ref_fa)
ref_prot_seq <- strsplit(as.vector(ref_fa), "")[[1]]

ref_prot_seq_df <- data.frame(Position=1:length(ref_prot_seq), OriAa=ref_prot_seq, IsOriAaMatch=TRUE)
MutPosInfo <- left_join(MutPosInfo, ref_prot_seq_df)
MutPosInfo$IsOriAaMatch[is.na(MutPosInfo$IsOriAaMatch)] <- FALSE
assertthat::assert_that(all(MutPosInfo$IsOriAaMatch))
MutPosInfo$IsOriAaMatch <- NULL

prot_mut_df <- read_delim(args$prot_mut, "\t", escape_double = FALSE, trim_ws = TRUE)

ann_pos_df <- inner_join(prot_mut_df, MutPosInfo[,c("MutName", "Position")])
mut_with_X_df <- ann_pos_df[ann_pos_df$MutAa=="X",]
mut_with_X_df <- mut_with_X_df %>% group_by(UniqProtID, MutName) %>% summarise(WithX=TRUE)

ann_prot_mut_df <- inner_join(prot_mut_df, MutPosInfo)
## Only sequence with D614G will be used
D614G_ID <- unique(ann_prot_mut_df$UniqProtID[ann_prot_mut_df$Position==614 & ann_prot_mut_df$MutAa=="G"])
ann_prot_mut_df <- ann_prot_mut_df[ann_prot_mut_df$UniqProtID %in% D614G_ID,]

ann_prot_mut_df <- ann_prot_mut_df[order(ann_prot_mut_df$UniqProtID, ann_prot_mut_df$MutName, ann_prot_mut_df$Position),]
mut_set_df <- ann_prot_mut_df %>% group_by(UniqProtID, MutName) %>% summarise(MutSetName=paste(sprintf("%s%s%s", OriAa, Position, MutAa), collapse = "~"))
uniq_mut_set_per_variant_df <- mut_set_df[!duplicated(mut_set_df[,c("MutName", "MutSetName")]),]

uniq_mut_set_df <- mut_set_df[!duplicated(mut_set_df$MutSetName),]
uniq_mut_set_df <- left_join(uniq_mut_set_df, ann_prot_mut_df)
uniq_mut_set_df <- uniq_mut_set_df[,c("MutSetName", "Position", "OriAa", "MutAa")]
write_tsv(uniq_mut_set_df, file.path(args$output, "UniqMutSet.tsv"))

uniq_mut_set_per_variant_df <- uniq_mut_set_per_variant_df[,c("MutName", "MutSetName")]
write_tsv(uniq_mut_set_df, file.path(args$output, "MutSet2Variant.tsv"))


find_mut_type_relationship <- function(uniq_mut_set_df){
  uniq_mut_type <- unique(uniq_mut_set_df$MutSetName)
  mut_num_df <- uniq_mut_set_df %>% group_by(MutSetName) %>% summarise(MutSiteNum=n())
  big_set <- c()
  small_set <- c()
  for(mut1 in uniq_mut_type){
    mut1 <- as.character(mut1)
    mut1_mut_site_num <- mut_num_df$MutSiteNum[mut_num_df$MutSetName==mut1]
    protential_mut2 <- mut_num_df$MutSetName[mut_num_df$MutSiteNum<mut1_mut_site_num]
    if(length(protential_mut2)==0) next()
    mut1_df <- uniq_mut_set_df[uniq_mut_set_df$MutSetName==mut1,]
    mut1_df$MutSetName <- NULL
    protential_mut2_df <- uniq_mut_set_df[uniq_mut_set_df$MutSetName%in%protential_mut2,]
    overlap_mut_pos <- inner_join(mut1_df, protential_mut2_df)
    overlap_mut_pos_num <- overlap_mut_pos %>% group_by(MutSetName) %>% summarise(OverlapMutSiteNum=n())
    overlap_mut_pos_num <- inner_join(overlap_mut_pos_num, mut_num_df)
    mut2_li <- overlap_mut_pos_num$MutSetName[overlap_mut_pos_num$OverlapMutSiteNum==overlap_mut_pos_num$MutSiteNum]
    if(length(mut2_li) == 0) next()
    big_set <- c(big_set, rep(mut1, length(mut2_li)))
    small_set <- c(small_set, mut2_li)
  }
  df <- data.frame(BigSet=big_set, SmallSet=small_set)
  return(df)
}
uniq_mut_set_df$MutSetName <- as.character(uniq_mut_set_df$MutSetName)
mut_type_relationship <- find_mut_type_relationship(uniq_mut_set_df)
write_tsv(mut_type_relationship, file.path(args$output, "MutSetRelationship.tsv"))

uniq_prot_mut_set_df <- mut_set_df[!duplicated(mut_set_df[,c("UniqProtID", "MutSetName")]),]
uniq_prot_mut_set_df$MutName <- NULL
uniq_prot_mut_set_num <- uniq_prot_mut_set_df %>% group_by(UniqProtID) %>% summarise(LabelNum=n())

prot_mut_label_li <- bplapply(unique(uniq_prot_mut_set_df$UniqProtID), function(uniq_id, uniq_prot_mut_set_df, mut_type_relationship){
  tmp_mut_num_df <- uniq_prot_mut_set_df[uniq_prot_mut_set_df$UniqProtID==uniq_id,]
  tmp_mut_label <- tmp_mut_num_df$MutSetName
  tmp_big_set <- intersect(tmp_mut_label, mut_type_relationship$BigSet)
  for(big_label in tmp_big_set){
    small_label_set <- mut_type_relationship$SmallSet[mut_type_relationship$BigSet==big_label]
    tmp_mut_label <- tmp_mut_label[!tmp_mut_label %in% small_label_set]
  }
  res <- data.frame(UniqProtID=tmp_mut_num_df$UniqProtID[1], MutSetName=tmp_mut_label)
  return(res)
}, mut_type_relationship=mut_type_relationship, uniq_prot_mut_set_df=uniq_prot_mut_set_df, BPPARAM=MulticoreParam(10))
prot_mut_label_df <- do.call(rbind, prot_mut_label_li)
prot_mut_label_df <- inner_join(prot_mut_label_df, uniq_mut_set_per_variant_df)
write_tsv(prot_mut_label_df, file.path(args$output, "UniqProtID2label.tsv"))

## Remove sequence with X in mutation set
prot_mut_label_df <- left_join(prot_mut_label_df, mut_with_X_df)
prot_mut_label_df$WithX[is.na(prot_mut_label_df$WithX)] <- FALSE
prot_mut_label_df <- prot_mut_label_df[prot_mut_label_df$WithX==FALSE,]
write_tsv(prot_mut_label_df, file.path(args$output, "UniqProtID2label.remove_X.tsv"))

seq_info_df <- read_delim(args$seq_info, "\t", escape_double = FALSE, trim_ws = TRUE)
seq_info_df <- seq_info_df[,c("IsolateName", "Year", "Month", "Day", "State", "Location")]
seq_info_df <- seq_info_df[!duplicated(seq_info_df),]
seq_info_iso_num <- seq_info_df %>% group_by(IsolateName) %>% summarise(Num=n())
seq_info_df <- seq_info_df[seq_info_df$IsolateName %in% seq_info_iso_num$IsolateName[seq_info_iso_num$Num==1],]

seq2uniq_df <- read_delim(args$seq2uniq, "\t", escape_double = FALSE, trim_ws = TRUE)
seq2uniq_df <- seq2uniq_df[!duplicated(seq2uniq_df),]
seq2uniq_iso_num <- seq2uniq_df %>% group_by(IsolateName) %>% summarise(Num=n())
seq2uniq_df <- seq2uniq_df[seq2uniq_df$IsolateName %in% seq2uniq_iso_num$IsolateName[seq2uniq_iso_num$Num==1],]

seq2uniq_info_df <- left_join(seq2uniq_df, seq_info_df)

seq2label_df <- left_join(seq2uniq_info_df, prot_mut_label_df)
seq2label_df <- na.omit(seq2label_df)
write_tsv(seq2label_df, file.path(args$output, "seq2label.tsv"))

mut_time_num <- seq2label_df %>% group_by(MutName, MutSetName, Year, Month) %>% summarise(IsolateNum=n())
mut_time_num <- mut_time_num[mut_time_num$Year>=2020,]
mut_time_num <- mut_time_num[mut_time_num$Month!="00",]
write_tsv(mut_time_num, file.path(args$output, "MutTimeNum.tsv"))
