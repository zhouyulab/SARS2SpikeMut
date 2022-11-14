library(readxl)
library(argparse)
library(dplyr)
library(readr)
library(Biostrings)
library(BiocParallel)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(plotrix)

parser <- ArgumentParser()
parser$add_argument("--mut-time", required=TRUE, type="character", dest = "mut_time", metavar="mut_time.tsv")
parser$add_argument("--uniq-mut-set", required=TRUE, type="character", dest = "uniq_mut_set", metavar="uniq_mut_set.tsv")
parser$add_argument("--seq2label", required=TRUE, type="character", dest = "seq2label", metavar="seq2label.tsv")
parser$add_argument("--seq2uniq", required=TRUE, type="character", dest = "seq2uniq", metavar="seq2uniq.tsv")
parser$add_argument("--seq-info", required=TRUE, type="character", dest = "seq_info", metavar="seq_info.tsv")
parser$add_argument("--mut-pos-info", required=TRUE, type="character", dest = "mut_pos_info", metavar="mut_pos_info.tsv")
parser$add_argument("--expr-mut-pos-info", required=TRUE, type="character", dest = "expr_mut_pos_info", metavar="expr_mut_pos_info.xlsx")
parser$add_argument("--mut-expr", required=TRUE, type="character", dest = "mut_expr", metavar="mut_expr.xlsx")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)

#args <- c(
#  "--mut-time", "~/mu01/project/cov2mut/analysis/variant_stat/MutTimeNum.tsv",
#  "--uniq-mut-set", "~/mu01/project/cov2mut/analysis/variant_stat/UniqMutSet.tsv",
#  "--seq2label", "~/mu01/project/cov2mut/analysis/variant_stat/seq2label.tsv",
#  "--mut-expr", "~/mu01/project/cov2mut/data/MutExprData.xlsx",
#  "--mut-pos-info", "~/mu01/project/cov2mut/data/VariantInfo.tsv",
#  "--expr-mut-pos-info", "~/mu01/project/cov2mut/data/MutInfo.xlsx",
#  "--seq2uniq", "~/mu01/project/cov2mut/analysis/seq_info/seq2uniq.tsv",
#  "--seq-info", "~/mu01/project/cov2mut/analysis/seq_info/seq_info.tsv",
#  "-o", "~/mu01/project/cov2mut/analysis/variant_stat"
#)
#args <- parser$parse_args(args)

seq_info_df <- read_delim(args$seq_info, "\t", escape_double = FALSE, trim_ws = TRUE)
seq_info_df <- seq_info_df[,c("IsolateName", "Year", "Month", "Day", "State", "Location")]
seq_info_df <- seq_info_df[!duplicated(seq_info_df),]
seq2uniq_df <- read_delim(args$seq2uniq, "\t", escape_double = FALSE, trim_ws = TRUE)
seq2uniq_df <- seq2uniq_df[!duplicated(seq2uniq_df),]
WT_df <- seq2uniq_df[seq2uniq_df$UniqProtID==seq2uniq_df$UniqProtID[seq2uniq_df$IsolateName=="hCoV-19/Wuhan/WIV04/2019"],]
WT_df <- left_join(WT_df, seq_info_df)
WT_df <- WT_df[WT_df$Month!="00" & WT_df$Day!="00",]
WT_date_num <- WT_df %>%group_by(Year, Month) %>% summarise(MutSetName="WT", IsolateNumPerMonth=n())
WT_date_num$Time <- sprintf("%s-%s", WT_date_num$Year, WT_date_num$Month)

mut_time_num <- read_delim(args$mut_time, "\t", escape_double = FALSE, trim_ws = TRUE)
uniq_mut_set_df <- read_delim(args$uniq_mut_set, "\t", escape_double = FALSE, trim_ws = TRUE)
seq2label_df <- read_delim(args$seq2label, "\t", escape_double = FALSE, trim_ws = TRUE)
mut_expr <- read_excel(args$mut_expr)
mut_expr <- mut_expr[!duplicated(mut_expr),]
mut_expr$Infectivity <- as.numeric(mut_expr$Infectivity)
expr_mut_pos_info <- read_excel(args$expr_mut_pos_info)
expr_mut_pos_info <- expr_mut_pos_info[!duplicated(expr_mut_pos_info),]
names(expr_mut_pos_info)[1] <- "ExprName"
names(mut_expr)[1] <- "ExprName"
expr_mut_num <- expr_mut_pos_info %>% group_by(ExprName) %>% summarise(ExprMutNum=n())

all(expr_mut_pos_info$ExprName %in% mut_expr$ExprName)
all(mut_expr$ExprName %in% expr_mut_pos_info$ExprName)
mut_expr$ExprName[!mut_expr$ExprName %in% expr_mut_pos_info$ExprName]

expr_mut_pos_info$ExprName[!expr_mut_pos_info$ExprName %in% mut_expr$ExprName]
expr_mut_pos_info <- expr_mut_pos_info[expr_mut_pos_info$ExprName %in% mut_expr$ExprName,]

MutPosInfo <- read_delim(args$mut_pos_info, "\t", escape_double = FALSE, trim_ws = TRUE)
MutPosInfo <- MutPosInfo[!duplicated(MutPosInfo),]

expr_mut_pos_info <- expr_mut_pos_info[order(expr_mut_pos_info$ExprName, expr_mut_pos_info$Position),]
expr2mut_set_df <- expr_mut_pos_info %>% group_by(ExprName) %>% summarise(MutSetName=paste(sprintf("%s%s%s", OriAa, Position, MutAa), collapse = "~"))

mut_time_num$Time <- sprintf("%s-%s", mut_time_num$Year, mut_time_num$Month)

seq2label_df$Date <- sprintf("%s-%s-%s", seq2label_df$Year, seq2label_df$Month, seq2label_df$Day)
seq2label_df <- seq2label_df[seq2label_df$Year>=2020,]
seq2label_df <- seq2label_df[seq2label_df$Month!="00",]
seq2label_df <- seq2label_df[seq2label_df$Day!="00",]
seq2label_df <- seq2label_df[order(seq2label_df$MutSetName ,seq2label_df$Year, seq2label_df$Month, seq2label_df$Day),]

label2num_df <- seq2label_df %>% group_by(MutSetName, Year, Month) %>% summarise(IsolateNumPerMonth=n())
label2num_df <- label2num_df[order(label2num_df$Year, label2num_df$Month),]
label2num_df$Time <- sprintf("%s-%s", label2num_df$Year, label2num_df$Month)

first_time_isolate <- seq2label_df[!duplicated(seq2label_df[,c("MutSetName")]),]
first_time_isolate <- first_time_isolate[,c("MutSetName","IsolateName", "Date", "State", "Location")]

min_time <- min(as.Date(seq2label_df$Date, format = "%Y-%m-%d"))
max_time <- max(as.Date(seq2label_df$Date, format = "%Y-%m-%d"))

time_li <- sort(unique(mut_time_num$Time))
mut_time_num$Time <- factor(mut_time_num$Time, levels = time_li)

total_time_num <- seq2label_df[!duplicated(seq2label_df$IsolateName),] %>% group_by(Year, Month) %>% summarise(IsolatePerTime=n())



for(mut_name in unique(mut_time_num$MutName)){
  tmp_mut_time_num <- mut_time_num[mut_time_num$MutName==mut_name,]
  tmp_mut_set_name_df <- tmp_mut_time_num %>% group_by(MutName, MutSetName) %>% summarise()
  tmp_mut_set_pos_df <- left_join(tmp_mut_set_name_df, uniq_mut_set_df)
  tmp_mut_set_pos_df$PosName <- sprintf("%s%s%s", tmp_mut_set_pos_df$OriAa, tmp_mut_set_pos_df$Position, tmp_mut_set_pos_df$MutAa)
  tmp_mut_pos_df <- MutPosInfo[MutPosInfo$MutName==mut_name,]
  tmp_mut_pos_df$PosName <- sprintf("%s%s%s", tmp_mut_pos_df$OriAa, tmp_mut_pos_df$Position, tmp_mut_pos_df$MutAa)
  tmp_mut_time_num_info <- left_join(tmp_mut_time_num, tmp_mut_set_pos_df)
  mut_pos_mat <- matrix("F", nrow = nrow(tmp_mut_set_name_df), ncol = nrow(tmp_mut_pos_df))
  rownames(mut_pos_mat) <- tmp_mut_set_name_df$MutSetName
  colnames(mut_pos_mat) <- tmp_mut_pos_df$PosName
  
  for(indx in 1:nrow(tmp_mut_time_num_info)){
    pos_indx <- which(tmp_mut_time_num_info$PosName[indx] == tmp_mut_pos_df$PosName)
    mut_indx <- which(tmp_mut_time_num_info$MutSetName[indx] == tmp_mut_set_name_df$MutSetName)
    if(tmp_mut_time_num_info$MutSetName[indx] %in% expr2mut_set_df$MutSetName){
      label <- "E"
    }else{
      label <- "T"
    }
    mut_pos_mat[mut_indx, pos_indx] <- label
  }

  tmp_mut_time_num_stat <- tmp_mut_time_num %>% group_by(MutSetName) %>% summarise(TotalIsolateNum=sum(IsolateNum), FirstTime=sort(Time)[1])
  tmp_mut_time_num_stat <- left_join(tmp_mut_time_num_stat, first_time_isolate)
  tmp_mut_time_num_stat <- left_join(tmp_mut_set_name_df,tmp_mut_time_num_stat)
  tmp_mut_time_num_stat$FirstTime <- factor(tmp_mut_time_num_stat$FirstTime, levels = time_li)
  mut_order <- order(tmp_mut_time_num_stat$FirstTime, tmp_mut_time_num_stat$Date)
  ht <- Heatmap(
    mut_pos_mat,
    name="Mutation",
    col = c("T"="black", "F"="white", "E"="red"),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    cluster_columns = FALSE, 
    show_row_names = FALSE, 
    row_order=mut_order,
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  ha <- rowAnnotation(
    FirstTime = anno_points(
      as.integer(tmp_mut_time_num_stat$FirstTime), 
      size = unit(1, "mm"),
      width = unit(0.8, "cm"), 
      ylim = c(1, length(time_li)),
      annotation_name_gp=gpar(fontsize = 5)
    ),
    Log10TotalIsolateNum = anno_barplot(
      log10(tmp_mut_time_num_stat$TotalIsolateNum), 
      width = unit(0.8, "cm"), 
      annotation_name_gp=gpar(fontsize = 5)
    ),
    gp=gpar(fontsize = 5),
    annotation_name_gp =gpar(fontsize = 5),
    annotation_legend_param=list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  
  inche_cm <- 2.54
  pdf(file.path(args$output, sprintf("%s.all_mut.pdf", mut_name)), width=10/inche_cm, height=20/inche_cm)
  print(ht + ha)
  dev.off()
  
  expr_mut_pos_mat <- mut_pos_mat[tmp_mut_time_num_stat$TotalIsolateNum>10,]
  expr_mut_time_num_stat <- tmp_mut_time_num_stat[tmp_mut_time_num_stat$TotalIsolateNum>10,]
  expr_mut_order <- order(expr_mut_time_num_stat$FirstTime, expr_mut_time_num_stat$Date)
  expr_ht <- Heatmap(
    expr_mut_pos_mat,
    name="Mutation",
    col = c("T"="black", "F"="white", "E"="red"),
    row_names_gp = gpar(fontsize = 6),
    row_title_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    row_dend_width = unit(5, "mm"),
    column_dend_height = unit(5, "mm"),
    cluster_columns = FALSE, 
    show_row_names = FALSE, 
    row_order=expr_mut_order,
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  expr_ha <- rowAnnotation(
    FirstTime = anno_points(
      as.integer(expr_mut_time_num_stat$FirstTime), 
      size = unit(1, "mm"),
      width = unit(0.8, "cm"), 
      ylim = c(1, length(time_li)),
      annotation_name_gp=gpar(fontsize = 5)
    ),
    Log10TotalIsolateNum = anno_barplot(
      log10(expr_mut_time_num_stat$TotalIsolateNum), 
      width = unit(0.8, "cm"), 
      annotation_name_gp=gpar(fontsize = 5)
    ),
    gp=gpar(fontsize = 5),
    annotation_name_gp =gpar(fontsize = 5),
    annotation_legend_param=list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      grid_width = unit(2, "mm"),
      grid_height = unit(2, "mm")
    )
  )
  
  inche_cm <- 2.54
  pdf(file.path(args$output, sprintf("%s.expr_mut.pdf", mut_name)), width=10/inche_cm, height=10/inche_cm)
  print(expr_ht + expr_ha)
  dev.off()
  
  expr_time_num <- tmp_mut_time_num[tmp_mut_time_num$MutSetName %in% expr_mut_time_num_stat$MutSetName,]
  expr_time_num <- left_join(expr_time_num, total_time_num)
  expr_time_num$IsolateRatio <- expr_time_num$IsolateNum / expr_time_num$IsolatePerTime
  expr_mut_time_num_stat$Label <- sprintf("#Isolates: %s\nFirst time: %s\nFirst isolate:%s", expr_mut_time_num_stat$TotalIsolateNum, expr_mut_time_num_stat$Date, expr_mut_time_num_stat$IsolateName)
  p <-  ggplot(expr_time_num, aes(x=Time, y=IsolateNum)) +
    geom_bar(stat="identity", fill="black") +
    geom_text(data=expr_mut_time_num_stat, mapping=aes(x=0, y=0, label=Label), size=1.2, color="black", hjust=0, vjust=0) +
    facet_wrap(~MutSetName, ncol = 6, scales = "free_y") +
    labs(y="#Total isolates", x="Time") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size = 6),
          title = element_text(family="ArialMT", size = 6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid = element_blank())
  ggsave(file.path(args$output, sprintf("%s.MutTypeTimeNum.pdf", mut_name)), p, height = 2+2*ceiling(length(unique(expr_time_num$MutSetName))/6), width = 25, units = "cm", colormodel = "cmyk")
  
  
  p <-  ggplot(expr_time_num, aes(x=Time, y=IsolateRatio)) +
    geom_bar(stat="identity", fill="black") +
    geom_text(data=expr_mut_time_num_stat, mapping=aes(x=0, y=0, label=Label), size=1.2, color="black", hjust=0, vjust=0) +
    facet_wrap(~MutSetName, ncol = 6, scales = "free_y") +
    labs(y="Isolate ratio", x="Time") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size = 6),
          title = element_text(family="ArialMT", size = 6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid = element_blank())
  ggsave(file.path(args$output, sprintf("%s.MutTypeTimeRatio.pdf", mut_name)), p, height = 2+2*ceiling(length(unique(expr_time_num$MutSetName))/6), width = 25, units = "cm", colormodel = "cmyk")
  
  p <-  ggplot(expr_time_num, aes(x=Time, y=IsolateRatio)) +
    geom_bar(stat="identity", fill="black") +
    geom_text(data=expr_mut_time_num_stat, mapping=aes(x=0, y=0, label=Label), size=1.2, color="black", hjust=0, vjust=0) +
    facet_wrap(~MutSetName, ncol = 6) +
    labs(y="Isolate ratio", x="Time") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size = 6),
          title = element_text(family="ArialMT", size = 6),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid = element_blank())
  ggsave(file.path(args$output, sprintf("%s.MutTypeTimeRatio.SameScale.pdf", mut_name)), p, height = 2+2*ceiling(length(unique(expr_time_num$MutSetName))/6), width = 25, units = "cm", colormodel = "cmyk")
  
  full_mut_set <- MutPosInfo[MutPosInfo$MutName==mut_name,]
  expr_overlap_mutation <- inner_join(full_mut_set, expr_mut_pos_info)
  expr_overlap_mutation_num <- expr_overlap_mutation %>% group_by(ExprName) %>% summarise(MutNum=n())
  expr_overlap_mutation_num <- inner_join(expr_mut_num, expr_overlap_mutation_num)
  match_expr_name <- expr_overlap_mutation_num$ExprName[expr_overlap_mutation_num$ExprMutNum==expr_overlap_mutation_num$MutNum]
  tmp_expr_info_df <- mut_expr[mut_expr$ExprName %in% match_expr_name,]
  tmp_expr_info_df <- left_join(tmp_expr_info_df, expr2mut_set_df)
  tmp_expr_info_df <- left_join(tmp_expr_info_df, tmp_mut_time_num_stat)
  tmp_expr_info_df$IsFound <- tmp_expr_info_df$MutSetName %in% tmp_mut_time_num_stat$MutSetName
  tmp_expr_info_df$Date <- as.Date(tmp_expr_info_df$Date, format = "%Y-%m-%d")
  tmp_expr_info_df$Escape <- as.numeric(tmp_expr_info_df$Escape)
  tmp_expr_info_df$Infectivity <- as.numeric(tmp_expr_info_df$Infectivity)
  tmp_expr_info_found_df <- tmp_expr_info_df[tmp_expr_info_df$IsFound,]
  tmp_expr_info_not_found_df <- tmp_expr_info_df[!tmp_expr_info_df$IsFound,]
  p <- ggplot() +
    geom_point(data=tmp_expr_info_found_df, mapping=aes(x=Escape, y=Infectivity, color=as.numeric(Date), size=log10(TotalIsolateNum)))+
    geom_point(data=tmp_expr_info_not_found_df, mapping=aes(x=Escape, y=Infectivity), color="black", size=2, shape=4) +
    geom_text(data=tmp_expr_info_df, mapping=aes(x=Escape, y=Infectivity, label=ExprName), size=1.2, color="black") +
    scale_color_gradientn(
      colours = rev(c('#FF0000','#FF7F00','#FFFF00','#00FF00','#00FFFF', '#0000FF')),
      breaks = as.numeric(as.Date(c("2020-01-01", "2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01", "2022-07-01"), format = "%Y-%m-%d")),
      labels = c("2020-01-01", "2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01", "2022-07-01"),
      limits=as.numeric(c(as.Date("2020-01-01", format = "%Y-%m-%d"), max_time))
      ) +
    scale_size_continuous(limits = c(0, 6), range = c(0, 4)) +
    scale_x_continuous(limits = c(min(mut_expr$Escape, na.rm = T), max(mut_expr$Escape, na.rm = T))) +
    scale_y_continuous(limits = c(min(mut_expr$Infectivity, na.rm = T), max(mut_expr$Infectivity, na.rm = T))) +
    labs(size="Log10IsolateNum", color="First time") +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size = 6),
          title = element_text(family="ArialMT", size = 6),
          axis.text = element_text(colour = "black"),
          legend.key.size = unit(2, "mm"),
          legend.title = element_text(family="ArialMT", size = 6),
          legend.text = element_text(family="ArialMT", size = 6),
          panel.grid = element_blank())
  ggsave(file.path(args$output, sprintf("%s.Expr.pdf", mut_name)), p, height = 8, width = 10, units = "cm", colormodel = "cmyk")
}



mut_expr <- left_join(mut_expr, expr2mut_set_df)
all_expr_data <- inner_join(mut_expr, first_time_isolate)
all_isolate_num <- seq2label_df[!duplicated(seq2label_df[,c("IsolateName", "MutSetName")]),] %>% group_by(MutSetName) %>% summarise(TotalIsolateNum=n())
all_expr_data <- inner_join(all_expr_data, all_isolate_num)
all_expr_data$Date <- as.Date(all_expr_data$Date, format = "%Y-%m-%d")
all_expr_data <- na.omit(all_expr_data)
all_expr_data$Label <- sapply(strsplit(all_expr_data$ExprName, "[+]"), function(x){return(x[length(x)])})
all_expr_data$Escape <- as.numeric(all_expr_data$Escape)
all_expr_data$Infectivity <- as.numeric(all_expr_data$Infectivity)
all_expr_data <- na.omit(all_expr_data)

inche_cm <- 2.54
pdf(file.path(args$output, "All.Expr.pie.pdf"), width=20/inche_cm, height=18/inche_cm)
plot(c(0.7, 3.5), c(0, 6),type="n",xlab="Infectivity",ylab="Escape",axes=TRUE)
# color_df <- data.frame(
#   Time=c("2019-12", sprintf("2020-%02d", 1:12), sprintf("2021-%02d", 1:9)),
#   Color=rev(c("#FF0500","#FD3F00","#FDA201","#FFBC03","#E5D005","#A7EB10","#A7EB10","#8EEE24","#5AF44F","#05FF94","#01FDAF","#02FCC7","#02FAE1","#04EFF3","#09CDF4","#0DABF6","#1288F8","#185DFA","#1E32FC","#2501FF","#2500FF","#2500FF"))
#   )
# color_df$Color <- as.character(color_df$Color)

color_df <- data.frame(Time=c("2019-12", sprintf("2020-%02d", 1:12), sprintf("2021-%02d", 1:12), sprintf("2022-%02d", 1:7)))
color_func <- colorRamp2(seq(1, 26, 5), rev(c('#FF0000','#FF7F00','#FFFF00','#00FF00','#00FFFF', '#0000FF')))
color_df$Color <- sapply(1:nrow(color_df), color_func)

for(i in 1:nrow(all_expr_data)){
  this_mut_name <- all_expr_data$MutSetName[i]
  label_name <- all_expr_data$Label[i]
  tmp_mouth_df <- label2num_df[label2num_df$MutSetName == this_mut_name,]
  tmp_mouth_df <- left_join(tmp_mouth_df, color_df)
  x <- all_expr_data$Infectivity[i]
  y <- all_expr_data$Escape[i]
  r <- log10(all_expr_data$TotalIsolateNum[i]) / 100
  floating.pie(x,y,tmp_mouth_df$IsolateNumPerMonth,radius=r,col=tmp_mouth_df$Color)
}
tmp_mouth_df <- left_join(WT_date_num, color_df)
floating.pie(
  mut_expr$Infectivity[mut_expr$ExprName=="WT"],
  mut_expr$Escape[mut_expr$ExprName=="WT"],
  tmp_mouth_df$IsolateNumPerMonth,radius=log10(sum(WT_date_num$IsolateNumPerMonth)) / 200,col=tmp_mouth_df$Color)
text(all_expr_data$Infectivity, all_expr_data$Escape, all_expr_data$Label, cex=0.2)
text(mut_expr$Infectivity[mut_expr$ExprName=="WT"], mut_expr$Escape[mut_expr$ExprName=="WT"], "WT", cex=0.2)
floating.pie(0.8,5.4, c(1),radius=1/100,col=c("black", "white"))
floating.pie(0.9,5.4, c(1),radius=2/100,col=c("black", "white"))
floating.pie(1,5.4, c(1),radius=3/100,col=c("black", "white"))
floating.pie(1.1,5.4, c(1),radius=4/100,col=c("black", "white"))
floating.pie(1.2,5.4, c(1),radius=5/100,col=c("black", "white"))
floating.pie(1.3,5.4, c(1),radius=6/100,col=c("black", "white"))
dev.off()

p <- ggplot() +
  geom_point(data=all_expr_data, mapping=aes(x=Escape, y=Infectivity, color=as.numeric(Date), size=log10(TotalIsolateNum)))+
  geom_text(data=all_expr_data, mapping=aes(x=Escape, y=Infectivity, label=Label), size=1.2, color="black") +
  scale_color_gradientn(
    colours = rev(c('#FF0000','#FF7F00','#FFFF00','#00FF00','#00FFFF', '#0000FF')),
    breaks = as.numeric(as.Date(c("2020-01-01", "2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01", "2022-07-01"), format = "%Y-%m-%d")),
    labels = c("2020-01-01", "2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01", "2022-07-01"),
    limits=as.numeric(c(as.Date(c("2020-01-01", "2021-09-01"), format = "%Y-%m-%d")))
  ) +
  scale_size_continuous(limits = c(0, 6), range = c(0, 4)) +
  scale_x_continuous(limits = c(min(all_expr_data$Escape, na.rm = T), max(all_expr_data$Escape, na.rm = T))) +
  scale_y_continuous(limits = c(min(all_expr_data$Infectivity, na.rm = T), max(all_expr_data$Infectivity, na.rm = T))) +
  labs(size="Log10IsolateNum", color="First time") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.key.size = unit(2, "mm"),
        legend.title = element_text(family="ArialMT", size = 6),
        legend.text = element_text(family="ArialMT", size = 6),
        panel.grid = element_blank())
p
ggsave(file.path(args$output, "All.Expr.pdf"), p, height = 8, width = 10, units = "cm", colormodel = "cmyk")

f_time_points <- file.path(args$output, "timepoint")
if(!dir.exists(f_time_points)){
  dir.create(f_time_points)
}
mut_expr <- left_join(mut_expr, expr2mut_set_df)
mut_expr$Escape <- as.numeric(mut_expr$Escape)
mut_expr$Infectivity <- as.numeric(mut_expr$Infectivity)
for(time in time_li){
  tmp_mut_time_num <- mut_time_num[mut_time_num$Time==time,]
  tmp_mut_time_num <- tmp_mut_time_num[!duplicated(tmp_mut_time_num$MutSetName),c("MutSetName", "Time", "Year", "Month", "IsolateNum")]
  tmp_mut_time_num$IsolateRatio <- tmp_mut_time_num$IsolateNum / total_time_num$IsolatePerTime[sprintf("%s-%s", total_time_num$Year, total_time_num$Month)==time]
  tmp_mut_expr <- left_join(mut_expr, tmp_mut_time_num)
  tmp_mut_expr_found <- tmp_mut_expr[!is.na(tmp_mut_expr$Time),]
  tmp_mut_expr_not_found <- tmp_mut_expr[is.na(tmp_mut_expr$Time),]
  p <- ggplot() +
    geom_point(data=tmp_mut_expr_found, mapping=aes(x=Escape, y=Infectivity, size=IsolateRatio), color="black")+
    geom_point(data=tmp_mut_expr_not_found, mapping=aes(x=Escape, y=Infectivity), color="grey70", size=0.5, shape=4) +
    geom_text(data=tmp_mut_expr, mapping=aes(x=Escape, y=Infectivity, label=ExprName), size=1.2, color="black") +
    scale_size_continuous(range = c(0, 6), limits = c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), label=c("0%", "20%", "40%", "60%", "80%", "100%")) +
    labs(size="IsolateRatio", title = time) +
    theme_bw() +
    theme(text = element_text(family="ArialMT", size = 6),
          title = element_text(family="ArialMT", size = 6),
          axis.text = element_text(colour = "black"),
          legend.key.size = unit(2, "mm"),
          legend.title = element_text(family="ArialMT", size = 6),
          legend.text = element_text(family="ArialMT", size = 6),
          panel.grid = element_blank())
  ggsave(file.path(f_time_points, sprintf("%s.Expr.pdf", time)), p, height = 8, width = 10, units = "cm", colormodel = "cmyk")
}
