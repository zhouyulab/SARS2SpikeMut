library(readxl)
library(argparse)
library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)

parser <- ArgumentParser()
parser$add_argument("--seq-info", required=TRUE, type="character", dest = "seq_info", metavar="seq_info.tsv")
parser$add_argument("--seq2uniq", required=TRUE, type="character", dest = "seq2uniq", metavar="seq2uniq.tsv")
parser$add_argument("--uniq-prot-label", required=TRUE, type="character", dest = "uniq_prot_label", metavar="uniq_prot_label.tsv")
parser$add_argument("--mut-expr", required=TRUE, type="character", dest = "mut_expr", metavar="mut_expr.xlsx")
parser$add_argument("-o", "--output", required=TRUE, type="character", dest = "output", metavar="output")

args <- commandArgs(TRUE)
args <- parser$parse_args(args)


mut_expr_df <- read_excel(args$mut_expr)
mut_expr_df <- mut_expr_df[,c("Mutant", "Escape", "Infectivity")]
names(mut_expr_df)[1] <- "MutType"

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

location_num_df <- seq2uniq_info_df %>% group_by(Location) %>% summarise(IsolateNum=n(), UniqProtNum=length(unique(UniqProtID)))
location_num_df <- location_num_df[order(location_num_df$IsolateNum, decreasing = T),]
write_tsv(location_num_df, file.path(args$output, "LocationNum.SourceData.tsv"))
location_num_df$Location <- factor(location_num_df$Location, levels = location_num_df$Location)
melt_location_num_df <- melt(location_num_df, "Location", variable.name = "Stat", value.name = "Number")
melt_location_num_df$Stat <- factor(melt_location_num_df$Stat, levels = c("UniqProtNum", "IsolateNum"), labels = c("#Unique proteins", "#Total isolates"))
p <- ggplot(melt_location_num_df, aes(x=Location, y=Number)) +
  geom_bar(stat="identity", fill="black") +
  facet_grid(Stat~., scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
ggsave(file.path(args$output, "LocationNum.pdf"), p, height = 6, width = 50, units = "cm", colormodel = "cmyk")

time_num_df <- seq2uniq_info_df %>% group_by(Year, Month) %>% summarise(IsolateNum=n(), UniqProtNum=length(unique(UniqProtID)))
time_num_df <- time_num_df[order(time_num_df$IsolateNum, decreasing = T),]
write_tsv(time_num_df, file.path(args$output, "TimeNum.SourceData.tsv"))
time_num_df$Label <- sprintf("%s-%s", time_num_df$Year, time_num_df$Month)
time_num_df$Label <- factor(time_num_df$Label, levels = time_num_df$Label)
melt_time_num_df <- melt(time_num_df, c("Label", "Year", "Month"), variable.name = "Stat", value.name = "Number")
melt_time_num_df$Stat <- factor(melt_time_num_df$Stat, levels = c("UniqProtNum", "IsolateNum"), labels = c("#Unique proteins", "#Total isolates"))
p <- ggplot(melt_time_num_df, aes(x=Label, y=Number)) +
  geom_bar(stat="identity", fill="black") +
  labs(x="Time", y="Number") +
  facet_grid(Stat~., scales = "free_y") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
ggsave(file.path(args$output, "TimeNum.pdf"), p, height = 6, width = 8, units = "cm", colormodel = "cmyk")



iso_num_per_seq_df <- seq2uniq_info_df %>% group_by(UniqProtID) %>% summarise(IsolateNum=n()) %>% group_by(IsolateNum) %>% summarise(UniqProtNum=n())
write_tsv(time_num_df, file.path(args$output, "IsolatePerUniqProt.SourceData.tsv"))
iso_num_per_seq_df$TotalIsoNum <- iso_num_per_seq_df$IsolateNum * iso_num_per_seq_df$UniqProtNum
iso_num_per_seq_level_df <- data.frame(
  IsolateNum=c("1", "2", "3", "4", "5-10", "11-20", "21-50", ">50"),
  UniqProtNum=c(
    iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum==1],
    iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum==2],
    iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum==3],
    iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum==4],
    sum(iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum>4 & iso_num_per_seq_df$IsolateNum<=10]),
    sum(iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum>10 & iso_num_per_seq_df$IsolateNum<=20]),
    sum(iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum>20 & iso_num_per_seq_df$IsolateNum<=50]),
    sum(iso_num_per_seq_df$UniqProtNum[iso_num_per_seq_df$IsolateNum>50])
    ),
  TotalIsolateNum=c(
    iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum==1],
    iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum==2],
    iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum==3],
    iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum==4],
    sum(iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum>4 & iso_num_per_seq_df$IsolateNum<=10]),
    sum(iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum>10 & iso_num_per_seq_df$IsolateNum<=20]),
    sum(iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum>20 & iso_num_per_seq_df$IsolateNum<=50]),
    sum(iso_num_per_seq_df$TotalIsoNum[iso_num_per_seq_df$IsolateNum>50])
  )
)
iso_num_per_seq_level_df$IsolateNum <- factor(iso_num_per_seq_level_df$IsolateNum, levels = iso_num_per_seq_level_df$IsolateNum)
p <- ggplot(iso_num_per_seq_level_df, aes(x=IsolateNum, y=UniqProtNum)) +
  geom_bar(stat="identity", fill="black") +
  labs(x="#Isolate", y="#Unique protein") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank())
ggsave(file.path(args$output, "IsolatePerUniqProt.pdf"), p, height = 5, width = 5, units = "cm", colormodel = "cmyk")

p <- ggplot(iso_num_per_seq_level_df, aes(x=IsolateNum, y=TotalIsolateNum)) +
  geom_bar(stat="identity", fill="black") +
  labs(x="#Isolate", y="#Total Isolate") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank())
ggsave(file.path(args$output, "TotalIsolatePerUniqProt.pdf"), p, height = 5, width = 5, units = "cm", colormodel = "cmyk")

p <- ggplot(iso_num_per_seq_df, aes(x=log10(IsolateNum), y=log10(UniqProtNum), size=log10(TotalIsoNum), alpha=log10(TotalIsoNum))) +
  geom_point() +
  scale_size_continuous(range = c(0, 2), breaks = c(0, 1, 2, 3, 4, 5, 6), labels = c("1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6")) +
  scale_alpha_continuous(range = c(0, 0.5), breaks = c(0, 1, 2, 3, 4, 5, 6), labels = c("1E0", "1E1", "1E2", "1E3", "1E4", "1E5", "1E6")) +
  scale_x_continuous(breaks = c(0, 2, 4, 6), labels = c("1E0", "1E2", "1E4", "1E6")) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), labels = c("1E0", "1E2", "1E4", "1E6")) +
  labs(x="#Isolate", y="#Unique protein", size="#Total isolate", alpha="#Total isolate") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        legend.title =  element_text(family="ArialMT", size = 6),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(0, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
ggsave(file.path(args$output, "IsolatePerUniqProt.point.pdf"), p, height = 5, width = 7, units = "cm", colormodel = "cmyk")

prot_label_df <- read_delim(args$uniq_prot_label, "\t", escape_double = FALSE, trim_ws = TRUE)



seq2label_df <- left_join(seq2uniq_info_df, prot_label_df)
write_tsv(seq2label_df, file.path(args$output, "seq2mut.SourceData.tsv"))

mut_label_per_uniq_seq_df <- seq2label_df %>% group_by(UniqProtID) %>% summarise(MutTypeNumPerUniqProt=sum(!is.na(unique(MutName))))
mut_label_per_uniq_seq_df <- left_join(mut_label_per_uniq_seq_df, seq2uniq_df %>% group_by(UniqProtID) %>% summarise(IsolateNum=n()))
mut_label_per_uniq_seq_df <- mut_label_per_uniq_seq_df %>% group_by(MutTypeNumPerUniqProt) %>% summarise(UniqProtNum=n(), TotalIsoNum=sum(IsolateNum))
write_tsv(mut_label_per_uniq_seq_df, file.path(args$output, "MutTypeNum.SourceData.tsv"))
melt_mut_label_per_uniq_seq_df <- melt(mut_label_per_uniq_seq_df, "MutTypeNumPerUniqProt", variable.name = "Stat", value.name = "Number")
melt_mut_label_per_uniq_seq_df$MutTypeNumPerUniqProt <- factor(melt_mut_label_per_uniq_seq_df$MutTypeNumPerUniqProt, levels = mut_label_per_uniq_seq_df$MutTypeNumPerUniqProt)
melt_mut_label_per_uniq_seq_df$Stat <- factor(melt_mut_label_per_uniq_seq_df$Stat, levels = c("UniqProtNum", "TotalIsoNum"), labels = c("#Unique proteins", "#Total isolates"))
p <- ggplot(melt_mut_label_per_uniq_seq_df, aes(x=MutTypeNumPerUniqProt, y=Number)) +
  geom_bar(stat="identity", fill="black") +
  facet_grid(Stat~., scales = "free_y") +
  labs(y="Number", x="#Mutation labels") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
ggsave(file.path(args$output, "MutTypeNum.pdf"), p, height = 6, width = 6, units = "cm", colormodel = "cmyk")

iso_num_in_loc_time <- seq2uniq_info_df %>% group_by(Location, Year, Month) %>% summarise(IsolateNumInLocTime=n())
mutation_loc_time <- seq2label_df %>% group_by(MutName, Location, Year, Month) %>% summarise(IsolateNum=n())
mutation_loc_time <- left_join(mutation_loc_time, iso_num_in_loc_time)
mutation_loc_time$IsolateRatio <- mutation_loc_time$IsolateNum / mutation_loc_time$IsolateNumInLocTime
write_tsv(mutation_loc_time, "MutationLocationTime.SourceData.tsv")

global_mutation_time <- seq2label_df %>% group_by(MutName, Year, Month) %>% summarise(IsolateNum=n())
global_mutation_time$Time <- sprintf("%s-%s", global_mutation_time$Year, global_mutation_time$Month)
global_mutation_time <- global_mutation_time[!is.na(global_mutation_time$MutName),]
iso_num_in_time <- seq2uniq_info_df %>% group_by(Year, Month) %>% summarise(IsolateNumInTime=n())
global_mutation_time <- left_join(global_mutation_time, iso_num_in_time)
global_mutation_time$IsolateRatio <- global_mutation_time$IsolateNum / global_mutation_time$IsolateNumInTime
global_mutation_time <- global_mutation_time[!is.na(global_mutation_time$Year),]
melt_global_mutation_time <- melt(global_mutation_time[,c("Time", "MutName", "IsolateNum", "IsolateRatio")], c("Time", "MutName"), variable.name = "Stat", value.name = "Number")
p <-  ggplot(global_mutation_time, aes(x=Time, y=IsolateNum)) +
  geom_bar(stat="identity", fill="black") +
  facet_wrap(~MutName) +
  labs(y="#Isolates", x="Time") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
ggsave(file.path(args$output, "MutTypeTimeNum.pdf"), p, height = 12, width = 25, units = "cm", colormodel = "cmyk")

p <-  ggplot(global_mutation_time, aes(x=Time, y=IsolateRatio)) +
  geom_bar(stat="identity", fill="black") +
  facet_wrap(~MutName) +
  labs(y="Isolate ratio", x="Time") +
  theme_bw() +
  theme(text = element_text(family="ArialMT", size = 6),
        title = element_text(family="ArialMT", size = 6),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank())
ggsave(file.path(args$output, "MutTypeTimeRatio.pdf"), p, height = 12, width = 25, units = "cm", colormodel = "cmyk")
