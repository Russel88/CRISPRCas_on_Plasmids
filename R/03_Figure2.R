library(ggplot2)

# Load data
load("Prepared.RData")

# Derep
all_plasmids <- all_plasmids[all_plasmids$Derep, ]

# Split into groups
all_plasmids$Group <- ifelse(all_plasmids$CRISPRCas > 0, "CRISPR-Cas", 
                             ifelse(all_plasmids$Cas_Orphan > 0 | all_plasmids$CRISPR_Orphan > 0, "Orphan CRISPR or cas", "No CRISPR or cas"))
all_plasmids$Group <- factor(all_plasmids$Group, levels = c("CRISPR-Cas", "Orphan CRISPR or cas", "No CRISPR or cas"))
proteo_plasmids <- all_plasmids[!is.na(all_plasmids$Phylum) & all_plasmids$Phylum == "Proteobacteria", ]
inc_plasmids <- all_plasmids[grepl("Inc|Col", all_plasmids$Inc), ]

# Size
size_quants <- quantile(all_plasmids$Length)
all_plasmids$SizeQ <- as.factor(ifelse(all_plasmids$Length < size_quants[2], paste0("[",round(size_quants[1]/1000),"-",round(size_quants[2]/1000),"kbp)"),
                             ifelse(all_plasmids$Length < size_quants[3],  paste0("[",round(size_quants[2]/1000),"-",round(size_quants[3]/1000),"kbp)"),
                                    ifelse(all_plasmids$Length < size_quants[4],  paste0("[",round(size_quants[3]/1000),"-",round(size_quants[4]/1000),"kbp)"), 
                                           paste0("[",round(size_quants[4]/1000),"-",round(size_quants[5]/1000),"kbp]")))))

all_plasmids_a <- aggregate(Length ~ SizeQ + Group, all_plasmids, function(x) sum(x>0))
all_plasmids_a$Percent <- all_plasmids_a$Length / nrow(all_plasmids) * 4

all_plasmids_a$SizeQ <- factor(all_plasmids_a$SizeQ, levels = levels(all_plasmids_a$SizeQ)[c(3,2,4,1)])

p <- ggplot(all_plasmids_a[all_plasmids_a$Group != "No CRISPR or cas", ], aes(SizeQ, Percent, fill = Group)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab("Size quantile")
p
ggsave(p, file = "Figures/Fig2_size.pdf", width = 12, height = 7, units = "cm")
write.csv(all_plasmids_a[, c("SizeQ", "Group", "Percent")], file = "Tables/Fig2_size.csv", quote = FALSE, row.names = FALSE)

p <- ggplot(all_plasmids, aes(Length, fill = Group)) +
    theme_bw() +
    geom_density() +
    scale_x_log10(labels = COEF::fancy_scientific) +
    xlab("Plasmid size (bp)") +
    ylab("Density") +
    facet_grid(Group ~ .) +
    theme(legend.position = "none")
p
ggsave(p, file = "Figures/Fig2_size2.pdf", width = 6, height = 14, units = "cm")
write.csv(all_plasmids[, c("Group", "Length")], file = "Tables/Fig2_size2.csv", quote = FALSE, row.names = FALSE)

# Mobility
proteo_plasmids$Mobility <- sub('^(\\w?)', '\\U\\1', proteo_plasmids$Mobility, perl=TRUE)
levels(proteo_plasmids$Group) <- rev(levels(proteo_plasmids$Group))

p <- ggplot(proteo_plasmids, aes(Group, fill = Mobility)) +
    theme_bw() +
    geom_bar(position = "fill") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = viridis::viridis(3)) +
    ylab(NULL) +
    xlab(NULL)
p
ggsave(p, file = "Figures/Fig2_mobility.pdf", width = 16, height = 5, units = "cm")
write.csv(proteo_plasmids[, c("Group", "Mobility")], file = "Tables/Fig2_mobility.csv", quote = FALSE, row.names = FALSE)

# Inc
inc_plasmids$CRISPRorCas <- inc_plasmids$CRISPRs > 0 | inc_plasmids$Cas > 0
inc_plasmids$Inc <- gsub("ColRNAI_rep_cluster_[0-9]*", "ColRNAI", inc_plasmids$Inc)

randomize <- function(x){
    df <- inc_plasmids
    df$CRISPRorCas <- sample(df$CRISPRorCas)
    res <- data.frame(table(unlist(lapply(as.character(df[df$CRISPRorCas, "Inc"]), function(x) strsplit(x, ",")))))
    return(res)
}

obs_inc <- data.frame(table(unlist(lapply(as.character(inc_plasmids[inc_plasmids$CRISPRorCas, "Inc"]),
                                         function(x) strsplit(x, ",")))))

set.seed(42)
rand_inc <- lapply(1:1000, randomize)

test <- Reduce(function(x, y) merge(x, y, by = "Var1", all = TRUE), rand_inc)
test[is.na(test)] <- 0

random_inc <- data.frame(Inc = test$Var1,
                         Freq = rowMeans(test[, -1]),
                         Freq.sd = apply(test[, -1], 1, sd))

combine_inc <- merge(obs_inc, random_inc, by.x = "Var1", by.y = "Inc", all = TRUE)
combine_inc <- combine_inc[!grepl("rep", combine_inc$Var1), ]
combine_inc[is.na(combine_inc$Freq.x), "Freq.x"] <- 0

combine_inc$Ratio <- combine_inc$Freq.x - combine_inc$Freq.y
combine_inc$Rand_max <- (combine_inc$Freq.y + combine_inc$Freq.sd)
combine_inc$Rand_min <- (combine_inc$Freq.y - combine_inc$Freq.sd)
combine_inc[combine_inc$Rand_min < 0, "Rand_min"] <- 0

combine_inc$Ratio_max <- combine_inc$Freq.x - combine_inc$Rand_min
combine_inc$Ratio_min <- combine_inc$Freq.x - combine_inc$Rand_max

combine_inc$Var1 <- as.character(combine_inc$Var1)

p <- ggplot(combine_inc, aes(Var1, Ratio, ymin = Ratio_min, ymax = Ratio_max)) +
    theme_bw() +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    ylab("Difference (Observed - Random)")
p
ggsave(p, file = "Figures/Fig2_inc_all.pdf", width = 16, height = 12, units = "cm")
write.csv(inc_plasmids[, c("Inc", "CRISPRorCas")], file = "Tables/Fig2_Inc.csv", quote = FALSE, row.names = FALSE)

# Only prevalent inc
inc_prev <- data.frame(table(unlist(lapply(as.character(inc_plasmids$Inc), function(x) strsplit(x, ",")))))
inc_prev <- inc_prev[!grepl("rep", inc_prev$Var1), ]

combine_inc <- merge(combine_inc, inc_prev, by = "Var1")

p <- ggplot(combine_inc[combine_inc$Freq >= 50, ],
            aes(Var1, Ratio/Freq, ymin = Ratio_min/Freq, ymax = Ratio_max/Freq)) +
    theme_bw() +
    geom_hline(yintercept = 0) +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    ylab("Difference (Observed - Random) %") +
    xlab(NULL) +
    scale_y_continuous(labels = scales::percent)
p
ggsave(p, file = "Figures/Fig2_inc.pdf", width = 12, height = 9, units = "cm")
