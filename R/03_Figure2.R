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
cc_median <- median(all_plasmids[all_plasmids$Group == "CRISPR-Cas", "Length"])
oc_median <- median(all_plasmids[all_plasmids$Group == "Orphan CRISPR or cas", "Length"])
nc_median <- median(all_plasmids[all_plasmids$Group == "No CRISPR or cas", "Length"])

library(mixtools)
m <- normalmixEM2comp(log10(all_plasmids[all_plasmids$Group == "No CRISPR or cas", "Length"]), sigsqrd = c(1,2), mu = c(3,5), lambda = 0.5)
param <- 10^m$mu

line_df <- data.frame(Group = c("CRISPR-Cas","Orphan CRISPR or cas","No CRISPR or cas","No CRISPR or cas"),
                      LengthS = paste(round(c(cc_median,oc_median,param[1],param[2])/1000),"kbp"),
                      Length = c(cc_median,oc_median,param[1],param[2]),
                      Place = c(2.1,1.4,1,1))

n <- table(all_plasmids$Group)
all_plasmids$Group <- factor(all_plasmids$Group, 
                                labels = paste0(levels(all_plasmids$Group), " (n=", n[levels(all_plasmids$Group)], ")"))
line_df$Group <- factor(line_df$Group, 
                             labels = paste0(levels(line_df$Group), " (n=", n[levels(line_df$Group)], ")"))

p <- ggplot(all_plasmids, aes(Length, fill = Group)) +
    theme_bw() +
    geom_density() +
    geom_vline(data = line_df, aes(xintercept = Length)) +
    geom_text(data = line_df, aes(label = LengthS, y = Place, x = Length), hjust = -0.2, size = 3) +
    scale_x_log10(labels = COEF::fancy_scientific) +
    xlab("Plasmid size (bp)") +
    ylab("Density") +
    facet_grid(Group ~ ., scales = "free") +
    theme(legend.position = "none")
p
ggsave(p, file = "Figures/Fig2_size2.pdf", width = 8, height = 16, units = "cm")
write.csv(all_plasmids[, c("Group", "Length")], file = "Tables/Fig2_size2.csv", quote = FALSE, row.names = FALSE)

# Mobility
proteo_plasmids$Mobility <- sub('^(\\w?)', '\\U\\1', proteo_plasmids$Mobility, perl=TRUE)
proteo_plasmids$Group <- factor(proteo_plasmids$Group, levels = c("No CRISPR or cas", "Orphan CRISPR or cas", "CRISPR-Cas"))
proteo_plasmids[!proteo_plasmids$Mobility == "Conjugative", "Mobility"] <- "Non-conjugative"

n <- table(proteo_plasmids$Group)
proteo_plasmids$Group <- factor(proteo_plasmids$Group, 
                                labels = paste0(levels(proteo_plasmids$Group), "\n(n=", n[levels(proteo_plasmids$Group)], ")"))

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

inc_prev <- data.frame(table(unlist(lapply(as.character(inc_plasmids$Inc), function(x) strsplit(x, ",")))))
inc_prev <- inc_prev[!grepl("rep", inc_prev$Var1), ]

combine_inc <- merge(combine_inc, inc_prev, by = "Var1")

combine_inc$IncN <- paste0(combine_inc$Var1, " (n=", combine_inc$Freq, ")")

p <- ggplot(combine_inc, aes(IncN, Ratio, ymin = Ratio_min, ymax = Ratio_max)) +
    theme_bw() +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    ylab("Difference (Observed - Random)") +
    xlab(NULL)
p
ggsave(p, file = "Figures/Fig2_inc_all_supp.pdf", width = 16, height = 16, units = "cm")
write.csv(inc_plasmids[, c("Inc", "CRISPRorCas")], file = "Tables/Fig2_Inc.csv", quote = FALSE, row.names = FALSE)

# Bar
all_inc <- data.frame(table(unlist(lapply(as.character(all_plasmids[all_plasmids$Derep, "Inc"]),
                                          function(x) strsplit(x, ",")))))

obs_inc$Group <- "CRISPR-Cas"
all_inc$Group <- "All"
merge_inc <- merge(obs_inc, all_inc, by = "Var1")
merge_inc$Perc <- merge_inc$Freq.x / merge_inc$Freq.y


merge_inc[is.na(merge_inc$Freq.x), "Freq.x"] <- 0
ncc_inc <- data.frame(Var1 = merge_inc$Var1, 
                      Freq = merge_inc$Freq.y - merge_inc$Freq.x,
                      Group = "No CRISPR-Cas")

merge_inc <- merge_inc[!grepl("rep", merge_inc$Var1), ]

merge_inc$Var1 <- paste0(merge_inc$Var1, " (n=", merge_inc$Freq.y, ")")

merge_inc <- merge_inc[merge_inc$Freq.y >= 50, ]

p <- ggplot(merge_inc, aes(Var1, Perc)) +
    theme_bw() +
    coord_flip() +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) +
    xlab(NULL) +
    ylab("Proportion with CRISPR or Cas loci")
p
ggsave(p, file = "Figures/Fig2_inc.pdf", width = 12, height = 9, units = "cm")
