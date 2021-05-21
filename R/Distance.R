library(ggplot2)

load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb_orf.m8")
apl_pl <- read.table("../Collect/Archaea_plsdb_orf.m8")

pl_pl <- rbind(pl_pl, apl_pl)

colnames(pl_pl) <- c("Spacer", "Target", "TStart", "TEnd", "ORF")
pl_pl <- pl_pl[, c("Spacer", "Target", "ORF")]

distance_df <- read.table("../Collect/bindash_sub.tab")
distance_df$V1 <- gsub("\\.[0-9]*$", "", distance_df$V1)
distance_df$V2 <- gsub("\\.[0-9]*$", "", distance_df$V2)
distance_df$Pair <- apply(distance_df, 1, function(x) paste(min(x[1:2]),max(x[1:2])))

# Get Focal names
pl_pl$CRISPR <- gsub("@.*", "", pl_pl$Spacer)
pl_pl$Source <- gsub("[._-][-,0-9]*$", "", pl_pl$CRISPR)
pl_pl$Target <- gsub("\\.[0-9]*$", "", pl_pl$Target)

# Derep
dereps <- gsub("\\.[0-9]*$", "", as.character(all_plasmids[all_plasmids$Derep, "Acc"]))
pl_pl <- pl_pl[pl_pl$Source %in% dereps & pl_pl$Target %in% dereps, ]

# Aggregate
pl_pl_agg <- aggregate(Spacer ~ Source + Target, data = pl_pl, function(x) length(unique(x)))
pl_pl_agg$Pair <- apply(pl_pl_agg, 1, function(x) paste(min(x[1:2]),max(x[1:2]))) 

# Merge
pl_pl_merged <- merge(pl_pl_agg, distance_df[, c("Pair", "V3")], by = "Pair", all.x = TRUE)
pl_pl_merged[is.na(pl_pl_merged$V3), "V3"] <- 1
pl_pl_merged <- pl_pl_merged[pl_pl_merged$Source != pl_pl_merged$Target, ]

# Possible pairs
all_derep <- gsub("\\.[0-9]*$","",all_plasmids[all_plasmids$Derep, "Acc"])
pairs <- sapply(gsub("\\.[0-9]*$","",cris_plasmid$Acc), function(x) sapply(all_derep, function(y) paste(min(x,y),max(x,y))))
pairsdf <- data.frame(Pair = as.vector(pairs))
pairsdf <- merge(pairsdf, distance_df, by = "Pair")
pairsdf[is.na(pairsdf$V3), "V3"] <- 1

rands <- sapply(1:1e5, function(x) mean(sample(pairsdf$V3, nrow(pl_pl_merged))))

p <- ggplot() +
    theme_bw() +
    geom_density(data = data.frame(Dist = rands), aes(Dist, ..scaled..), fill = "grey") +
    geom_vline(xintercept = mean(pl_pl_merged$V3), color = "darkred", size = 2) +
    xlim(c(0.98, 1)) +
    xlab("Plasmid-plasmid kmer distance") +
    ylab("Density")
p
ggsave(p, file = "Figures/Fig4_distance.pdf", width = 8, height = 7, units = "cm")

save.image("Distance.RData")
