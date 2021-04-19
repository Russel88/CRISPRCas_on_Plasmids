library(ggplot2)

load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb.m8")
apl_pl <- read.table("../Collect/Archaea_plsdb.m8")

pl_pl <- rbind(pl_pl, apl_pl)

colnames(pl_pl) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")

distance_df <- read.table("../Collect/bindash_sub.tab")
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

save.image("Distance.RData")
