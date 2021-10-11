library(eulerr)
library(ggplot2)

load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb_orf.m8")
pl_ph <- read.table("../Collect/PLSDB_IMG.m8")
hs_pl <- read.table("../Collect/Hosts_plsdb_orf.m8")
hs_ph <- read.table("../Collect/Hosts_IMG.m8")
apl_pl <- read.table("../Collect/Archaea_plsdb_orf.m8")
apl_ph <- read.table("../Collect/PLSDB_IMG.m8")
ahs_pl <- read.table("../Collect/Hosts_plsdb_orf.m8")
ahs_ph <- read.table("../Collect/Hosts_IMG.m8")

pl_pl <- rbind(pl_pl, apl_pl)
pl_ph <- rbind(pl_ph, apl_ph)
hs_pl <- rbind(hs_pl, ahs_pl)
hs_ph <- rbind(hs_ph, ahs_ph)

colnames(pl_pl) <- c("Spacer", "Target", "TStart", "TEnd", "ORF")
colnames(pl_ph) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")
colnames(hs_pl) <- c("Spacer", "Target", "TStart", "TEnd", "ORF")
colnames(hs_ph) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")

pl_pl$Source_Type <- "Plasmid"
pl_ph$Source_Type <- "Plasmid"
hs_pl$Source_Type <- "Chromosome"
hs_ph$Source_Type <- "Chromosome"

pl_pl$Target_Type <- "Plasmid"
pl_ph$Target_Type <- "Phage"
hs_pl$Target_Type <- "Plasmid"
hs_ph$Target_Type <- "Phage"

# Combine
all_spacers <- rbind(pl_pl[, c("Spacer", "Target", "Source_Type", "Target_Type")], 
                     pl_ph[, c("Spacer", "Target", "Source_Type", "Target_Type")], 
                     hs_pl[, c("Spacer", "Target", "Source_Type", "Target_Type")], 
                     hs_ph[, c("Spacer", "Target", "Source_Type", "Target_Type")])

# Get Focal names
all_spacers$CRISPR <- gsub("@.*", "", all_spacers$Spacer)
all_spacers$Source <- gsub("[._-][-,0-9]*$", "", all_spacers$CRISPR)
all_spacers$Target <- gsub("\\.[0-9]*$", "", all_spacers$Target)

all_spacers <- all_spacers[all_spacers$Source != all_spacers$Target, ]

# Derep
dereps <- gsub("\\.[0-9]*$", "", c(as.character(all_plasmids[all_plasmids$Derep, "Acc"]), as.character(all_hosts[all_hosts$Derep, "Acc"])))
all_spacers <- all_spacers[all_spacers$Source %in% dereps &
                               (all_spacers$Target %in% dereps | all_spacers$Target_Type == "Phage"), ]

# Remove those from untyped orphans
orphans <- c(as.character(cris_plasmid[is.na(cris_plasmid$Prediction) & is.na(cris_plasmid$Subtype), "CRISPR"]),
             as.character(cris_host[is.na(cris_host$Prediction) & is.na(cris_host$Subtype), "CRISPR"]))

all_spacers <- all_spacers[!gsub("@[0-9]*$", "", all_spacers$Spacer) %in% orphans, ]

# Remove mini-arrays
all_spacers <- all_spacers[all_spacers$CRISPR %in% c(as.character(cris_plasmid$CRISPR), as.character(cris_host$CRISPR)), ]

# Aggregate
all_spacers_agg <- aggregate(Spacer ~ Source + Target + Source_Type + Target_Type, data = all_spacers, function(x) length(unique(x)))

pl_spacers_agg <- all_spacers_agg[all_spacers_agg$Source_Type == "Plasmid" & all_spacers_agg$Target_Type == "Plasmid", ]
all_plasmids$Acc2 <- gsub(".[0-9]*$", "", all_plasmids$Acc)

pl_spacers_agg <- merge(pl_spacers_agg, all_plasmids[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Acc2")], by.x = "Source", by.y = "Acc2")
pl_spacers_agg <- merge(pl_spacers_agg, all_plasmids[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Acc2")], by.x = "Target", by.y = "Acc2")

pl_spacers_agg <- pl_spacers_agg[!is.na(pl_spacers_agg$Species.x) & !is.na(pl_spacers_agg$Species.y), ]

# Matches to different levels
species <- nrow(pl_spacers_agg[pl_spacers_agg$Species.x == pl_spacers_agg$Species.y, ])
genus <- nrow(pl_spacers_agg[pl_spacers_agg$Genus.x == pl_spacers_agg$Genus.y, ]) - species
family <- nrow(pl_spacers_agg[pl_spacers_agg$Family.x == pl_spacers_agg$Family.y, ]) - species - genus
order <- nrow(pl_spacers_agg[pl_spacers_agg$Order.x == pl_spacers_agg$Order.y, ]) - species - genus - family
class <- nrow(pl_spacers_agg[pl_spacers_agg$Class.x == pl_spacers_agg$Class.y, ]) - species - genus - family - order
phylum <- nrow(pl_spacers_agg[pl_spacers_agg$Phylum.x == pl_spacers_agg$Phylum.y, ]) - species  - genus - family - order - class
kingdom <- nrow(pl_spacers_agg[pl_spacers_agg$Kingdom.x == pl_spacers_agg$Kingdom.y, ]) - species  - genus - family - order - class - phylum

set.seed(42)
species_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Species.x) == pl_spacers_agg$Species.y, ])))
genus_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Genus.x) == pl_spacers_agg$Genus.y, ]))) - species_rand
family_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Family.x) == pl_spacers_agg$Family.y, ]))) - species_rand - genus_rand
order_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Order.x) == pl_spacers_agg$Order.y, ]))) - species_rand - genus_rand - family_rand
class_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Class.x) == pl_spacers_agg$Class.y, ]))) - species_rand - genus_rand - family_rand - order_rand
phylum_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Phylum.x) == pl_spacers_agg$Phylum.y, ]))) - species_rand - genus_rand - family_rand - order_rand - class_rand
kingdom_rand <- median(sapply(1:100, function(x) nrow(pl_spacers_agg[sample(pl_spacers_agg$Kingdom.x) == pl_spacers_agg$Kingdom.y, ]))) - species_rand - genus_rand - family_rand - order_rand - class_rand - phylum_rand

# Combine and plot
df <- data.frame(Value = c(species, genus, family, order, class, phylum, kingdom,
                   species_rand, genus_rand, family_rand, order_rand, class_rand, phylum_rand, kingdom_rand),
                 Group = c(rep("Observed", 7), rep("Random", 7)),
                 Rank = factor(rep(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"), 2), 
                               levels = rev(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"))))
df$Prop <- df$Value / nrow(pl_spacers_agg)

sum(df[df$Group == "Observed", "Prop"])

p <- ggplot(df, aes(Group, Prop, fill = Rank)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(name = "Matches confined to:", values = viridis::viridis(7)) +
    xlab(NULL) +
    ylab("Plasmid-plasmid spacer-target proportion")
p    
ggsave(p, file = "Figures/Fig4_tax.pdf", units = "cm", width = 10, height = 8)
