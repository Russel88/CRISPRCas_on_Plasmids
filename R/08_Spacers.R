
load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb.m8")
pl_ph <- read.table("../Collect/PLSDB_IMG.m8")
hs_pl <- read.table("../Collect/Hosts_plsdb.m8")
hs_ph <- read.table("../Collect/Hosts_IMG.m8")
apl_pl <- read.table("../Collect/Archaea_plsdb.m8")
apl_ph <- read.table("../Collect/PLSDB_IMG.m8")
ahs_pl <- read.table("../Collect/Hosts_plsdb.m8")
ahs_ph <- read.table("../Collect/Hosts_IMG.m8")

pl_pl <- rbind(pl_pl, apl_pl)
pl_ph <- rbind(pl_ph, apl_ph)
hs_pl <- rbind(hs_pl, ahs_pl)
hs_ph <- rbind(hs_ph, ahs_ph)

colnames(pl_pl) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")
colnames(pl_ph) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")
colnames(hs_pl) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")
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
all_spacers <- rbind(pl_pl, pl_ph, hs_pl, hs_ph)

# Get Focal names
all_spacers$CRISPR <- gsub("@.*", "", all_spacers$Spacer)
all_spacers$Source <- gsub("[._-][-,0-9]*$", "", all_spacers$CRISPR)
all_spacers$Target <- gsub("\\.[0-9]*$", "", all_spacers$Target)

all_spacers[all_spacers$Source == all_spacers$Target, "Target_Type"] <- "Self"

# Derep
dereps <- gsub("\\.[0-9]*$", "", c(as.character(all_plasmids[all_plasmids$Derep, "Acc"]), as.character(all_hosts[all_hosts$Derep, "Acc"])))
all_spacers <- all_spacers[all_spacers$Source %in% dereps &
                               (all_spacers$Target %in% dereps | all_spacers$Target_Type == "Phage"), ]

# Aggregate
all_spacers_agg <- aggregate(Spacer ~ Source + Target + Source_Type + Target_Type, data = all_spacers, function(x) length(unique(x)))
colnames(all_spacers_agg)[ncol(all_spacers_agg)] <- "Spacers"

# Write plasmid target count table
plas_count <- all_spacers_agg[all_spacers_agg$Source_Type == "Plasmid", ]
all_plasmids$Acc <- gsub("\\.[0-9]*$", "", all_plasmids$Acc)
plas_count <- merge(plas_count, all_plasmids[, c("Acc", "Cell")], by.x = "Source", by.y = "Acc", all.x = TRUE)
plas_count <- merge(plas_count, all_plasmids[, c("Acc", "Cell")], by.x = "Target", by.y = "Acc", all.x = TRUE)
plas_count[!is.na(plas_count$Cell.x) & !is.na(plas_count$Cell.y) & 
               plas_count$Cell.x == plas_count$Cell.y, "Target_Type"] <- "Plasmid_in_Cell"
plas_count[plas_count$Source == plas_count$Target, "Target_Type"] <- "Self"

plas_count <- aggregate(Target ~ Source + Target_Type, data = plas_count, function(x) length(unique(x)))
pl_many_cast <- reshape2::dcast(plas_count, Source~Target_Type, value.var = "Target", fill = 0)

write.table(pl_many_cast, file = "Tables/Plasmid_target_count.tab", row.names = FALSE, quote = FALSE, sep = "\t")
