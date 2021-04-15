
load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb.m8")
colnames(pl_pl) <- c("Spacer", "Target", "Identity", "Alignment", "MM", "Gap", "QStart", "QEnd", "TStart", "TEnd", "Eval", "Bit")

# Get Focal names
pl_pl$CRISPR <- gsub("@.*", "", pl_pl$Spacer)
pl_pl$Source <- gsub("[._-][-,0-9]*$", "", pl_pl$CRISPR)

pl_pl$Target <- gsub("\\.[0-9]*$", "", pl_pl$Target)

# Aggregate
pl_pl_agg <- aggregate(Eval ~ Source + Target, data = pl_pl, function(x) sum(x < 1))
colnames(pl_pl_agg)[3] <- "Weight"

# Self-targeting pairs
nrow(pl_pl_agg[pl_pl_agg$Source == pl_pl_agg$Target, ])

# Cross-targeting pairs
(sum(paste(pl_pl_agg$Source, pl_pl_agg$Target) %in% paste(pl_pl_agg$Target, pl_pl_agg$Source)) - .Last.value) / 2
crosstarget <- pl_pl_agg[which(paste(pl_pl_agg$Source, pl_pl_agg$Target) %in% paste(pl_pl_agg$Target, pl_pl_agg$Source)), ]
write.table(crosstarget, file = "Tables/crosstarget.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# Get nodes
nodes <- unique(c(pl_pl_agg$Source, pl_pl_agg$Target))
assem <- sapply(nodes, function(x) cell[grepl(x, cell$Plasmids), "Cell"][1])
node_df <- data.frame(Id = nodes,
                      Cell = gsub(".[0-9]$","",assem))
node_df <- merge(node_df, tax, by = "Cell", all.x = TRUE)
node_df$CRISPR <- node_df$Id %in% pl_pl_agg$Source

# Derep
drep_plsdb$Acc <- gsub(".[0-9]$", "", drep_plsdb$Acc)
pl_pl_agg <- pl_pl_agg[pl_pl_agg$Source %in% drep_plsdb$Acc & pl_pl_agg$Target %in% drep_plsdb$Acc,]
node_df <- node_df[node_df$Id %in% drep_plsdb$Acc, ]
node_df <- node_df[node_df$Id %in% c(pl_pl_agg$Source, pl_pl_agg$Target), ]

write.table(node_df, file = "Network_nodes.csv", row.names = FALSE, quote = FALSE, sep = ",")
write.table(pl_pl_agg, file = "Network_edges.csv", row.names = FALSE, quote = FALSE, sep = ",")
