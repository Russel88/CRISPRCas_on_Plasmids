
load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb_orf.m8")
colnames(pl_pl) <- c("Spacer", "Target", "TStart", "TEnd", "ORF")

# Get Focal names
pl_pl$CRISPR <- gsub("@.*", "", pl_pl$Spacer)
pl_pl$Source <- gsub("[._-][-,0-9]*$", "", pl_pl$CRISPR)
pl_pl$Target <- gsub("\\.[0-9]*$", "", pl_pl$Target)

# Remove mini-arrays
pl_pl <- pl_pl[pl_pl$CRISPR %in% cris_plasmid$CRISPR, ]

# Aggregate
pl_pl_agg <- aggregate(ORF ~ Source + Target, data = pl_pl, function(x) sum(!is.na(x)))
colnames(pl_pl_agg)[3] <- "Weight"

# Self-targeting pairs
selfs <- nrow(pl_pl_agg[pl_pl_agg$Source == pl_pl_agg$Target, ])
pl_pl_agg <- pl_pl_agg[pl_pl_agg$Source != pl_pl_agg$Target, ]

# Cross-targeting pairs
crosstarget <- pl_pl_agg[which(paste(pl_pl_agg$Source, pl_pl_agg$Target) %in% paste(pl_pl_agg$Target, pl_pl_agg$Source)), ]

# Get nodes
nodes <- unique(c(pl_pl_agg$Source, pl_pl_agg$Target))
assem <- sapply(nodes, function(x) cell[grepl(x, cell$Plasmids), "Cell"][1])
node_df <- data.frame(Id = nodes,
                      Cell = gsub(".[0-9]$","",assem))
node_df <- merge(node_df, tax, by = "Cell", all.x = TRUE)
node_df$CRISPR <- node_df$Id %in% pl_pl_agg$Source

# Derep
drep <- gsub(".[0-9]$", "", all_plasmids[all_plasmids$Derep, "Acc"])
pl_pl_agg <- pl_pl_agg[pl_pl_agg$Source %in% drep & pl_pl_agg$Target %in% drep,]
node_df <- node_df[node_df$Id %in% drep, ]
node_df <- node_df[node_df$Id %in% c(pl_pl_agg$Source, pl_pl_agg$Target), ]

write.table(node_df, file = "Figures/Network_nodes.csv", row.names = FALSE, quote = FALSE, sep = ",")
write.table(pl_pl_agg, file = "Figures/Network_edges.csv", row.names = FALSE, quote = FALSE, sep = ",")
