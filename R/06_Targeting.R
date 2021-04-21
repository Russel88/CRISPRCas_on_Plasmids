library(ggplot2)

load("Prepared.RData")

# Load data
pl_pl <- read.table("../Collect/PLSDB_plsdb_orf.m8")
hs_pl <- read.table("../Collect/Hosts_plsdb_orf.m8")
apl_pl <- read.table("../Collect/Archaea_plsdb_orf.m8")
ahs_pl <- read.table("../Collect/Hosts_plsdb_orf.m8")

pl_pl <- rbind(pl_pl, apl_pl)
hs_pl <- rbind(hs_pl, ahs_pl)

colnames(pl_pl) <- c("Spacer", "Target", "TStart", "TEnd", "ORF")
colnames(hs_pl) <- c("Spacer", "Target", "TStart", "TEnd", "ORF")

pl_pl$Source_Type <- "Plasmid"
hs_pl$Source_Type <- "Chromosome"

# Combine
all_spacers <- rbind(pl_pl[, c("Spacer", "Target", "Source_Type")], 
                     hs_pl[, c("Spacer", "Target", "Source_Type")])

# Get Focal names
all_spacers$CRISPR <- gsub("@.*", "", all_spacers$Spacer)
all_spacers$Source <- gsub("[._-][-,0-9]*$", "", all_spacers$CRISPR)
all_spacers$Target <- gsub("\\.[0-9]*$", "", all_spacers$Target)

all_spacers <- all_spacers[all_spacers$Source != all_spacers$Target, ]

# Derep
dereps <- gsub("\\.[0-9]*$", "", c(as.character(all_plasmids[all_plasmids$Derep, "Acc"]), 
                                   as.character(all_hosts[all_hosts$Derep, "Acc"])))
all_spacers <- all_spacers[all_spacers$Source %in% dereps & all_spacers$Target %in% dereps, ]

# Aggregate
all_spacers_unq <- aggregate(CRISPR ~ Spacer + Source_Type + Source + Target, data = all_spacers, function(x) 1)

# Columns for merging
all_spacers_unq$CRISPR <- gsub("@[0-9]*$", "", all_spacers_unq$Spacer)
all_plasmids$Source <- gsub("\\.[0-9]*$", "", all_plasmids$Acc)
all_hosts$Source <- gsub("\\.[0-9]*$", "", all_hosts$Acc)

# Merge Souce and Target info
all_spacers_unq <- merge(all_spacers_unq, rbind(all_plasmids[, c("Source", "Phylum")],
                                                all_hosts[, c("Source", "Phylum")]), by = "Source")
all_spacers_unq <- merge(all_spacers_unq, all_plasmids[, c("Source", "Mobility", "Phylum")], by.x = "Target", by.y = "Source")

# Only Proteos
all_spacers_unq <- all_spacers_unq[!is.na(all_spacers_unq$Phylum.x) & all_spacers_unq$Phylum.x == "Proteobacteria" &
                                       !is.na(all_spacers_unq$Phylum.y) & all_spacers_unq$Phylum.y == "Proteobacteria", ]


# Spacer by target matrix
all_spacers_unq_mat <- reshape2::dcast(Spacer ~ Mobility, data = all_spacers_unq, fun.aggregate = length, value.var = "Target")
rownames(all_spacers_unq_mat) <- all_spacers_unq_mat$Spacer
all_spacers_unq_mat <- all_spacers_unq_mat[, -1]

# Normalize to 1 for each spacer
all_spacers_unq_mat <- apply(all_spacers_unq_mat, 1, function(x) x/sum(x))

# Combine
spacer_mob <- data.frame(Spacer = c(colnames(all_spacers_unq_mat),colnames(all_spacers_unq_mat)),
                         N = c(all_spacers_unq_mat[1, ], colSums(all_spacers_unq_mat[2:3, ])),
                         Con = c(rep("Conjugative", ncol(all_spacers_unq_mat)),
                                 rep("Non-conjugative", ncol(all_spacers_unq_mat))))

spacer_info <- all_spacers_unq[, c("Spacer", "Source_Type")]
spacer_info <- spacer_info[!duplicated(spacer_info$Spacer), ]

spacer_mob <- merge(spacer_mob, spacer_info[, c("Spacer", "Source_Type")], by = "Spacer")

df_agg <- aggregate(N ~ Source_Type + Con, spacer_mob, sum)

# Add proteos in general
tabcon <- table(all_plasmids[all_plasmids$Derep & !is.na(all_plasmids$Phylum) &all_plasmids$Phylum == "Proteobacteria", "Mobility"])

df_agg <- rbind(df_agg, data.frame(Source_Type = c("All plasmids", "All plasmids"),
                                   Con = c("Conjugative", "Non-conjugative"),
                                   N = c(tabcon[1], sum(tabcon[2:3]))))

p <- ggplot(df_agg, aes(Source_Type, N, fill = Con)) +
    theme_bw() +
    geom_bar(stat = "identity", position = "fill") +
    coord_flip() +
    ylab("Percentage") +
    xlab("Spacer source") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual("Mobility of target plasmid", values = viridis::viridis(2)) +
    geom_vline(xintercept = 1.5)
p
ggsave(p, file = "Figures/Fig4_target_mobility.pdf", width = 16, height = 6, units = "cm")
