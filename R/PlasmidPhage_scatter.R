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

all_spacers[all_spacers$Source == all_spacers$Target, "Target_Type"] <- "Self"

# Derep
dereps <- gsub("\\.[0-9]*$", "", c(as.character(all_plasmids[all_plasmids$Derep, "Acc"]), as.character(all_hosts[all_hosts$Derep, "Acc"])))
all_spacers <- all_spacers[all_spacers$Source %in% dereps &
                               (all_spacers$Target %in% dereps | all_spacers$Target_Type == "Phage"), ]

# Aggregate
all_spacers_agg <- aggregate(Spacer ~ Source + Target + Source_Type + Target_Type, data = all_spacers, function(x) length(unique(x)))

# Target count
match_count <- aggregate(Spacer ~ Target_Type + Source + Source_Type, all_spacers_agg, function(x) length(unique(x)))

# Total spacers
total_spacer <- aggregate(Repeats+1 ~ Acc, rbind(cris_plasmid[, c("Acc", "Repeats")],
                                                    cris_host[, c("Acc", "Repeats")]), sum)
total_spacer$Source <- gsub("\\.[0-9]*$", "", total_spacer$Acc)
colnames(total_spacer)[2] <- "Spacers"

# Should not be counts, but percentage of all spacers matching phage or plasmid
match_count <- merge(match_count, total_spacer, by = "Source")
match_count$Prop <- match_count$Spacer / match_count$Spacers

plasmid_count <- match_count[match_count$Target_Type == "Plasmid", c("Source", "Source_Type", "Prop")]
phage_count <- match_count[match_count$Target_Type == "Phage", c("Source", "Source_Type", "Prop")]

plot_count <- merge(phage_count, plasmid_count, by = "Source")

p <- ggplot(plot_count, aes(Prop.x, Prop.y, color = Source_Type.x)) +
    theme_bw() +
    geom_density2d() +
    ylab("Proportion of spacers targeting plasmids") +
    xlab("Proportion of spacers targeting phages") +
    geom_smooth(method = "lm", formula = y~x, size = 2) +
    scale_x_log10(labels = scales::percent) +
    scale_y_log10(labels = scales::percent) +
    scale_color_manual(name = "Source", values = c("blue2", "red2"))+
    theme(legend.position = c(0.18,0.85))
p
ggsave(p, file = "Figures/Fig3_density.pdf", width = 11, height = 9, units = "cm")

p <- ggplot(plot_count, aes(Prop.x, Prop.y, color = Source_Type.x)) +
    theme_bw() +
    geom_density2d() +
    ylab("Proportion of spacers targeting plasmids") +
    xlab("Proportion of spacers targeting phages") +
    scale_color_manual(name = "Source", values = c("blue2", "red2"))+
    theme(legend.position = c(0.80,0.85)) +
    coord_cartesian(ylim = c(0,0.22), xlim = c(0,0.22))
p
ggsave(p, file = "Figures/Fig3_density_no_log.pdf", width = 11, height = 9, units = "cm")
