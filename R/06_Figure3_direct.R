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

# Remove those from untyped orphans
orphans <- c(as.character(cris_plasmid[is.na(cris_plasmid$Prediction) & is.na(cris_plasmid$Subtype), "CRISPR"]),
             as.character(cris_host[is.na(cris_host$Prediction) & is.na(cris_host$Subtype), "CRISPR"]))

all_spacers <- all_spacers[!gsub("@[0-9]*$", "", all_spacers$Spacer) %in% orphans, ]

# Only count spacers once in each target group (CRISPR is used as a stand-in variable for de-replicating)
all_spacers_unq <- aggregate(CRISPR ~ Spacer + Source_Type + Target_Type, data = all_spacers, function(x) 1)

pl_spacers_unq <- all_spacers_unq[all_spacers_unq$Source_Type == "Plasmid" & all_spacers_unq$Target_Type != "Self", ]
hs_spacers_unq <- all_spacers_unq[all_spacers_unq$Source_Type == "Chromosome", ]

# For direct comparison
these <- intersect(cris_host$Cell, cris_plasmid$Cell)
cris_host <- cris_host[cris_host$Cell %in% these, ]
cris_plasmid <- cris_plasmid[cris_plasmid$Cell %in% these, ]

# Total spacers
sum(cris_plasmid[!cris_plasmid$CRISPR %in% orphans, "Repeats"] - 1)
sum(cris_host[!cris_host$CRISPR %in% orphans, "Repeats"] - 1)

# Subset direct
pl_spacers_unq <- pl_spacers_unq[gsub("@[0-9]*$", "", pl_spacers_unq$Spacer) %in% cris_plasmid$CRISPR, ]
hs_spacers_unq <- hs_spacers_unq[gsub("@[0-9]*$", "", hs_spacers_unq$Spacer) %in% cris_host$CRISPR, ]

# Plot
pdf(file = "Figures/Fig3_plasmid_euler_direct.pdf", width = 6, height = 4)
plot(euler(list(Phage = pl_spacers_unq[pl_spacers_unq$Target_Type == "Phage", "Spacer"],
                Plasmid = pl_spacers_unq[pl_spacers_unq$Target_Type == "Plasmid", "Spacer"])),
     fill = c("orange", "cyan3"), alpha = 0.8,
     quantities = list(type = c("counts", "percent")))
dev.off()

pdf(file = "Figures/Fig3_host_euler_direct.pdf", width = 6, height = 4)
plot(euler(list(Phage = hs_spacers_unq[hs_spacers_unq$Target_Type == "Phage", "Spacer"],
                Plasmid = hs_spacers_unq[hs_spacers_unq$Target_Type == "Plasmid", "Spacer"])),
     fill = c("orange", "cyan3"), alpha = 0.8,
     quantities = list(type = c("counts", "percent")))
dev.off()

# Per subtype and class
all_spacers_unq2 <- all_spacers_unq[all_spacers_unq$Target_Type != "Self", ]

all_spacers_unq2$CRISPR <- gsub("@[0-9]*$", "", all_spacers_unq2$Spacer)
all_spacers_unq2 <- merge(all_spacers_unq2, rbind(cris_host[, c("CRISPR", "Subtype", "Prediction", "Class")], 
                                                  cris_plasmid[, c("CRISPR", "Subtype", "Prediction", "Class")]), by = "CRISPR")
all_spacers_unq2$SubType_final <- ifelse(!is.na(all_spacers_unq2$Prediction), 
                                       as.character(all_spacers_unq2$Prediction), 
                                       as.character(all_spacers_unq2$Subtype))

all_spacers_unq2$Source_Type <- factor(all_spacers_unq2$Source_Type, levels = c("Plasmid", "Chromosome"))

# Type
all_spacers_unq_type <- aggregate(Spacer ~ SubType_final + Source_Type + Target_Type, data = all_spacers_unq2, function(x) length(unique(x)))

pl_types <- aggregate(Spacer ~ SubType_final, all_spacers_unq_type[all_spacers_unq_type$Source_Type == "Plasmid", ], sum)
hs_types <- aggregate(Spacer ~ SubType_final, all_spacers_unq_type[all_spacers_unq_type$Source_Type == "Chromosome", ], sum)

types <- union(pl_types[rev(order(pl_types$Spacer)), "SubType_final"][1:6],
          hs_types[rev(order(hs_types$Spacer)), "SubType_final"][1:6])

all_spacers_unq_type$SubType <- ifelse(all_spacers_unq_type$SubType_final %in% types, all_spacers_unq_type$SubType_final, "Other")
all_spacers_unq_type$SubType <- relevel(as.factor(all_spacers_unq_type$SubType), "Other")

p <- ggplot(all_spacers_unq_type, aes(SubType, Spacer, fill = Target_Type)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_grid(. ~ Source_Type, scales = "free") +
    xlab(NULL) +
    ylab("Unique spacers") +
    scale_fill_manual(values = c("orange", "cyan3"))
p
ggsave(p, file = "Figures/Fig3_subtypes_direct.pdf", width = 14, height = 7, units = "cm")
write.csv(all_spacers_unq_type[, c("SubType_final", "Spacer", "Target_Type", "Source_Type")], file = "Tables/Fig3_subtypes_direct.csv", quote = FALSE, row.names = FALSE)

# Class
all_spacers_unq_class <- aggregate(Spacer ~ Class + Source_Type + Target_Type, data = all_spacers_unq2, function(x) length(unique(x)))

pl_class <- aggregate(Spacer ~ Class, all_spacers_unq_class[all_spacers_unq_class$Source_Type == "Plasmid", ], sum)
hs_class <- aggregate(Spacer ~ Class, all_spacers_unq_class[all_spacers_unq_class$Source_Type == "Chromosome", ], sum)

classes <- union(pl_class[rev(order(pl_class$Spacer)), "Class"][1:6],
               pl_class[rev(order(pl_class$Spacer)), "Class"][1:6])

all_spacers_unq_class$Classes <- ifelse(all_spacers_unq_class$Class %in% classes, as.character(all_spacers_unq_class$Class), "Other")
all_spacers_unq_class$Classes <- relevel(as.factor(all_spacers_unq_class$Classes), "Other")

p <- ggplot(all_spacers_unq_class, aes(Classes, Spacer, fill = Target_Type)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_grid(. ~ Source_Type, scales = "free") +
    xlab(NULL) +
    ylab("Unique spacers") +
    scale_fill_manual(values = c("orange", "cyan3"))
p
ggsave(p, file = "Figures/Fig3_class_direct.pdf", width = 14, height = 7, units = "cm")
write.csv(all_spacers_unq_class[, c("Class", "Spacer", "Target_Type", "Source_Type")], file = "Tables/Fig3_class_direct.csv", quote = FALSE, row.names = FALSE)
