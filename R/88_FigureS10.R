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

# Summarise by CRISPR
all_spacers_unq <- aggregate(Spacer ~ CRISPR + Source_Type + Target_Type, data = all_spacers, function(x) length(unique(x)))

# Add N repeats to df
cris_plasmid <- cris_plasmid[, c("CRISPR", "Repeats", "Acc", "Class", "Prediction")]
cris_host <- cris_host[, c("CRISPR", "Repeats", "Acc", "Class", "Prediction")]

cris_plasmid$Origin <- "Plasmid"
cris_host$Origin <- "Chromosome"

cris <- rbind(cris_plasmid, cris_host)

cris_target_pl <- cris
cris_target_ph <- cris

cris_target_pl$Target_Type <- "Plasmid"
cris_target_ph$Target_Type <- "Phage"

cris_target <- rbind(cris_target_pl, cris_target_ph)

cris_target$Merge <- paste(cris_target$CRISPR, cris_target$Target_Type)
all_spacers_unq$Merge <- paste(all_spacers_unq$CRISPR, all_spacers_unq$Target_Type)

all_spacers_merged <- merge(all_spacers_unq, cris_target, by = "Merge", all = TRUE)
all_spacers_merged[is.na(all_spacers_merged$Spacer), "Spacer"] <- 0

# Subtype
all_spacers_merged$Prediction <- as.character(all_spacers_merged$Prediction)
all_spacers_merged[grepl("Hybrid", all_spacers_merged$Prediction), "Prediction"] <- "Hybrid"
all_spacers_merged$Spacers <- all_spacers_merged$Repeats - 1

subtype <- aggregate(Spacer ~ Prediction + Origin + Target_Type.y, all_spacers_merged, sum) 
subtype_rp <- aggregate(Spacers ~ Prediction + Origin + Target_Type.y, all_spacers_merged, sum) 

subtype$Spacers <- subtype_rp$Spacers
subtype$Fraction <- subtype$Spacer / subtype$Spacers

subtype$Prediction <- paste0(subtype$Prediction, " (n=", subtype$Spacers, ")")
subtype$Prediction <- as.factor(subtype$Prediction)
subtype$Prediction <- factor(subtype$Prediction, levels = rev(levels(subtype$Prediction)))

p <- ggplot(subtype, aes(Prediction, Fraction, fill = Target_Type.y)) +
    theme_bw() +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5, preserve = "single"), width = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    facet_grid(Origin~., scales = "free", space = "free") +
    coord_flip() +
    theme(legend.position = c(0.8, 0.93)) +
    ylab("Fraction of spacers with a match") +
    xlab("Subtype") +
    scale_fill_discrete(name = "Target")
p
ggsave(p, file = "Figures/Spacer_fraction_subtypes.pdf", units = "cm", width = 10, height = 18)

# class
classes <- aggregate(Spacer ~ Class + Origin + Target_Type.y, all_spacers_merged, sum) 
classes_rp <- aggregate(Spacers ~ Class + Origin + Target_Type.y, all_spacers_merged, sum) 

classes$Spacers <- classes_rp$Spacers
classes$Fraction <- classes$Spacer / classes$Spacers

classes$Class <- paste0(classes$Class, " (n=", classes$Spacers, ")")
classes$Class <- as.factor(classes$Class)
classes$Class <- factor(classes$Class, levels = rev(levels(classes$Class)))

p <- ggplot(classes, aes(Class, Fraction, fill = Target_Type.y)) +
    theme_bw() +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5, preserve = "single"), width = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    facet_grid(Origin~., scales="free", space="free") +
    coord_flip() +
    theme(legend.position = c(0.8, 0.93)) +
    ylab("Fraction of spacers with a match") +
    xlab("Class") +
    scale_fill_discrete(name = "Target")
p
ggsave(p, file = "Figures/Spacer_fraction_classes.pdf", units = "cm", width = 13, height = 20)
