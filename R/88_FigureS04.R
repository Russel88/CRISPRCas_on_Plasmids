library(ggplot2)

# Load data
load("Prepared.RData")

# For direct comparison
these <- intersect(c(cris_host$Cell, cas_host$Cell), c(cris_plasmid$Cell, cas_plasmid$Cell))
cris_host <- cris_host[cris_host$Cell %in% these, ]
cris_plasmid <- cris_plasmid[cris_plasmid$Cell %in% these, ]
cas_host <- cas_host[cas_host$Cell %in% these, ]
cas_plasmid <- cas_plasmid[cas_plasmid$Cell %in% these, ]

##### Distribution Cas ##### 
cas_plasmid_d <- cas_plasmid[, c("Acc", "Prediction", "Orphan", "Cell")]
cas_host_d <- cas_host[, c("Acc", "Prediction", "Orphan", "Cell")]

cas_plasmid_d$Origin <- "Plasmid"
cas_host_d$Origin <- "Chromosome"

cas <- rbind(cas_plasmid_d, cas_host_d)

##### Distribution CRISPR ##### 
cris_plasmid_d <- cris_plasmid[, c("Acc", "Repeats", "Distance", "Prediction", "Cell", "Subtype")]
cris_host_d <- cris_host[, c("Acc", "Repeats", "Distance", "Prediction", "Cell", "Subtype")]

cris_plasmid_d$Origin <- "Plasmid"
cris_host_d$Origin <- "Chromosome"

cris <- rbind(cris_plasmid_d, cris_host_d)
cris$Prediction <- as.character(cris$Prediction)

cris$Orphan <- is.na(cris$Prediction)
cris$Prediction <- ifelse(is.na(cris$Subtype), "Other", cris$Subtype)

##### Distribution CRISPR-Cas ##### 
cas_sub <- cas[, c("Prediction", "Origin", "Cell", "Orphan", "Acc")]
cris_sub <- cris[cris$Orphan, c("Prediction", "Origin", "Cell", "Acc")]
cris_sub$Orphan <- "Orphan CRISPR"

criscas <- rbind(cas_sub, cris_sub)
colnames(criscas)[4] <- "System"

criscas$Prediction <- as.character(criscas$Prediction)
criscas[grepl("Hybrid", criscas$Prediction), "Prediction"] <- "Other"

criscas$Type <- gsub("-.*", "", criscas$Prediction)
criscas$SubType <- gsub(".*-", "", criscas$Prediction)

criscas_a <- aggregate(Prediction ~ Type + SubType + Origin, data = criscas, function(x) sum(!is.na(x)))
criscas_aT <- aggregate(Prediction ~ Type + Origin, data = criscas, function(x) sum(!is.na(x)))
criscas_aT$SubType <- "Total"
criscas_a <- rbind(criscas_a, criscas_aT)

criscas_a$SubType <- factor(criscas_a$SubType, levels = c("A", "A1", "A2", "A3", "B", "B1", "C", "D", "E", "F", "F_T", "G", "J", "K", "Other", "Total"))
criscas_a$Type <- factor(criscas_a$Type, levels = c("Other", "VI", "V", "IV", "III", "II", "I"))

ns <- aggregate(Prediction ~ Origin, criscas_a[criscas_a$SubType == "Total", ], sum)

# Only plasmids
p <- ggplot(criscas_a[criscas_a$Origin == "Plasmid", ], aes(SubType, Type, size = Prediction, color = Type)) +
    geom_point() +
    geom_text(aes(label = Prediction), 
              colour = "black", 
              size = 3,
              vjust = -1.7) +
    ylab("Type") +
    xlab("Subtype") +
    theme_bw() +
    scale_size_area(max_size = 10) +
    theme(legend.position = "none",
          panel.border = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_manual(values = c("grey", "gold", "purple", "red2", "blue2", "green2"))
p
ggsave(p, file = "Figures/Fig1_subtypes_pl_direct.pdf", width = 12, height = 8, units = "cm")

# Only chromosomes
p <- ggplot(criscas_a[criscas_a$Origin == "Chromosome", ], aes(SubType, Type, size = Prediction, color = Type)) +
    geom_point() +
    geom_text(aes(label = Prediction), 
              colour = "black", 
              size = 3,
              vjust = -1.7) +
    ylab("Type") +
    xlab("Subtype") +
    theme_bw() +
    scale_size_area(max_size = 10) +
    theme(legend.position = "none",
          panel.border = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_manual(values = c("grey", "cyan", "gold", "purple", "red2", "blue2", "green2"))
p
ggsave(p, file = "Figures/Fig1_subtypes_hs_direct.pdf", width = 12, height = 9, units = "cm")

# Bar
criscas_agg <- aggregate(Acc ~ Prediction + Origin, data = criscas, function(x) sum(!is.na(x)))
criscas_aggT <- aggregate(Acc ~ Origin, data = criscas, function(x) sum(!is.na(x)))
criscas_agg <- merge(criscas_agg, criscas_aggT, by = "Origin")
criscas_agg$Perc <- criscas_agg$Acc.x / criscas_agg$Acc.y

criscas_agg$Prediction <- as.factor(criscas_agg$Prediction)
criscas_agg$Prediction <- factor(criscas_agg$Prediction, levels = rev(levels(criscas_agg$Prediction)))
criscas_agg$Prediction <- relevel(criscas_agg$Prediction, "Other")

p <- ggplot(criscas_agg, aes(Prediction, Acc.x, fill = Origin)) +
    theme_bw() +
    geom_bar(stat = "identity", position = position_dodge(0.5, preserve = "single"), width = 0.5) +
    coord_flip() +
    ylab("#") +
    xlab("Subtype")
p
ggsave(p, file = "Figures/Fig1_subtypes_direct_bar.pdf", width = 12, height = 12, units = "cm")
