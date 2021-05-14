library(ggplot2)

# Load data
load("Prepared.RData")

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
ggsave(p, file = "Figures/Fig1_subtypes.pdf", width = 12, height = 8, units = "cm")
write.csv(criscas_a[criscas_a$Origin == "Plasmid", c("Type", "SubType", "Prediction")], file = "Tables/Fig1_subtypes.csv", quote = FALSE, row.names = FALSE)

# Only chromosomes
p <- ggplot(criscas_a[criscas_a$Origin == "Chromosome", ], aes(SubType, Type, size = Prediction, color = Type)) +
    geom_point() +
    geom_text(aes(label = Prediction), 
              colour = "black", 
              size = 3,
              vjust = -1.1) +
    ylab("Type") +
    xlab("Subtype") +
    theme_bw() +
    scale_size_area() +
    theme(legend.position = "none",
          panel.border = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_manual(values = c("grey", "cyan", "gold", "purple", "red2", "blue2", "green2"))
p
ggsave(p, file = "Figures/Fig1_subtypes_chromosome.pdf", width = 12, height = 8, units = "cm")
write.csv(criscas_a[criscas_a$Origin == "Chromosome", c("Type", "SubType", "Prediction")], file = "Tables/Fig1_subtypes_chromosomes.csv", quote = FALSE, row.names = FALSE)

##### Prevalence #####
prev <- data.frame(Origin = c(rep("Plasmid", 3),
                              rep("Chromosome", 3)),
                   Type = rep(c("CRISPR-Cas", "Orphan CRISPR", "Orphan Cas operon"), 2),
                   Prevalence = c(aggregate(CRISPRCas ~ 1, all_plasmids[all_plasmids$Derep, ], function(x) mean(x > 0))[[1]],
                                  aggregate(CRISPR_Orphan ~ 1, all_plasmids[all_plasmids$Derep, ], function(x) mean(x > 0))[[1]],
                                  aggregate(Cas_Orphan ~ 1, all_plasmids[all_plasmids$Derep, ], function(x) mean(x > 0))[[1]],
                                  aggregate(CRISPRCas ~ 1, all_hosts[all_hosts$Derep, ], function(x) mean(x > 0))[[1]],
                                  aggregate(CRISPR_Orphan ~ 1, all_hosts[all_hosts$Derep, ], function(x) mean(x > 0))[[1]],
                                  aggregate(Cas_Orphan ~ 1, all_hosts[all_hosts$Derep, ], function(x) mean(x > 0))[[1]]))

p <- ggplot(prev, aes(Type, Prevalence, fill = Type)) +
    theme_bw() +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = scales::percent) +
    facet_grid(. ~ Origin, scales = "free") +
    coord_flip()
p
ggsave(p, file = "Figures/Fig1_prevalence_1.pdf", width = 20, height = 6, units = "cm")
write.csv(prev[, c("Type", "Origin", "Prevalence")], file = "Tables/Fig1_prevalence_1.csv", quote = FALSE, row.names = FALSE)

prevN <- data.frame(Origin = c(rep("Plasmid", 3),
                              rep("Chromosome", 3)),
                   Type = rep(c("CRISPR-Cas", "Orphan CRISPR", "Orphan Cas operon"), 2),
                   Prevalence = c(aggregate(CRISPRCas/Length ~ 1, all_plasmids[all_plasmids$Derep, ], function(x) mean(x*1e6))[[1]],
                                  aggregate(CRISPR_Orphan/Length ~ 1, all_plasmids[all_plasmids$Derep, ], function(x) mean(x*1e6))[[1]],
                                  aggregate(Cas_Orphan/Length ~ 1, all_plasmids[all_plasmids$Derep, ], function(x) mean(x*1e6))[[1]],
                                  aggregate(CRISPRCas/Length ~ 1, all_hosts[all_hosts$Derep, ], function(x) mean(x*1e6))[[1]],
                                  aggregate(CRISPR_Orphan/Length ~ 1, all_hosts[all_hosts$Derep, ], function(x) mean(x*1e6))[[1]],
                                  aggregate(Cas_Orphan/Length ~ 1, all_hosts[all_hosts$Derep, ], function(x) mean(x*1e6))[[1]]))

p <- ggplot(prevN, aes(Origin, Prevalence, fill = Type)) +
    theme_bw() +
    geom_bar(stat="identity") +
    coord_flip() +
    ylab("# per Mbp")
p
ggsave(p, file = "Figures/Fig1_prevalence_2.pdf", width = 16, height = 5, units = "cm")
write.csv(prevN[, c("Type", "Origin", "Prevalence")], file = "Tables/Fig1_prevalence_2.csv", quote = FALSE, row.names = FALSE)
