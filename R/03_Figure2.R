library(ggplot2)

# Load data
load("Prepared.RData")

# Combine
cas_sub <- cas_plasmid_d[, c("Acc", "Length", "Orphan", "Prediction" ,"Phylum", "Class", "Order", "Family", "Genus")]
cris_sub <- cris_plasmid_d[is.na(cris_plasmid_d$Prediction), c("Acc", "Length", "Prediction", "Phylum", "Class", "Order", "Family", "Genus")]
cris_sub$Orphan <- "Orphan CRISPR"

criscas <- rbind(cas_sub, cris_sub)
colnames(criscas)[3] <- "Group"

##### Mobility which taxa ######
mob <- read.table("plasmid_mob.tab")
colnames(mob) <- c("Acc", "Inc", "Mobility")
mob$Acc <- gsub("\\.fna","",mob$Acc)
mob <- mob[mob$Acc %in% drep$Acc, ]

sum(grepl("Inc", mob$Inc)) / nrow(mob)

criscas <- merge(criscas, mob, by = "Acc", all.x = TRUE)

p <- ggplot(criscas, aes(Phylum, fill = Mobility)) +
    theme_bw() +
    geom_bar() +
    coord_flip()
p

p <- ggplot(criscas, aes(Phylum, fill = grepl("Inc", Inc))) +
    theme_bw() +
    geom_bar() +
    coord_flip()
p

p <- ggplot(criscas, aes(Phylum, fill = grepl("rep_cluster", Inc))) +
    theme_bw() +
    geom_bar() +
    coord_flip()
p

mobt <- merge(mob, all_plasmids, by = "Acc", all.x = TRUE)

p <- ggplot(mobt, aes(Phylum, fill = Mobility)) +
    theme_bw() +
    geom_bar(position="fill") +
    coord_flip()
p

p <- ggplot(mobt, aes(Phylum, fill = grepl("Inc", Inc))) +
    theme_bw() +
    geom_bar(position="fill") +
    coord_flip()
p

p <- ggplot(mobt, aes(Phylum, fill = grepl("rep_cluster", Inc))) +
    theme_bw() +
    geom_bar(position="fill") +
    coord_flip()
p

#### Mobility #### 
mobt <- mobt[!mobt$Acc %in% criscas$Acc, ]

criscas_prot <- criscas[criscas$Phylum == "Proteobacteria", ]
mobt_prot <- mobt[mobt$Phylum == "Proteobacteria", ]
criscas_mob <- aggregate(Acc ~ Group + Mobility + Genus, data = criscas_prot, length)
plsdb_mob <- aggregate(Acc ~ Mobility + Genus, data = mobt_prot, length)
plsdb_mob$Group <- "No CRISPR or Cas"
plsdb_mob <- plsdb_mob[,c("Group", "Mobility", "Acc", "Genus")]

criscas_mob <- rbind(criscas_mob, plsdb_mob)
#criscas_mob <- criscas_mob[grepl("proteo", criscas_mob$Class), ]
#View(aggregate(Acc ~ Genus, mobt_prot, length))
#criscas_mob <- criscas_mob[criscas_mob$Genus %in% c("Escherichia", "Klebsiella", "Salmonella", 
#                                                    "Acinetobacter", "Rhizobium","Enterobacter"), ]



criscas_mob_agg <- aggregate(Acc ~Genus, data = criscas_mob, sum) 

criscas_mob$Group <- factor(criscas_mob$Group, levels = c("No CRISPR or Cas", "Orphan CRISPR", "Orphan Cas", "CRISPR-Cas"))

criscas_mob <- criscas_mob[criscas_mob$Genus %in% criscas_mob_agg[criscas_mob_agg$Acc > 50, "Genus"],]

p <- ggplot(criscas_mob, aes(Group, Acc, fill = Mobility)) +
    theme_bw() +
    geom_bar(stat="identity", position = "fill") +
    coord_flip() +
    scale_y_continuous(labels=scales::percent) +
    ylab("Proportion") +
    facet_wrap(~Genus)
p
ggsave(p, file = "Figures/Fig2b.pdf", width = 20, height = 6, units = "cm")

##### Inc ######
mob$Inc <- as.character(mob$Inc)
criscas$Inc <- as.character(criscas$Inc)
criscas$Prediction <- as.character(criscas$Prediction)
criscas[is.na(criscas$Prediction), "Prediction"] <- "Orphan"

all_mob <- data.frame(table(unlist(lapply(mob$Inc, function(x) strsplit(x, ",")))))
cc_mob <- data.frame(table(unlist(lapply(criscas[criscas$Group == "CRISPR-Cas","Inc"], function(x) strsplit(x, ",")))))

cc_inc_type <- lapply(unique(criscas$Prediction), function(x) data.frame(table(unlist(lapply(criscas[criscas$Prediction == x,"Inc"], function(x) strsplit(x, ",")))),
                                                                         Subtype = x))

cc_inc_type <- do.call(rbind, cc_inc_type)
cc_inc_type <- cc_inc_type[grepl("Inc", cc_inc_type$Var1), ]
View(xtabs(Freq ~ as.character(Subtype) + as.character(Var1), cc_inc_type))

sum(all_mob[grepl("Inc", all_mob$Var1),  "Freq"])/sum(all_mob$Freq)

#cas_plasmid_d[cas_plasmid_d$Acc %in% criscas[grepl("IncF", criscas$Inc), "Acc"],]
#cas_plasmid_d[cas_plasmid_d$Acc %in% criscas[grepl("IncH", criscas$Inc), "Acc"],]

all_mob <- all_mob[!grepl("rep_cluster", all_mob$Var1), ]
cc_mob <- cc_mob[!grepl("rep_cluster", cc_mob$Var1), ]
all_mob <- all_mob[all_mob$Var1 != "-", ]
cc_mob <- cc_mob[cc_mob$Var1 != "-", ]

all_mob$Freq <- all_mob$Freq / sum(all_mob$Freq)
cc_mob$Freq <- cc_mob$Freq / sum(cc_mob$Freq)

all_mob$Group <- "All"
cc_mob$Group <- "CRISPR-Cas"
m_mob <- rbind(all_mob, cc_mob)
#m_mob <- m_mob[!grepl("rep_cluster", m_mob$Var1), ]
#m_mob <- m_mob[m_mob$Var1 != "-", ]

p <- ggplot(m_mob, aes(Var1, Freq)) +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines")) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_grid(.~Group, scales = "free") +
    xlab("Inc group") +
    ylab("Frequency") +
    scale_y_continuous(labels=scales::percent)
p
ggsave(p, file = "Figures/Fig2c.pdf", width = 16, height = 12, units = "cm")

all_mob <- data.frame(table(unlist(lapply(mob$Inc, function(x) strsplit(x, ",")))))
cc_mob <- data.frame(table(unlist(lapply(criscas[criscas$Group == "CRISPR-Cas","Inc"], function(x) strsplit(x, ",")))))
all_mob$Freq <- all_mob$Freq / sum(all_mob$Freq)
cc_mob$Freq <- cc_mob$Freq / sum(cc_mob$Freq)
cc_mob <- cc_mob[cc_mob$Freq > 0.0066, ]
m_mob <- merge(all_mob, cc_mob, by = "Var1", all.x = TRUE)
m_mob[is.na(m_mob$Freq.x), "Freq.x"] <- 0
m_mob[is.na(m_mob$Freq.y), "Freq.y"] <- 0

m_mob[log2(m_mob$Freq.y/m_mob$Freq.x) > 5,]

hist(log2(m_mob$Freq.y / m_mob$Freq.x))

##### Plasmid size #####
criscas <- rbind(cas_sub, cris_sub)
colnames(criscas)[3] <- "Group"

size_plasmid_no <- size_plasmid[!size_plasmid$Acc %in% criscas$Acc, ]
size_plasmid_no$Group <- "No CRISPR or Cas"

criscas_len <- rbind(criscas[, c("Acc", "Length", "Group")], size_plasmid_no)
criscas_len$Group <- factor(criscas_len$Group, levels = rev(c("No CRISPR or Cas", "Orphan CRISPR", "Orphan Cas", "CRISPR-Cas")))

p <- ggplot(criscas_len, aes(Length)) +
    theme_bw() +
    geom_density(alpha = 0.6, fill = "grey") +
    ylab("Density") +
    xlab("Length (bp)") +
    scale_x_log10(labels = COEF::fancy_scientific) +
    facet_grid(Group ~ .)
p
ggsave(p, file = "Figures/Fig2a.pdf", width = 8, height = 14, units = "cm")



