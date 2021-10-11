library(ggplot2)

# Load
load("Prepared.RData")

# CRISPR Length
cris_plasmid$Origin <- "Plasmid"
cris_host$Origin <- "Chromosome"

cris <- rbind(cris_plasmid[, c("CRISPR","Repeats","Prediction", "Origin")], 
              cris_host[, c("CRISPR","Repeats","Prediction", "Origin")])

cris$Orphan <- ifelse(is.na(cris$Prediction), "Orphan", "Cas-associated")

cris$Origin <- factor(cris$Origin, levels = c("Plasmid", "Chromosome"))

p <- ggplot(cris, aes(Orphan, Repeats, fill = Origin)) +
    theme_bw() +
    geom_boxplot() +
    xlab(NULL) +
    scale_y_log10() +
    scale_fill_manual(values = c("purple", "darkgreen"))
p
ggsave(p, file = "Figures/CRISPR_Orphan.pdf", units = "cm", width = 12, height = 6)

# Stats
library(MASS)          
fit <- glm.nb(Repeats ~ Orphan*Origin, data = cris)
summary(fit)

# CRISPR Frequency
d_plasmids <- all_plasmids[all_plasmids$Derep & all_plasmids$CRISPRs > 0, ]
d_hosts <- all_hosts[all_hosts$Derep & all_hosts$CRISPRs > 0, ]

d_plasmids$Origin <- "Plasmid"
d_hosts$Origin <- "Chromosome"

dd <- rbind(d_plasmids[, c("CRISPRs", "Origin")], d_hosts[, c("CRISPRs", "Origin")])

dd$Origin <- factor(dd$Origin, levels = c("Plasmid", "Chromosome"))

p <- ggplot(dd, aes(CRISPRs)) +
    theme_bw() +
    geom_bar(position = "dodge", width = 0.9) +
    facet_grid(Origin~., scales = "free") +
    xlab("CRISPR arrays per replicon(n)") +
    ylab("Count") 
p
ggsave(p, file = "Figures/CRISPR_Frequency.pdf", units = "cm", width = 10, height = 6)
