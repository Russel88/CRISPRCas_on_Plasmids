library(ggplot2)

# Load
load("Prepared.RData")

# CRISPR Length
cris_plasmid$Origin <- "Plasmid"
cris_host$Origin <- "Chromosome"

cris <- rbind(cris_plasmid[, c("CRISPR","Repeats","Prediction", "Origin")], 
              cris_host[, c("CRISPR","Repeats","Prediction", "Origin")])

cris$Orphan <- ifelse(is.na(cris$Prediction), "Orphan", "Cas-associated")

p <- ggplot(cris, aes(Orphan, Repeats, fill = Origin)) +
    theme_bw() +
    geom_boxplot() +
    coord_flip() +
    xlab(NULL) +
    scale_y_log10()
p
ggsave(p, file = "Figures/CRISPR_Orphan.pdf", units = "cm", width = 18, height = 6)

# CRISPR Frequency
d_plasmids <- all_plasmids[all_plasmids$Derep & all_plasmids$CRISPRs > 0, ]
d_hosts <- all_hosts[all_hosts$Derep & all_hosts$CRISPRs > 0, ]

d_plasmids$Origin <- "Plasmid"
d_hosts$Origin <- "Chromosome"

dd <- rbind(d_plasmids[, c("CRISPRs", "Origin")], d_hosts[, c("CRISPRs", "Origin")])

p <- ggplot(dd, aes(CRISPRs)) +
    theme_bw() +
    geom_bar(position = "dodge", width = 0.9) +
    facet_grid(Origin~., scales = "free") +
    xlab("CRISPR arrays per replicon(n)") +
    ylab("Count") 
p
ggsave(p, file = "Figures/CRISPR_Frequency.pdf", units = "cm", width = 10, height = 6)
