# Load
criscas_plasmid <- read.table("../Collect/CRISPRCas_plasmid.tab")
criscas_host <- read.table("../Collect/CRISPRCas_host.tab")
colnames(criscas_plasmid) <- c("Acc", "CRISPR", "Operon", "Distance")
colnames(criscas_host) <- c("Acc", "CRISPR", "Operon", "Distance")

criscas_all <- rbind(criscas_host, criscas_plasmid)

load("Prepared.RData")

cris_plasmid$Type <- "Plasmid"
cris_host$Type <- "Chromosome"

criscas <- rbind(cris_plasmid[, c("CRISPR","Prediction", "Type")], cris_host[, c("CRISPR","Prediction", "Type")])

criscas_merged <- merge(criscas_all, criscas, by = "CRISPR")

criscas_merged$Typed <- ifelse(is.na(criscas_merged$Prediction), "No", "Yes")

p <- ggplot(criscas_merged, aes(Distance+1, fill = Typed)) +
    theme_bw() +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = 1e4) +
    xlab("Cas - CRISPR distance (bp)") +
    ylab("Count") +
    scale_x_log10(labels = COEF::fancy_scientific) +
    facet_grid(Type~., scales = "free")
p

ggsave(p, file = "Figures/CRISPRCas_Distance.pdf", units = "cm", width = 12, height = 10)
