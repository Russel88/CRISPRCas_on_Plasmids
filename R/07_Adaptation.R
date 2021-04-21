library(ggplot2)

# Load data
load("Prepared.RData")

##### Distribution Cas ##### 
cas_plasmid_d <- cas_plasmid[, c("Acc", "Prediction", "Orphan", "Cell", "Genes")]
cas_host_d <- cas_host[cas_host$Derep, c("Acc", "Prediction", "Orphan", "Cell", "Genes")]

cas_plasmid_d$Origin <- "Plasmid"
cas_host_d$Origin <- "Chromosome"

cas <- rbind(cas_plasmid_d, cas_host_d)

# Cas1 or Cas2
cas$Adapt <- grepl("Cas[12]_", cas$Genes)

# Total
adaptTotal <- aggregate(Adapt ~ Origin, cas, mean)
adaptTotalN <- aggregate(Adapt ~ Origin, cas, length)
adaptTotal$Type <- "Total"
adaptTotal$N <- adaptTotalN$Adapt

# Per type
cas$Type <- gsub("-.*", "", cas$Prediction)
adaptType <- aggregate(Adapt ~ Origin + Type, cas, mean)
adaptTypeN <- aggregate(Adapt ~ Origin + Type, cas, length)
adaptType$N <- adaptTypeN$Adapt

# Plot
adaptType <- adaptType[adaptType$Type %in% c("I","II","III","IV","V"), ]

adaptAll <- rbind(adaptType, adaptTotal)
adaptAll$Type <- factor(adaptAll$Type, levels = c("Total","I","II","III","IV","V"))

adaptAll$N <- paste("n=", adaptAll$N)

p <- ggplot(adaptAll, aes(Origin, Adapt, label = N)) +
    theme_bw() +
    geom_bar(stat = "identity", fill = "darkred") +
    geom_text(aes(y = 0.03), hjust = 0) +
    facet_wrap(.~Type) +
    coord_flip(ylim = c(0,1)) +
    ylab(NULL) +
    xlab(NULL) +
    scale_y_continuous(labels = scales::percent)
p
ggsave(p, file = "Figures/Fig_adaptation.pdf", width = 18, height = 8, units = "cm")
