load("Prepared.RData")

# Number of dereplicated plasmids
nrow(drep_plsdb)
sum(all_plasmids$Derep) - nrow(drep_plsdb)

# Number of CRISPR-Cas on plasmids
table(cas_plasmid[cas_plasmid$Derep,  "Orphan"])
nrow(cris_plasmid[is.na(cris_plasmid$Operon), ])

# % CRISPR-Cas in hosts
mean(apply(all_hosts[all_hosts$Derep & !is.na(all_hosts$Kingdom) & all_hosts$Kingdom != "Archaea", c("CRISPRs", "Cas")], 1, sum) > 0)
mean(apply(all_hosts[all_hosts$Derep & !is.na(all_hosts$Kingdom) & all_hosts$Kingdom == "Archaea", c("CRISPRs", "Cas")], 1, sum) > 0)

# CRISPR length - orphan/non-orphan
agg <- aggregate(Repeats ~ is.na(Operon), data = cris_plasmid[cris_plasmid$Derep, c("Repeats", "Operon")], mean)
agg[2,2] / agg[1,2]
library(MASS)          
summary(glm.nb(Repeats ~ is.na(Operon), data = cris_plasmid[cris_plasmid$Derep, ]))

# CRISPR length - plasmid/hosts
mean(cris_plasmid[cris_plasmid$Derep, "Repeats"]) / mean(cris_host[cris_host$Derep, "Repeats"])
summary(glm.nb(Repeats ~ Plasmid, data = data.frame(Repeats = c(cris_plasmid[cris_plasmid$Derep, "Repeats"],
                                                                cris_host[cris_host$Derep, "Repeats"]),
                                                    Plasmid = c(rep("Yes", sum(cris_plasmid$Derep)),
                                                                rep("No", sum(cris_host$Derep))))))

# Adaptataion genes presence (Only cas 1 and 2 !!!!!!!!!!)
adapt_p <- table(cas_plasmid[grepl("^I-", cas_plasmid$Prediction), "Adaptation"] == "0%")
adapt_h <- table(cas_host[grepl("^I-", cas_host$Prediction), "Adaptation"] == "0%")
adapt_p[1] / sum(adapt_p)
adapt_h[1] / sum(adapt_h)

# Multiple CRISPRs
mult_p <- table(all_plasmids[all_plasmids$Derep, "CRISPRs"])
mult_h <- table(all_hosts[all_hosts$Derep, "CRISPRs"])
sum(mult_p[3:length(mult_p)]) / sum(mult_p[2:length(mult_p)])
sum(mult_h[3:length(mult_h)]) / sum(mult_h[2:length(mult_h)])

# Proteobacterial plasmids
nr_known_host_p <- all_plasmids[all_plasmids$Derep & !is.na(all_plasmids$Phylum),]

nrow(nr_known_host_p[nr_known_host_p$Phylum == "Proteobacteria", ])/nrow(nr_known_host_p)
nrow(nr_known_host_p[nr_known_host_p$Class %in% c("Gammaproteobacteria", "Bacilli", "Alphaproteobacteria"), ])/nrow(nr_known_host_p)
