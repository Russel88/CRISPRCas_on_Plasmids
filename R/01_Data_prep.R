library(ape)
library(foreach)

# Load
cas_plasmid <- read.table("../Collect/Cas_plasmid.tab", sep="\t")
cas_host <- read.table("../Collect/Cas_host.tab", sep="\t")
colnames(cas_plasmid) <- c("Acc", "Operon", "Start", "End", "Prediction", "Interference", "Adaptation", "Best type", "Score", "Genes", "Pos", "Evals", "Cov_seq", "Cov_hmm")
colnames(cas_host) <- c("Acc", "Operon", "Start", "End", "Prediction", "Interference", "Adaptation", "Best type", "Score", "Genes", "Pos", "Evals", "Cov_seq", "Cov_hmm")

cris_plasmid <- read.table("../Collect/CRISPR_plasmid.tab", sep=" ")
cris_host <- read.table("../Collect/CRISPR_host.tab", sep=" ")
colnames(cris_plasmid) <- c("Acc", "CRISPR", "Start", "End", "Repeats", "Subtype", "Probability")
colnames(cris_host) <- c("Acc", "CRISPR", "Start", "End", "Repeats", "Subtype", "Probability")

criscas_plasmid <- read.table("../Collect/CRISPRCas_plasmid.tab")
criscas_host <- read.table("../Collect/CRISPRCas_host.tab")
colnames(criscas_plasmid) <- c("Acc", "CRISPR", "Operon", "Distance")
colnames(criscas_host) <- c("Acc", "CRISPR", "Operon", "Distance")

mobility <- read.table("../Collect/mobility.tab", sep = "\t")
colnames(mobility) <- c("Acc", "Inc", "Mobility")

# Remove duplicated plasmid's (NC_002128.1) host with different acc 
cris_host <- cris_host[cris_host$Acc != "NC_002695.1",]
cas_host <- cas_host[cas_host$Acc != "NC_002695.1",]
criscas_host <- criscas_host[criscas_host$Acc != "NC_002695.1",]

# Remove Mini-arrays
cris_plasmid <- cris_plasmid[cris_plasmid$Repeats > 2, ]
cris_host <- cris_host[cris_host$Repeats > 2, ]

criscas_plasmid <- criscas_plasmid[criscas_plasmid$CRISPR %in% cris_plasmid$CRISPR, ]
criscas_host <- criscas_host[criscas_host$CRISPR %in% cris_host$CRISPR, ]

# Low scoring probs are NA
cris_plasmid[is.na(cris_plasmid$Probability) | cris_plasmid$Probability < 0.75, "Subtype"] <- NA
cris_host[is.na(cris_host$Probability) | cris_host$Probability < 0.75, "Subtype"] <- NA

# Fix ambiguous
cas_plasmid[cas_plasmid$Operon == "NZ_CP031032.1@4", "Prediction"] <- "I-B"
cas_plasmid <- cas_plasmid[cas_plasmid$Operon != "NZ_CP031032.1@1", ]
cas_plasmid[cas_plasmid$Operon %in% c("KY563785.1@1",
                                      "NC_007959.1@1",
                                      "NC_017201.1@2",
                                      "NZ_CP016453.1@7",
                                      "NZ_CP017574.1@5",
                                      "NZ_CP028348.1@10",
                                      "NC_017958.1@11",
                                      "NZ_CP014527.1@2",
                                      "NZ_CP029356.1@6"), "Prediction"] <- "II-C"
cas_plasmid[cas_plasmid$Operon %in% c("NZ_CP040855.1@2",
                                      "NZ_CP047143.1@1"), "Prediction"] <- "II-A"
cas_plasmid[cas_plasmid$Prediction == "Ambiguous", "Prediction"] <- "II-C"

cas_host$Prediction <- as.character(cas_host$Prediction)
cas_host[cas_host$Operon == "NZ_AP017968.1@15", "Prediction"] <- "VI-C"
cas_host[cas_host$Operon == "NC_007643.1@18", "Prediction"] <- "III-D"
cas_host[cas_host$Operon %in% c("CP011381.2@9",
                                "NZ_CP011374.1@8"), "Prediction"] <- "III-B"
cas_host[cas_host$Operon == "NZ_CP018335.1@41", "Prediction"] <- "I-B"
cas_host[cas_host$Operon %in% c("NZ_LR215037.1@6",
                                "NZ_CP040908.1@"), "Prediction"] <- "II-A"
cas_host <- cas_host[cas_host$Operon != "NC_013385.1@28", ]
cas_host[cas_host$Operon == "NZ_CP046973.1@58", "Prediction"] <- "I-D"
cas_host <- cas_host[cas_host$Operon != "NZ_CP046973.1@57", ]
cas_host[cas_host$Operon == "NZ_CP011381.2@9", "Prediction"] <- "III-B"
cas_host[cas_host$Operon == "NZ_CP040272.1@32", "Prediction"] <- "V-J"
cas_host[cas_host$Operon == "NZ_CP043876.1@9 ", "Prediction"] <- "I-A"
cas_host[cas_host$Prediction == "Ambiguous", "Prediction"] <- "II-C"

# Remove variant info for V-F
cas_plasmid[grepl("V-F", cas_plasmid$Prediction), "Prediction"] <- "V-F"
cas_host[grepl("V-F", cas_host$Prediction), "Prediction"] <- "V-F"

cris_plasmid$Subtype <- as.character(cris_plasmid$Subtype)
cris_host$Subtype <- as.character(cris_host$Subtype)
cris_plasmid[grepl("V-F", cris_plasmid$Subtype), "Subtype"] <- "V-F"
cris_host[grepl("V-F", cris_host$Subtype), "Subtype"] <- "V-F"

# Cell data
cell_plsdb <- read.table("../Collect/cell_plsdb.tab", sep="\t", header=FALSE, stringsAsFactors = FALSE)
cell_archaea <- read.table("../Collect/cell_archaea.tab", sep="\t", header=FALSE, stringsAsFactors = FALSE)
cell <- rbind(cell_plsdb, cell_archaea)
colnames(cell) <- c("Cell", "Chromosomes", "Plasmids", "Unplaced")
cell$Cell <- gsub("GCA","GC",gsub("GCF","GC",cell$Cell))

tax <- read.table("../Collect/tax.tab", sep="\t")
colnames(tax) <- c("Cell", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax$Cell <- gsub("GCA","GC",gsub("GCF","GC",tax$Cell))
tax <- tax[!duplicated(tax$Cell), ]

drep_plsdb <- read.table("../Collect/PLSDB_drep.csv", sep=",", header=TRUE)
drep_archaea <- read.table("../Collect/Archaea_drep.csv", sep=",", header=TRUE)

repeat_sim <- read.table("../Collect/Repeat_sim.csv", sep=",", header=FALSE, stringsAsFactors = FALSE)
colnames(repeat_sim) <- c("First", "Second", "Sim")

size_plasmid <- read.table("../Collect/size_plasmid.tab")
colnames(size_plasmid) <- c("Acc", "Length")

size_host <- read.table("../Collect/size_host.tab")
colnames(size_host) <- c("Acc", "Length")

# New columns
drep_plsdb$Acc <- gsub("\\.fna", "", drep_plsdb$genome)
drep_archaea$Acc <- gsub("\\.fna", "", drep_archaea$genome)

# Remove unplaced contigs
cell <- cell[, 1:3]

# Expand cell data for plasmids and chromosomes separately
cell_plasmid <- cell[cell$Plasmids != "", c("Cell", "Plasmids")]

cell_plasmid <- foreach(k = cell_plasmid$Cell, .combine = rbind) %do% {
    data.frame(Cell = k,
               Acc = strsplit(cell_plasmid[cell_plasmid == k, "Plasmids"], ",")[[1]])
} 
cell_plasmid <- cell_plasmid[!duplicated(cell_plasmid$Acc), ]
cell_plasmid$Cell <- gsub("\\.[0-9]*$", "", cell_plasmid$Cell)


cell_host <- cell[cell$Chromosomes != "", c("Cell", "Chromosomes")]

cell_host <- foreach(k = cell_host$Cell, .combine = rbind) %do% {
    data.frame(Cell = k,
               Acc = strsplit(cell_host[cell_host == k, "Chromosomes"], ",")[[1]])
} 
cell_host <- cell_host[!duplicated(cell_host$Acc), ]
cell_host$Cell <- gsub("\\.[0-9]*$", "", cell_host$Cell)

# Add taxonomy
cell_plasmid <- merge(cell_plasmid, tax, by = "Cell", all.x = TRUE)
cell_host <- merge(cell_host, tax, by = "Cell", all.x = TRUE)

# Add clustering
drep <- rbind(drep_plsdb, drep_archaea)

cris_plasmid_final <- cris_plasmid
cas_plasmid_final <- cas_plasmid

cris_plasmid_final$Derep <- cris_plasmid_final$Acc %in% drep$Acc
cas_plasmid_final$Derep <- cas_plasmid_final$Acc %in% drep$Acc

# Determine if distant CRISPRs have similar repeats
repeat_sim <- repeat_sim[order(repeat_sim$Sim, decreasing = TRUE), ]

## Plasmids
criscas_plasmid_distant <- criscas_plasmid[criscas_plasmid$Distance > 1e4, ]
criscas_plasmid_near <- criscas_plasmid[criscas_plasmid$Distance <= 1e4, ]

criscas_plasmid_distant$Operon <- "Orphan"
for(k in seq_len(nrow(criscas_plasmid_distant))){
    tmp <- repeat_sim[repeat_sim$First == criscas_plasmid_distant[k, "CRISPR"] & 
                          repeat_sim$Second %in% criscas_plasmid_near[criscas_plasmid_near$Acc == criscas_plasmid_distant[k, "Acc"], "CRISPR"], ]
    if(nrow(tmp) > 0 & tmp[1, "Sim"] >= 0.85){
        criscas_plasmid_distant[k, "Operon"] <- as.character(criscas_plasmid_near[criscas_plasmid_near$CRISPR == tmp[1, "Second"], "Operon"])
    }
}
criscas_plasmid_final <- rbind(criscas_plasmid_near, criscas_plasmid_distant[criscas_plasmid_distant$Operon != "Orphan", ] )

cris_plasmid_final <- merge(cris_plasmid_final, criscas_plasmid_final[, c("CRISPR", "Operon", "Distance")], by = "CRISPR", all.x = TRUE)


## Hosts
criscas_host_distant <- criscas_host[criscas_host$Distance > 1e4, ]
criscas_host_near <- criscas_host[criscas_host$Distance <= 1e4, ]

criscas_host_distant$Operon <- "Orphan"
for(k in seq_len(nrow(criscas_host_distant))){
    tmp <- repeat_sim[repeat_sim$First == criscas_host_distant[k, "CRISPR"] & 
                          repeat_sim$Second %in% criscas_host_near[criscas_host_near$Acc == criscas_host_distant[k, "Acc"], "CRISPR"], ]
    if(nrow(tmp) > 0 & tmp[1, "Sim"] >= 0.85){
        criscas_host_distant[k, "Operon"] <- as.character(criscas_host_near[criscas_host_near$CRISPR == tmp[1, "Second"], "Operon"])
    }
}
criscas_host_final <- rbind(criscas_host_near, criscas_host_distant[criscas_host_distant$Operon != "Orphan", ] )

cris_host_final <- merge(cris_host, criscas_host_final[, c("CRISPR", "Operon", "Distance")], by = "CRISPR", all.x = TRUE)

# Add subtypes to CRISPRs
cris_plasmid_final <- merge(cris_plasmid_final, cas_plasmid_final[, c("Prediction", "Operon")], by = "Operon", all.x = TRUE)
cris_host_final <- merge(cris_host_final, cas_host[, c("Prediction", "Operon")], by = "Operon", all.x = TRUE)

# Remove untyped orphans
# cris_plasmid_final <- cris_plasmid_final[!(is.na(cris_plasmid_final$Prediction) & is.na(cris_plasmid_final$Subtype)), ]
# cris_host_final <- cris_host_final[!(is.na(cris_host_final$Prediction) & is.na(cris_host_final$Subtype)), ]

# Add CRISPRs to Cas
cas_host_final <- cas_host
cas_plasmid_final$Orphan <- ifelse(cas_plasmid_final$Operon %in% na.omit(cris_plasmid_final[, "Operon"]), 
                                   "CRISPR-Cas", "Orphan Cas")
cas_host_final$Orphan <- ifelse(cas_host_final$Operon %in% na.omit(cris_host_final[, "Operon"]), 
                                   "CRISPR-Cas", "Orphan Cas")

# Add cell information
cris_plasmid_final <- merge(cris_plasmid_final, cell_plasmid, by = "Acc", all.x = TRUE)
cas_plasmid_final <- merge(cas_plasmid_final, cell_plasmid, by = "Acc", all.x = TRUE)

cris_host_final <- merge(cris_host_final, cell_host, by = "Acc", all.x = TRUE)
cas_host_final <- merge(cas_host_final, cell_host, by = "Acc", all.x = TRUE)

# Add derep info to host tables
cris_host_final$Derep <- cris_host_final$Cell %in% unique(cell_plasmid[cell_plasmid$Acc %in% drep$Acc, "Cell"])
cas_host_final$Derep <- cas_host_final$Cell %in% unique(cell_plasmid[cell_plasmid$Acc %in% drep$Acc, "Cell"])

# Add derep info to cris cas distance
criscas_plasmid$Derep <- criscas_plasmid$Acc %in% drep$Acc
criscas_host$Derep <- criscas_host$Acc %in% cell_host[cell_host$Cell %in% unique(cell_plasmid[cell_plasmid$Acc %in% drep$Acc, "Cell"]), "Acc"]

# Add size
cas_plasmid_final <- merge(cas_plasmid_final, size_plasmid, by = "Acc")
cris_plasmid_final <- merge(cris_plasmid_final, size_plasmid, by = "Acc")

# Tables also with non-CRISPR-Cas plasmids and hosts
## Plasmids
cris_plasmid_count <- as.data.frame(table(cris_plasmid$Acc))
colnames(cris_plasmid_count) <- c("Acc", "CRISPRs")
cas_plasmid_count <- as.data.frame(table(cas_plasmid$Acc))
colnames(cas_plasmid_count) <- c("Acc", "Cas")
criscas_plasmid_count <- as.data.frame(table(cas_plasmid_final[cas_plasmid_final$Orphan == "CRISPR-Cas", "Acc"]))
colnames(criscas_plasmid_count) <- c("Acc", "CRISPRCas")
casorphan_plasmid_count <- as.data.frame(table(cas_plasmid_final[cas_plasmid_final$Orphan == "Orphan Cas", "Acc"]))
colnames(casorphan_plasmid_count) <- c("Acc", "Cas_Orphan")
crisorphan_plasmid_count <- as.data.frame(table(cris_plasmid_final[is.na(cris_plasmid_final$Prediction), "Acc"]))
colnames(crisorphan_plasmid_count) <- c("Acc", "CRISPR_Orphan")

all_plasmids <- size_plasmid
all_plasmids$Derep <- all_plasmids$Acc %in% drep$Acc
all_plasmids <- merge(all_plasmids, cell_plasmid, by = "Acc", all.x = TRUE)
all_plasmids <- merge(all_plasmids, cris_plasmid_count, by = "Acc", all.x = TRUE)
all_plasmids <- merge(all_plasmids, cas_plasmid_count, by = "Acc", all.x = TRUE)
all_plasmids <- merge(all_plasmids, criscas_plasmid_count, by = "Acc", all.x = TRUE)
all_plasmids <- merge(all_plasmids, casorphan_plasmid_count, by = "Acc", all.x = TRUE)
all_plasmids <- merge(all_plasmids, crisorphan_plasmid_count, by = "Acc", all.x = TRUE)

all_plasmids[is.na(all_plasmids$CRISPRs), "CRISPRs"] <- 0
all_plasmids[is.na(all_plasmids$Cas), "Cas"] <- 0
all_plasmids[is.na(all_plasmids$CRISPRCas), "CRISPRCas"] <- 0
all_plasmids[is.na(all_plasmids$Cas_Orphan), "Cas_Orphan"] <- 0
all_plasmids[is.na(all_plasmids$CRISPR_Orphan), "CRISPR_Orphan"] <- 0

## Hosts
cris_host_count <- as.data.frame(table(cris_host$Acc))
colnames(cris_host_count) <- c("Acc", "CRISPRs")
cas_host_count <- as.data.frame(table(cas_host$Acc))
colnames(cas_host_count) <- c("Acc", "Cas")
criscas_host_count <- as.data.frame(table(cas_host_final[cas_host_final$Orphan == "CRISPR-Cas", "Acc"]))
colnames(criscas_host_count) <- c("Acc", "CRISPRCas")
casorphan_host_count <- as.data.frame(table(cas_host_final[cas_host_final$Orphan == "Orphan Cas", "Acc"]))
colnames(casorphan_host_count) <- c("Acc", "Cas_Orphan")
crisorphan_host_count <- as.data.frame(table(cris_host_final[is.na(cris_host_final$Prediction), "Acc"]))
colnames(crisorphan_host_count) <- c("Acc", "CRISPR_Orphan")

all_hosts <- merge(size_host, cell_host, by = "Acc")
all_hosts$Derep <- all_hosts$Acc %in% cell_host[cell_host$Cell %in% unique(cell_plasmid[cell_plasmid$Acc %in% drep$Acc, "Cell"]), "Acc"]
all_hosts <- merge(all_hosts, cris_host_count, by = "Acc", all.x = TRUE)
all_hosts <- merge(all_hosts, cas_host_count, by = "Acc", all.x = TRUE)
all_hosts <- merge(all_hosts, criscas_host_count, by = "Acc", all.x = TRUE)
all_hosts <- merge(all_hosts, casorphan_host_count, by = "Acc", all.x = TRUE)
all_hosts <- merge(all_hosts, crisorphan_host_count, by = "Acc", all.x = TRUE)

all_hosts[is.na(all_hosts$CRISPRs), "CRISPRs"] <- 0
all_hosts[is.na(all_hosts$Cas), "Cas"] <- 0
all_hosts[is.na(all_hosts$CRISPRCas), "CRISPRCas"] <- 0
all_hosts[is.na(all_hosts$Cas_Orphan), "Cas_Orphan"] <- 0
all_hosts[is.na(all_hosts$CRISPR_Orphan), "CRISPR_Orphan"] <- 0

# Clean
cris_plasmid <- cris_plasmid_final
cas_plasmid <- cas_plasmid_final
cris_host <- cris_host_final
cas_host <- cas_host_final

rm(cris_plasmid_final, cas_plasmid_final, cris_host_final, cas_host_final, tmp, k, drep_archaea,
   criscas_host_distant, criscas_host_near, cell_archaea, cell_host, cell_plasmid, cell_plsdb, criscas_plasmid_final,
   criscas_host_final, criscas_plasmid_distant, criscas_plasmid_near,
   cris_host_count, cris_plasmid_count, cas_plasmid_count, cas_host_count, criscas_host_count, criscas_plasmid_count,
   casorphan_host_count, casorphan_plasmid_count, crisorphan_host_count, crisorphan_plasmid_count)

# Orphan plasmid CRISPRs
orphan_cris <- cris_plasmid[is.na(cris_plasmid$Operon) & !is.na(cris_plasmid$Cell), ]
orphan_cris_hosts <- cris_host[cris_host$Cell %in% orphan_cris$Cell, ]

# Trees
tree_bac <- read.tree("../Collect/plsdb.tree")
tree_arc <- read.tree("../Collect/archaea.tree")

tree_arc$tip.label <- gsub("GC[FA]_", "GC_", gsub("\\.fna$", "", tree_arc$tip.label))
tree_bac$tip.label <- gsub("GC[FA]_", "GC_", gsub("\\.fna$", "", tree_bac$tip.label))

tree_bac <- drop.tip(tree_bac, tip = tree_bac$tip.label[!tree_bac$tip.label %in% tax$Cell])
tree_arc <- drop.tip(tree_arc, tip = tree_arc$tip.label[!tree_arc$tip.label %in% tax$Cell])

tree <- bind.tree(tree_bac, tree_arc, where = "root")

tree <- drop.tip(tree, tip = which(duplicated(tree$tip.label)))

# Derep
cris_plasmid <- cris_plasmid[cris_plasmid$Derep, ]
cris_host <- cris_host[cris_host$Derep, ]
cas_plasmid <- cas_plasmid[cas_plasmid$Derep, ]
cas_host <- cas_host[cas_host$Derep, ]

# Add mobility
all_plasmids <- merge(all_plasmids, mobility, by = "Acc")

# Save
save(cris_plasmid, cris_host, cas_plasmid, cas_host, all_plasmids, all_hosts,
     cell, tax, tree, tree_bac, tree_arc,
     file = "Prepared.RData")
