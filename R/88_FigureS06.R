library(ggplot2)
library(ape)

# Load data
load("Prepared.RData")

##### Prevalence Class #####
plasmid_count <- table(all_plasmids[all_plasmids$Derep, "Class"])
host_count <- table(all_hosts[all_hosts$Derep, "Class"])

# For direct comparison
these <- intersect(c(cris_host$Cell, cas_host$Cell), c(cris_plasmid$Cell, cas_plasmid$Cell))
cris_host <- cris_host[cris_host$Cell %in% these, ]
cris_plasmid <- cris_plasmid[cris_plasmid$Cell %in% these, ]
cas_host <- cas_host[cas_host$Cell %in% these, ]
cas_plasmid <- cas_plasmid[cas_plasmid$Cell %in% these, ]

cas_plasmid <- cas_plasmid[, c("Acc", "Prediction", "Cell")]
cris_plasmid <- cris_plasmid[, c("Acc", "Subtype", "Cell", "Prediction")]
cris_plasmid <- cris_plasmid[is.na(cris_plasmid$Prediction), ]
cris_plasmid[is.na(cris_plasmid$Subtype), "Subtype"] <- "Ambiguous"
cris_plasmid$Prediction <- cris_plasmid$Subtype
cris_and_cas_plasmid <- rbind(cas_plasmid, cris_plasmid[, c("Acc", "Prediction", "Cell")])
cris_and_cas_plasmid <- merge(cris_and_cas_plasmid, tax, by = "Cell")
cris_and_cas_plasmid$Prediction <- ifelse(is.na(cris_and_cas_plasmid$Prediction), "Other", as.character(cris_and_cas_plasmid$Prediction))
cris_and_cas_plasmid$Prediction <- ifelse(grepl("Hybrid|Ambiguous", cris_and_cas_plasmid$Prediction), "Other", cris_and_cas_plasmid$Prediction)
prev_plasmid_ClassP <- aggregate(Acc ~ Class + Prediction, cris_and_cas_plasmid, function(x) length(unique(x)))
prev_plasmid_ClassP <- merge(prev_plasmid_ClassP, data.frame(plasmid_count), by.x = "Class", by.y = "Var1")
prev_plasmid_ClassP$Fraction <- prev_plasmid_ClassP$Acc / prev_plasmid_ClassP$Freq

cas_host <- cas_host[, c("Acc", "Prediction", "Cell")]
cris_host <- cris_host[, c("Acc", "Subtype", "Cell", "Prediction")]
cris_host <- cris_host[is.na(cris_host$Prediction), ]
cris_host[is.na(cris_host$Subtype), "Subtype"] <- "Ambiguous"
cris_host$Prediction <- cris_host$Subtype
cris_and_cas_host <- rbind(cas_host, cris_host[, c("Acc", "Prediction", "Cell")])
cris_and_cas_host <- merge(cris_and_cas_host, tax, by = "Cell")
cris_and_cas_host$Prediction <- ifelse(is.na(cris_and_cas_host$Prediction), "Other", as.character(cris_and_cas_host$Prediction))
cris_and_cas_host$Prediction <- ifelse(grepl("Hybrid|Ambiguous", cris_and_cas_host$Prediction), "Other", cris_and_cas_host$Prediction)
prev_host_ClassP <- aggregate(Acc ~ Class + Prediction, cris_and_cas_host, function(x) length(unique(x)))
prev_host_ClassP <- merge(prev_host_ClassP, data.frame(host_count), by.x = "Class", by.y = "Var1")
prev_host_ClassP$Fraction <- prev_host_ClassP$Acc / prev_host_ClassP$Freq

classes <- names(plasmid_count[plasmid_count >= 10])

# Fig 1 heatmap
mat_plasmid <- reshape2::dcast(Class ~ Prediction, value.var = "Fraction", fill = 0, 
                data = prev_plasmid_ClassP[prev_plasmid_ClassP$Class %in% classes, c(1,2,5)])
mat_host <- reshape2::dcast(Class ~ Prediction, value.var = "Fraction", fill = 0, 
                data = prev_host_ClassP[prev_host_ClassP$Class %in% classes, c(1,2,5)])

prev_plasmid_ClassP <- prev_plasmid_ClassP[prev_plasmid_ClassP$Class %in% classes, c(1,2,5)]
prev_host_ClassP <- prev_host_ClassP[prev_host_ClassP$Class %in% classes, c(1,2,5)]
prev_plasmid_ClassP$Pair <- paste(prev_plasmid_ClassP$Class, prev_plasmid_ClassP$Prediction)
prev_host_ClassP$Pair <- paste(prev_host_ClassP$Class, prev_host_ClassP$Prediction)

prev_all <- merge(prev_plasmid_ClassP, prev_host_ClassP, by = "Pair", all = TRUE)
prev_all[is.na(prev_all$Fraction.x), "Fraction.x"] <- 0
prev_all[is.na(prev_all$Fraction.y), "Fraction.y"] <- 0
prev_all$Diff <- log2((prev_all$Fraction.x) / (prev_all$Fraction.y))
prev_all$Class <- ifelse(is.na(prev_all$Class.x), as.character(prev_all$Class.y), as.character(prev_all$Class.x))
prev_all$Prediction <- ifelse(is.na(prev_all$Prediction.x), prev_all$Prediction.y, prev_all$Prediction.x)

prev_mat <- reshape2::dcast(Class ~ Prediction, value.var = "Diff", data = prev_all)
rownames(prev_mat) <- prev_mat$Class
prev_mat <- prev_mat[, -1]

prev_mat <- prev_mat[, colnames(prev_mat)[c(1:19,21:23,20)]]

min(unlist(prev_mat)[is.finite(unlist(prev_mat))], na.rm = TRUE)
max(unlist(prev_mat)[is.finite(unlist(prev_mat))], na.rm = TRUE)

prev_mat[prev_mat == -Inf] <- -7
prev_mat[prev_mat == Inf] <- 3

prev_mat <- rbind(prev_mat[!rownames(prev_mat) == "Halobacteria", ],
                  prev_mat[rownames(prev_mat) == "Halobacteria", ])

pheatmap::pheatmap(prev_mat,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   color = c("darkblue","blue3","blue","cyan4","cyan3","cyan","grey","grey","pink","darkred"),
                   breaks = seq(-7, 3, 1),
                   na_col = "white",
                   gaps_col = c(8,10,14,19,21,22),
                   width = 8,
                   height = 4.5,
                   filename = "Figures/Fig1_heatmap_direct.pdf")

# Stats
plasmid_subtypes <- reshape2::dcast(Acc ~ Prediction, data = cris_and_cas_plasmid)
host_subtypes <- reshape2::dcast(Acc ~ Prediction, data = cris_and_cas_host)
missing_host <- colnames(plasmid_subtypes)[!colnames(plasmid_subtypes) %in% colnames(host_subtypes)]
missing_plasmid <- colnames(host_subtypes)[!colnames(host_subtypes) %in% colnames(plasmid_subtypes)]
for(i in missing_host){
    host_subtypes[, i] <- 0
}
for(i in missing_plasmid){
    plasmid_subtypes[, i] <- 0
}
plasmid_subtypes <- plasmid_subtypes[, order(colnames(plasmid_subtypes))]
host_subtypes <- host_subtypes[, order(colnames(host_subtypes))]

mat_subtypes <- rbind(plasmid_subtypes, host_subtypes)
group <- c(rep("Plasmid", nrow(plasmid_subtypes)), rep("Chromosome", nrow(host_subtypes)))

library(indicspecies)

set.seed(42)
indval <- multipatt(mat_subtypes[, -1], group, control = how(nperm=9999), func = "IndVal")
summary(indval, alpha = 0.05/(ncol(mat_subtypes)-1))
