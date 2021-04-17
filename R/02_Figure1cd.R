library(ggplot2)
library(ape)

# Load data
load("Prepared.RData")

##### Prevalence Class #####
plasmid_count <- table(all_plasmids[all_plasmids$Derep, "Class"])
host_count <- table(all_hosts[all_hosts$Derep, "Class"])

cris_and_cas_plasmid <- rbind(cas_plasmid[, c("Acc", "Prediction", "Cell")], cris_plasmid[, c("Acc", "Prediction", "Cell")])
cris_and_cas_plasmid <- merge(cris_and_cas_plasmid, tax, by = "Cell")
cris_and_cas_plasmid$Prediction <- ifelse(is.na(cris_and_cas_plasmid$Prediction), "Other", as.character(cris_and_cas_plasmid$Prediction))
cris_and_cas_plasmid$Prediction <- ifelse(grepl("Hybrid|Ambiguous", cris_and_cas_plasmid$Prediction), "Other", cris_and_cas_plasmid$Prediction)
prev_plasmid_Class <- aggregate(Acc ~ Class, cris_and_cas_plasmid, function(x) length(unique(x)))
prev_plasmid_Class <- merge(prev_plasmid_Class, data.frame(plasmid_count), by.x = "Class", by.y = "Var1")
prev_plasmid_Class$Fraction <- prev_plasmid_Class$Acc / prev_plasmid_Class$Freq
prev_plasmid_ClassP <- aggregate(Acc ~ Class + Prediction, cris_and_cas_plasmid, function(x) length(unique(x)))
prev_plasmid_ClassP <- merge(prev_plasmid_ClassP, data.frame(plasmid_count), by.x = "Class", by.y = "Var1")
prev_plasmid_ClassP$Fraction <- prev_plasmid_ClassP$Acc / prev_plasmid_ClassP$Freq

cris_and_cas_host <- rbind(cas_host[, c("Acc", "Prediction", "Cell")], cris_host[, c("Acc", "Prediction", "Cell")])
cris_and_cas_host <- merge(cris_and_cas_host, tax, by = "Cell")
cris_and_cas_host$Prediction <- ifelse(is.na(cris_and_cas_host$Prediction), "Other", as.character(cris_and_cas_host$Prediction))
cris_and_cas_host$Prediction <- ifelse(grepl("Hybrid|Ambiguous", cris_and_cas_host$Prediction), "Other", cris_and_cas_host$Prediction)
prev_host_Class <- aggregate(Acc ~ Class, cris_and_cas_host, function(x) length(unique(x)))
prev_host_Class <- merge(prev_host_Class, data.frame(host_count), by.x = "Class", by.y = "Var1")
prev_host_Class$Fraction <- prev_host_Class$Acc / prev_host_Class$Freq
prev_host_ClassP <- aggregate(Acc ~ Class + Prediction, cris_and_cas_host, function(x) length(unique(x)))
prev_host_ClassP <- merge(prev_host_ClassP, data.frame(host_count), by.x = "Class", by.y = "Var1")
prev_host_ClassP$Fraction <- prev_host_ClassP$Acc / prev_host_ClassP$Freq

prev_plasmid_Class$Origin <- "Plasmid"
prev_host_Class$Origin <- "Chromosome"

prevClass <- rbind(prev_plasmid_Class, prev_host_Class)

prevClass <- prevClass[prevClass$Class %in% names(plasmid_count[plasmid_count >= 10]), ]

# Agglom tree
coph <- cophenetic(as.phylo(tree))

tax_class <- tax[, c("Cell","Kingdom", "Phylum", "Class")]
tax_class <- tax_class[tax_class$Class %in% prevClass$Class, ]

classes <- unique(prevClass$Class)

dd <- matrix(0, ncol = length(classes), nrow = length(classes))
rownames(dd) <- colnames(dd) <- classes
for(i in seq_along(classes)){
    for(j in seq_along(classes)){
        if(j > i){
            k <- median(coph[rownames(coph) %in% tax_class[tax_class$Class == classes[i], "Cell"],
                             colnames(coph) %in% tax_class[tax_class$Class == classes[j], "Cell"]])
            dd[i, j] <- k
            dd[j, i] <- k
        }
    }
}

tree_class <- nj(dd)
tree_class <- ladderize(root(tree_class, outgroup = c("Halobacteria", "Methanosarcinia")))
write.tree(tree_class, file = "Figures/Fig1c_iTOL_tree.nwk")

# Database abundance
class_abund <- as.data.frame(plasmid_count)
these <- class_abund[class_abund$Freq >= 10, "Var1"]

# Fig 1 tree iTOL files
handle <- file("Figures/Fig1_tree_iTOL_1.txt")
writeLines(c("DATASET_GRADIENT
SEPARATOR SPACE
DATASET_LABEL Database
COLOR #ff0000
STRIP_WIDTH 25
COLOR_MIN #ffffff
COLOR_MAX #ff0000
DATA"), handle)
close(handle)
write.table(class_abund[class_abund$Freq >= 10, ], row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ", 
            append = TRUE, file = "Figures/Fig1_tree_iTOL_1.txt")

handle <- file("Figures/Fig1_tree_iTOL_2.txt")
writeLines(c("DATASET_SIMPLEBAR
SEPARATOR COMMA
DATASET_LABEL,label 1
COLOR,#bb00bb
DATASET_SCALE,0.05,0.1,0.15,0.2
WIDTH,200
MARGIN,20
DATA"), handle)
close(handle)
write.table(prev_plasmid_Class[prev_plasmid_Class$Class %in% these, c(1,4)], row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ",", 
            append = TRUE, "Figures/Fig1_tree_iTOL_2.txt")

handle <- file("Figures/Fig1_tree_iTOL_3.txt")
writeLines(c("DATASET_SIMPLEBAR
SEPARATOR COMMA
DATASET_LABEL,label 1
COLOR,#00bb00
DATASET_SCALE,0.2,0.4,0.6,0.8,1
WIDTH,200
MARGIN,20
DATA"), handle)
close(handle)
write.table(prev_host_Class[prev_host_Class$Class %in% these, c(1,4)], row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ",", 
            append = TRUE, "Figures/Fig1_tree_iTOL_3.txt")

# Fig 1 heatmap
mat_plasmid <- reshape2::dcast(Class ~ Prediction, value.var = "Fraction", fill = 0, 
                data = prev_plasmid_ClassP[prev_plasmid_ClassP$Class %in% prevClass$Class, c(1,2,5)])
mat_host <- reshape2::dcast(Class ~ Prediction, value.var = "Fraction", fill = 0, 
                data = prev_host_ClassP[prev_host_ClassP$Class %in% prevClass$Class, c(1,2,5)])

prev_plasmid_ClassP <- prev_plasmid_ClassP[prev_plasmid_ClassP$Class %in% prevClass$Class, c(1,2,5)]
prev_host_ClassP <- prev_host_ClassP[prev_host_ClassP$Class %in% prevClass$Class, c(1,2,5)]
prev_plasmid_ClassP$Pair <- paste(prev_plasmid_ClassP$Class, prev_plasmid_ClassP$Prediction)
prev_host_ClassP$Pair <- paste(prev_host_ClassP$Class, prev_host_ClassP$Prediction)

prev_all <- merge(prev_plasmid_ClassP, prev_host_ClassP, by = "Pair", all = TRUE)
prev_all[is.na(prev_all$Fraction.x), "Fraction.x"] <- 0
prev_all[is.na(prev_all$Fraction.y), "Fraction.y"] <- 0
prev_all$Diff <- (prev_all$Fraction.x) - (prev_all$Fraction.y)
prev_all$Class <- ifelse(is.na(prev_all$Class.x), as.character(prev_all$Class.y), as.character(prev_all$Class.x))
prev_all$Prediction <- ifelse(is.na(prev_all$Prediction.x), prev_all$Prediction.y, prev_all$Prediction.x)

prev_mat <- reshape2::dcast(Class ~ Prediction, value.var = "Diff", data = prev_all)
rownames(prev_mat) <- prev_mat$Class
prev_mat <- prev_mat[, -1]

prev_mat <- prev_mat[rev(tree_class$tip.label[tree_class$edge[tree_class$edge[,2] <= length(tree_class$tip.label), 2]]),
         colnames(prev_mat)[c(1:20,22:29,21)]]

pheatmap::pheatmap(prev_mat,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   color = c("blue4","blue2","cyan3","cyan","lightblue","red", "red4"),
                   breaks = c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1),
                   na_col = "white",
                   gaps_col = c(8,11,15,20,25,28),
                   filename = "Figures/Fig1_heatmap.pdf",
                   width = 8,
                   height = 4.5)
