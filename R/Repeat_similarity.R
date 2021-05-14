# Load
repeat_sim <- read.table("../Collect/Repeat_sim.csv", sep=",", header=FALSE, stringsAsFactors = FALSE)
colnames(repeat_sim) <- c("First", "Second", "Sim")

# Derep
load("Prepared.RData")

library(ggplot2)

p <- ggplot(repeat_sim, aes(Sim)) +
    theme_bw() +
    geom_density(color = "red") +
    geom_vline(xintercept = 0.85) +
    xlab("Sequence identity") +
    ylab("Density") +
    scale_x_continuous(labels = scales::percent)
p

ggsave(p, file = "Figures/Repeat_sim.pdf", units = "cm", width = 8, height = 5)
