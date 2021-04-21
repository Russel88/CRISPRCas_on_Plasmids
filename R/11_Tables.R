load("Prepared.RData")

write.table(cas_plasmid, file = "Tables/Cas_plasmid.tab", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(cas_host, file = "Tables/Cas_host.tab", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(cris_plasmid, file = "Tables/Crispr_plasmid.tab", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(cris_host, file = "Tables/Crispr_host.tab", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(all_plasmids, file = "Tables/Plasmids.tab", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(all_hosts, file = "Tables/Hosts.tab", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(cell, file = "Tables/Plasmid_Host_Connection.tab", quote = FALSE, row.names = FALSE, sep = "\t")

