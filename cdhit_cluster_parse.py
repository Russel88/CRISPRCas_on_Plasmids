#!/usr/bin/env python3

# With the .clstr output from CDHIT make a tab-delimited table

import re
import sys
f = open(sys.argv[1], "r")
f1 = f.readlines()
for line in f1:
    if ">Cluster" in line:
        clust = line.replace(" ", "_").replace(">", "")
    else:
        print(clust.replace("\n"," ") + re.sub(".*>", "", re.sub("\.\.\..*", "", line)).replace("\n", ""))
f.close()
