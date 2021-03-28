#!/usr/bin/env python
from tools.genefamily import read_gene_families, split_gene_families

split_gene_families("data/testgf.txt", num=3)


#gfs = read_gene_families("data/testgf.txt", "data/testgf.pep", "data/testgf.cds")
#gfs = read_gene_families("data/testgf.txt")

#for gf in gfs:
#    gf.get_id()
    #gf.get_prot()
    #gf.get_cds()
#    gf.obtain_msa()
#    gf.build_phylo()