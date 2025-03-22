import os
from Bio import Entrez

from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence

if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    print("getting gene seqnece")
    #arr = ncbi_get_gene_sequence("NC_067412.1", 85934660, 85966394, "1")
    #print(len(arr))

    for i in range(15):
        print("got one")
        arr = ncbi_get_gene_sequence("NC_067412.1", 85934660, 85966394, "1")

    print("done")

    # it's actually very fast so it's not so bad