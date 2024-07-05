import os
from Bio import Entrez


def missing_elasmobranch_genes(general_folder: str):

    for gene in os.listdir(general_folder):

        gene_path = os.path.join(general_folder, gene)
        elasmobranch_path = os.path.join(gene_path, "elasmobranchii")
        if len(os.listdir(elasmobranch_path)) == 0:
            print(gene)


if __name__ == "__main__":

    path = r"C:\Users\tonyx\Downloads\NCBI_exon_pull_results (10)"
    for gene in os.listdir(path):
        gene_path = os.path.join(path, gene)

        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            if len(os.listdir(taxon_path)) < 1:
                print(gene + " " + taxon)


    #Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    #gene_table_file = Entrez.efetch(db="gene", id="118884076", rettype="gene_table", retmode="text")

    #print(str(gene_table_file.read()))
    # missing_elasmobranch_genes(r'C:\Users\tonyx\Downloads\NCBI_exon_pull_results (10)')
