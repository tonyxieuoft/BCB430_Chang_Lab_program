import os
from Bio import Entrez


def missing_elasmobranch_genes(general_folder: str):

    for gene in os.listdir(general_folder):

        gene_path = os.path.join(general_folder, gene)
        elasmobranch_path = os.path.join(gene_path, "elasmobranchii")
        if len(os.listdir(elasmobranch_path)) == 0:
            print(gene)


if __name__ == "__main__":

    path = r"C:\Users\tonyx\Downloads\NCBI_exons_filtered"
    sum = 0
    less = 0
    none = 0
    elasmobranch_sum = 0
    genes = 0
    for gene in os.listdir(path):
        gene_path = os.path.join(path, gene)

        genes += 1

        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)
            sum += len(os.listdir(taxon_path))

            if taxon == "elasmobranchii":
                elasmobranch_sum += len(os.listdir(taxon_path))

            if len(os.listdir(taxon_path)) == 0:
                print(gene + " " + taxon)
                none += 1
            if len(os.listdir(taxon_path)) < 5:
                #print(gene + " " + taxon)
                less += 1


    print("sum: " + str(sum))
    print("none: " + str(none))
    print("less: " + str(less))
    print("elasmo_sum: " + str(elasmobranch_sum))
    print("genes: " + str(genes))


    #Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    #gene_table_file = Entrez.efetch(db="gene", id="118884076", rettype="gene_table", retmode="text")

    #print(str(gene_table_file.read()))
    # missing_elasmobranch_genes(r'C:\Users\tonyx\Downloads\NCBI_exon_pull_results (10)')
