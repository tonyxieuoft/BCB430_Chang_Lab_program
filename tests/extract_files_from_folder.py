import os

if __name__ == "__main__":

    in_path = r"C:\Users\tonyx\Downloads\NCBI_exon_pull_results5"
    out_dir = r"C:\Users\tonyx\Downloads\hypanus"

    for gene in os.listdir(in_path):
        gene_path = os.path.join(in_path, gene)

        for taxa in os.listdir(gene_path):
            taxa_path = os.path.join(gene_path, taxa)

            for species in os.listdir(taxa_path):
                species_path = os.path.join(taxa_path, species)

                for file in os.listdir(species_path):

                    f1 = open(os.path.join(species_path, file), "r")
                    f2 = open(os.path.join(out_dir, gene + "_" + file), "w")

                    f2.write(f1.read())
