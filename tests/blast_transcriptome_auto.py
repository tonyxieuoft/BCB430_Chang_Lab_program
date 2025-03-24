import os

if __name__ == "__main__":

    reference_cds_path = input("Enter reference full cds file")
    trinity_path = input("Enter path to trinity")

    output_dir = input("Enter path to a new directory (not yet created to store the output")

    for gene in os.listdir(reference_cds_path):
        gene_path = os.path.join(reference_cds_path, gene)

        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            for species in os.listdir(taxon_path):
                species_path = os.path.join(taxon_path, species)

                for transcript in os.listdir(species):

                    os.system()


