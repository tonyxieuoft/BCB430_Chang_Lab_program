import os
from typing import List

from After_BLAST2.concatenate_exons import concatenate_exons
from Prepare_For_BLAST.get_longest_transcript import get_longest_transcript


def concatenate_gene_results(paths: List[str], save_path):

    encountered_species = {}

    for general_directory in paths:

        for gene in os.listdir(general_directory):
            # create (or append to) the gene file
            gene_save_file = open(os.path.join(save_path, gene) + ".fas", "a")
            gene_path = os.path.join(general_directory, gene)

            for taxa in os.listdir(gene_path):
                taxa_path = os.path.join(gene_path, taxa)

                for species in os.listdir(taxa_path):
                    species_path = os.path.join(taxa_path, species)

                    if species not in encountered_species or \
                            gene not in encountered_species[species]:

                        if species not in encountered_species:
                            encountered_species[species] = {gene: True}
                        else:
                            encountered_species[species][gene] = True

                        file_to_use = get_longest_transcript(species_path)
                        if file_to_use != "":
                            to_write = concatenate_exons(
                                os.path.join(species_path, file_to_use))
                            gene_save_file.write(to_write)

            gene_save_file.close()

if __name__ == "__main__":

    exon_pull_path1 = r"C:\Users\tonyx\Downloads\NCBI_length - Copy"
    #exon_pull_path2 = r"C:\Users\tonyx\Downloads\blast_results1"
    # blast_path = r"C:\Users\tonyx\Downloads\blast_test_again4"
    save_path = r"C:\Users\tonyx\Downloads\alignments_opt_taxa"
    concatenate_gene_results([exon_pull_path1], save_path)
