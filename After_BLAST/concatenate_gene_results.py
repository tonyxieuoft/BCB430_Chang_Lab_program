"""
THis is a simple script that combines predictions for the same gene into the same file (for alignment purposes).
"""

import os
from typing import List

from After_BLAST.concatenate_exons import concatenate_exons
from Prepare_For_BLAST.get_longest_transcript import get_longest_transcript


def concatenate_gene_results(paths: List[str], save_path):
    """
    Combine predictions for the same gene into a single file. This function requires two inputs:
    1) paths: a number of paths each pointing to BLAST output directory for which genes can be collected
    and combined into single alignment-like files
    2) save_path: Path to which the combined gene result files can be outputted.
    """

    gene_to_encountered_species = {}

    for general_directory in paths:

        for gene in os.listdir(general_directory):
            # create (or append to) the gene file
            gene_save_file = open(os.path.join(save_path, gene) + ".fas", "a")
            gene_path = os.path.join(general_directory, gene)

            # it's different for every gene, no? but we prioritize the
            # exon puller results regardless
            if gene not in gene_to_encountered_species:
                gene_to_encountered_species[gene] = {}

            # this actually also helps negate species duplicates for overlapping
            # assigned_taxa
            encountered_species = gene_to_encountered_species[gene]

            # here, we go into the taxon-level folder
            for taxa in os.listdir(gene_path):
                taxa_path = os.path.join(gene_path, taxa)

                # here, we go into the species-level folder
                for species in os.listdir(taxa_path):
                    species_path = os.path.join(taxa_path, species)

                    # if we haven't already added a gene for the species into the file, we do so now
                    if species not in encountered_species or \
                            gene not in encountered_species[species]:

                        # this updates the hash table according to the gene-species combinations that we have seen
                        if species not in encountered_species:
                            encountered_species[species] = {gene: True}
                        else:
                            encountered_species[species][gene] = True

                        # IF ALL REFERENCE TRANSCRIPTS ARE TO BE INCLUDED

                        for transcript in os.listdir(species_path):
                            to_write = concatenate_exons(os.path.join(species_path, transcript))
                            gene_save_file.write(to_write)

                        # IF ONLY THE LONGEST ONE IS TO BE INCLUDED
                        # file_to_use = get_longest_transcript(species_path)
                        # if file_to_use != "":
                        #    to_write = concatenate_exons(
                        #         os.path.join(species_path, file_to_use))
                        #    gene_save_file.write(to_write)

            gene_save_file.close()


if __name__ == "__main__":

    exon_pull_path1 = r"C:\Users\tonyx\Downloads\NCBI_length - Copy"
    #exon_pull_path2 = r"C:\Users\tonyx\Downloads\blast_results1"
    # blast_path = r"C:\Users\tonyx\Downloads\blast_test_again4"
    save_path = r"C:\Users\tonyx\Downloads\alignments_opt_taxa"
    concatenate_gene_results([exon_pull_path1], save_path)
