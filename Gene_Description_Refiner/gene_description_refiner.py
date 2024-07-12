import os
from typing import Dict

from After_BLAST2.concatenate_exons import concatenate_exons
from Basic_Tools.lists_and_files import file_to_list, list_to_string
from Gene_Description_Refiner.get_gene_names_from_accessions import \
    get_gene_names_from_accessions
from Gene_Description_Refiner.homology_search_for_accessions import homology_search

LESS_CUTOFF = 5


def make_refined_gene_name_file(names_dict: Dict, original_file: str,
                                refined_path: str) -> None:

    new_query_file = open(refined_path, "w")

    original_arr = file_to_list(original_file)
    for line in original_arr:
        gene_name = line.split("\t")[0].split(":")[1].upper()
        existing_queries = line.split("\t")

        if gene_name in names_dict:
            for new_query in names_dict[gene_name]:
                if new_query not in existing_queries:
                    existing_queries.append(new_query)

        new_query_file.write(list_to_string(existing_queries, "\t") + "\n")

    new_query_file.close()


def get_random_query_sequence(taxon_path: str) -> str or None:

    random_species = os.listdir(taxon_path)[0]
    random_species_path = os.path.join(taxon_path, random_species)
    random_file = os.listdir(random_species_path)[0]
    random_file_path = os.path.join(random_species_path, random_file)

    return concatenate_exons(random_file_path)


def get_genes_needing_refinement(exon_pull_path, temp_homology_directory):

    search_requests = []

    for gene in os.listdir(exon_pull_path):

        # make a directory that will contain query sequences for the gene
        query_path = os.path.join(temp_homology_directory, gene)
        os.mkdir(query_path)

        gene_path = os.path.join(exon_pull_path, gene)

        # keeps track of the lacking taxa (different taxa should have different
        # query sequences)
        gene_info = {"Gene": gene, "None": [], "Less": [], "Good": []}

        # go through each taxon
        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            # get number of reference sequences for the taxon
            available_species = os.listdir(taxon_path)
            if len(available_species) == 0:
                gene_info["None"].append(taxon)
            elif len(available_species) < LESS_CUTOFF:
                gene_info["Less"].append(taxon)
            else:
                gene_info["Good"].append(taxon)

        if len(gene_info["Less"]) + len(gene_info["Good"]) == 0:
            # nothing is pulled out for the gene
            print("Gene: " + gene + " seems to have no ref. seq.")

        else:
            # choose a random query sequence to use and put it in a file
            for less_taxon in gene_info["Less"]:

                # query sequence is from the taxon (since it has at least one)
                new_filename = os.path.join(query_path, less_taxon + ".fas")
                query_file = open(new_filename, "w")
                query_seq = get_random_query_sequence(os.path.join(gene_path, less_taxon))
                query_file.write(query_seq)
                query_file.close()

            if len(gene_info["None"]) > 0:

                query_file = open(os.path.join(query_path, "none.fas"), "w")

                # if there's a taxa with a lot of reference sequences
                # just pull from available pool
                if len(gene_info["Good"]) > 0:
                    query_file.write(get_random_query_sequence(
                        os.path.join(gene_path, gene_info["Good"][0])
                    ))
                else:
                    query_file.write(get_random_query_sequence(
                        os.path.join(gene_path, gene_info["Less"][0])
                    ))
                query_file.close()

            search_requests.append(gene_info)

    return search_requests


def gene_description_refiner(exon_pull_path: str, temp_homology_directory: str,
                             original_query_file: str, refined_query_filename: str):

    # holds query sequences that will be used for the homology search

    # starting with an NCBI exon pull results folder, obtain genes that are
    # missing reference sequences
    search_requests = get_genes_needing_refinement(exon_pull_path, temp_homology_directory)

    # conduct parallel BLAST searches, identify any accessions that might have
    # gene names of interest
    ids_of_interest, gene_efetch_order = homology_search(search_requests, temp_homology_directory)

    # the query sequences are no longer needed
    # shutil.rmtree(queries_path)

    # efetch all ids of interest, parse the feature table
    new_names = get_gene_names_from_accessions(ids_of_interest, gene_efetch_order)

    # based on the pulled out gene names, make a new gene name query file
    make_refined_gene_name_file(new_names, original_query_file, refined_query_filename)


if __name__ == "__main__":
    # exon_path = r"C:\Users\tonyx\Downloads\refine_test"
    exon_path = r"C:\Users\tonyx\Downloads\NCBI_exon_pull_results_non_refined"
    og_query_file = r"C:\Users\tonyx\Downloads\gene_queries - Copy.txt"
    gene_description_refiner(exon_path, r"C:\Users\tonyx\Downloads", og_query_file)
