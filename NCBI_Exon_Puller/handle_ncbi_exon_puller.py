import os
from typing import List, Tuple
from Bio import Entrez

from Basic_Tools.lists_and_files import file_to_list

from NCBI_Exon_Puller.ncbi_exon_puller_v2 import ncbi_exon_puller_v2

NCBI_CALL_ATTEMPTS = 5


class UserInputException(Exception):
    pass


def parse_gene_input_line(gene_line: str) -> Tuple[List[str], List[str]]:
    """
    Parses a single line of the user-inputted gene search queries text file
    :param gene_line: the user-inputted line to parse
    :return: a list of gene name queries and a list of gene description queries
    """

    gene_queries = []
    description_queries = []

    # each query is seperated by commas
    for search_parameter in gene_line.split("\t"):
        # each query is either g or d, then a colon, the the query keywords
        search_specs = search_parameter.split(":")

        # gene name
        if search_specs[0].strip("\"").strip() == "g":
            gene_queries.append(search_specs[1].strip("\"").strip().upper())
        # gene description
        elif search_specs[0].strip("\"").strip() == "d":
            description_queries.append(search_specs[1].strip("\"").strip().upper())

    return gene_queries, description_queries


def build_search_query(gene_queries: List[str], description_queries: List[str],
                       taxon: str) -> str:
    """
    Builds the complete search query to input into NCBI for one gene and one
    taxon.

    :param gene_queries: a list of gene name queries (official gene name + other
    acronyms it may go by)
    :param description_queries: a list of gene description queries (product
    name, etc)
    :param taxon: the taxon to search for
    :return: the correctly formatted search query in string form
    """
    divider = ""
    search_query = "("
    # appends the gene name queries to the search
    for gene_query in gene_queries:
        search_query += divider + "\"" + gene_query + "\"[gene]"
        divider = " OR "

    # appends the gene description queries to the search
    for description_query in description_queries:
        search_query += divider + "\"" + description_query + "\"[ALL]"
        divider = " OR "

    search_query += ") AND \"" + taxon + "\"[orgn]"

    return search_query


def handle_ncbi_exon_puller(save_path, genes_filepath, taxon_filepath, exons_or_full):

    # parse user input (every gene line in genes_lines corresponds to a set of
    # queries for a particular gene) - see documentation for more
    genes_lines = file_to_list(genes_filepath)
    taxa = file_to_list(taxon_filepath)

    # follows the format g:gene,d:descriptor,... where each query is separated
    # by commas and specified as a gene name or gene descriptor by "g" or "d",
    # followed by a colon
    for gene_line in genes_lines:
        # obtain gene names and descriptions to search
        gene_queries, description_queries = parse_gene_input_line(gene_line)
        # the 'name' of the gene used to name the folder is the first
        # query in the line
        gene_name = gene_line.split("\t")[0].split(":")[1].strip("\"").strip()

        # create gene folder
        gene_folder = os.path.join(save_path, gene_name)
        os.mkdir(gene_folder)

        for taxon in taxa:

            # create taxon folder
            taxon_folder = os.path.join(gene_folder, taxon)
            os.mkdir(taxon_folder)

            # build the query
            search_query = build_search_query(gene_queries,description_queries, taxon)
            print("searching with query: " + search_query + "...")

            # everything else is done by the use case
            ncbi_exon_puller_v2(search_query, gene_queries, description_queries,
                             taxon_folder, gene_name, exons_or_full)


if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    user_specified_directory = r'C:\Users\tonyx\Downloads\NCBI_length'
    #genes_filepath = r'C:\Users\tonyx\Downloads\gene_queries - Copy.txt'
    genes_filepath = r'C:\Users\tonyx\Downloads\refined_query_file_no_gene_name.txt'
    taxon_filepath = r'C:\Users\tonyx\Downloads\phototransduction_taxa.txt'

    handle_ncbi_exon_puller(user_specified_directory, genes_filepath,
                            taxon_filepath, "exons")
