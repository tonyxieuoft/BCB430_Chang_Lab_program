import os

from Basic_Tools.lists_and_files import file_to_list


def enter_gene_filepath() -> str:

    print("Enter a valid file path containing names and descriptions for GENES "
          "of interest. These will be used as keywords to query the NCBI "
          "database.")

    valid_gene_query_file = False
    gene_query_filepath = ""
    while not valid_gene_query_file:
        print("Filepath:")
        gene_query_filepath = input()
        if not os.path.isfile(gene_query_filepath):
            print("Invalid filepath. Please enter again.")
        else:
            valid_gene_query_file = True
            gene_arr = file_to_list(gene_query_filepath)
            for gene_line in gene_arr:
                for search_query in gene_line.split("\t"):
                    search_components = search_query.split(":")
                    if len(search_components) != 2 or \
                            search_components[0].strip("\"").strip() not in ["g", "d"]:
                        valid_gene_query_file = False

            if not valid_gene_query_file:
                print("Invalid query command(s). Please preface each "
                      "query with either the marker 'g' or 'd' "
                      "indicating the nature of the query as a gene "
                      "name or gene descriptor. Separate markers from "
                      "queries with a colon ':'.")

    return gene_query_filepath


def get_generic_filepath():

    valid_file = False
    filepath = ""
    while not valid_file:
        print("Filepath:")
        filepath = input()

        valid_file = os.path.isfile(filepath)
        if not valid_file:
            print("Invalid filepath. Please enter again.")

    return filepath


def get_generic_directory():

    valid_dir = False
    dir_path = ""
    while not valid_dir:
        print("Directory path:")
        dir_path = input()

        valid_dir = os.path.isdir(dir_path)
        if not valid_dir:
            print("Invalid directory path. Please enter again.")

    return dir_path
