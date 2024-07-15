from typing import Dict, List

from Bio import Entrez

from Basic_Tools.lists_and_files import list_to_string


def parse_feature_table(ft_table: str, gene_efetch_order: List) -> Dict:

    if not gene_efetch_order:
        return {}

    names_dict = {}
    for gene in gene_efetch_order:
        if gene not in names_dict:
            names_dict[gene] = []

    ft_arr = ft_table.split("\n")
    gene_efetch_counter = 0

    line_counter = 0
    # general loop - checks for "feature fasta heading"
    while line_counter < len(ft_arr) and ft_arr[line_counter].strip():

        line_counter += 1  # this is for the fasta heading for the feature
        curr_gene = gene_efetch_order[gene_efetch_counter]

        keywords = ft_arr[line_counter].split()
        # loops through, stopping either at hitting an empty line or a "gene"
        while len(keywords) != 0 and \
                not (len(keywords) == 3 and
                     keywords[0].isdigit() and
                     keywords[1].isdigit() and
                     keywords[2] == "gene"):

            line_counter += 1
            keywords = ft_arr[line_counter].split()

            # if we've found a "gene" feature (denoted by digits in front first)
        if len(keywords) == 3 and keywords[0].isdigit() and \
                keywords[1].isdigit() and keywords[2] == "gene":

            potential_queries = []
            has_like_in_name = False

            line_counter += 1
            keywords = ft_arr[line_counter].split()

            # iterate through the gene feature itself, looking for gene + desc
            # once we hit the next "number", we are onto the next feature
            while len(keywords) > 0 and not keywords[0].isdigit():
                potential_query = ""
                #if keywords[0] == "gene":
                #    potential_query = "g:" + keywords[1]
                if keywords[0] == "gene_desc":
                    gene_description = list_to_string(keywords[1:], " ")
                    if gene_description.find("-like") == -1:
                        potential_query = "d:" + gene_description
                    else:
                        has_like_in_name = True

                if potential_query and \
                        potential_query not in names_dict[curr_gene] and \
                        potential_query.count(":") == 1:
                    potential_queries.append(potential_query)

                line_counter += 1
                keywords = ft_arr[line_counter].split()

            if not has_like_in_name:
                names_dict[curr_gene] += potential_queries

            # we got what we're looking for, onto the next feature
            while ft_arr[line_counter].split():
                line_counter += 1

        # this is for the space
        line_counter += 1

        gene_efetch_counter += 1

    return names_dict


def get_gene_names_from_accessions(ids, gene_efetch_order):

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    efetch_handle = Entrez.efetch(db='nuccore', id=ids, rettype="ft")

    return parse_feature_table(str(efetch_handle.read()), gene_efetch_order)
