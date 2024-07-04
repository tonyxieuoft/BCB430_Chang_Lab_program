import os
import time
from typing import Dict, List

from selenium.common import WebDriverException
from selenium.webdriver import Keys

from After_BLAST.concatenate_exons import concatenate_exons
from selenium import webdriver

from selenium.webdriver.common.by import By
from Basic_Tools.driver_tools import get_element, DriverTimeoutException, WAIT_CONSTANT
from Basic_Tools.basic_dictionaries import print_dict
from Basic_Tools.lists_and_files import file_to_list, list_to_string

from Bio import Entrez


LESS_CUTOFF = 5


def get_desc_table(driver: webdriver, timer_limit: int):

    counter = 0
    while counter < timer_limit:
        try:
            # tries finding the description table
            return driver.find_element(By.ID, "dscTable"), True
        except WebDriverException:
            # if not present, looks to see if error message has popped up
            try:
                return driver.find_element(By.ID, "noResMsg"), False
            except WebDriverException:
                # if not, keep waiting
                time.sleep(WAIT_CONSTANT)
                counter += WAIT_CONSTANT

    if counter >= timer_limit:
        DriverTimeoutException("timeout error")


def parse_feature_table(ft_table: str, gene_efetch_order: List) -> Dict:

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

            line_counter += 1
            keywords = ft_arr[line_counter].split()

            # iterate through the gene feature itself, looking for gene + desc
            while len(keywords) > 0 and not keywords[0].isdigit():
                potential_query = ""
                if keywords[0] == "gene":
                    potential_query = "g:" + keywords[1]
                elif keywords[0] == "gene_desc":
                    potential_query = "d:" + list_to_string(keywords[1:], " ")

                if potential_query and \
                        potential_query not in names_dict[curr_gene] and \
                        potential_query.count(":") == 1:
                    names_dict[curr_gene].append(potential_query)

                line_counter += 1
                keywords = ft_arr[line_counter].split()

            # we got what we're looking for, onto the next feature
            while ft_arr[line_counter].split():
                line_counter += 1

        # this is for the space
        line_counter += 1

        gene_efetch_counter += 1

    return names_dict


def make_refined_query_file(names_dict: Dict, original_file: str,
                            save_path: str) -> None:

    new_query_path = os.path.join(save_path, "refined_query_file.txt")
    new_query_file = open(new_query_path, "w")

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


def homology_search(search_requests: List[Dict],
                    queries_path: str) -> Dict:

    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    driver = webdriver.Chrome(chrome_options)

    gene_blast_order = []
    # access the gene search page, which will have the refseq data
    # rate limiting factor
    driver.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome")
    for request in search_requests:
        for key in request.keys():
            if key == "None" or key == "Less":
                for taxon in request[key]:
                    driver.find_element(By.XPATH, "//*[text()[. = 'Somewhat similar sequences (blastn)']]").click()

                    if key == "Less":
                        driver.find_element(By.NAME, "QUERYFILE").\
                            send_keys(os.path.join(queries_path, request["Gene"], taxon + ".fas"))
                    else:
                        driver.find_element(By.NAME, "QUERYFILE"). \
                            send_keys(os.path.join(queries_path, request["Gene"], "none.fas"))
                    #  time.sleep(0.5)
                    while True:

                        driver.find_element(By.NAME, "EQ_MENU"). \
                            click()
                        driver.find_element(By.NAME, "EQ_MENU"). \
                            clear()
                        driver.find_element(By.NAME, "EQ_MENU"). \
                            send_keys(taxon)
                        get_element(driver, By.CLASS_NAME, "ui-ncbiautocomplete-options", 40)
                        driver.find_element(By.XPATH, "//input[@value='BLAST']").\
                            click()
                        try:
                            alert = driver.switch_to.alert
                            print(alert.text)
                            alert.accept()
                            print("alert came up")
                            time.sleep(0.5)
                        except Exception:
                            break

                    print("got to the end!")
                    gene_blast_order.append([request["Gene"], key])
                    driver.execute_script("window.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome')")

                    driver.switch_to.window(driver.window_handles[-1])

    entries = ""
    gene_efetch_order = []
    num_blast_entries = len(driver.window_handles) - 1

    for i in range(num_blast_entries):

        driver.switch_to.window(driver.window_handles[i])
        desc_table, success = get_desc_table(driver, 40)

        if success:

            desc_table_rows = desc_table.find_elements(By.XPATH, ".//tr")[1:]
            for desc_table_row in desc_table_rows:

                if gene_blast_order[i][1] == "Less":
                    cutoff = 0.0
                else:
                    cutoff = 1e-50

                if float(get_element(desc_table_row, By.CLASS_NAME, "c9", 40).text) <= cutoff:
                    entries += desc_table_row.find_element(By.CSS_SELECTOR, ".c12.l.lim").\
                                   find_element(By.XPATH, ".//a").text + ","
                    gene_efetch_order.append(gene_blast_order[i][0])

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    efetch_handle = Entrez.efetch(db='nuccore', id=entries, rettype="ft")

    return parse_feature_table(str(efetch_handle.read()), gene_efetch_order)

    # now just write a function to pull from the gene table
    # efetch results are in order of the ids inputted in -> easiy match to gene
    # invariant is separation via >, then, a series of two numbers -> category,
    # then category description -> look for "gene" and take from "gene",
    # "gene desc"


def gene_description_refiner(exon_pull_path: str, save_dir: str,
                             original_query_file: str):

    search_requests = []

    queries_path = os.path.join(save_dir, "homology_search")
    os.mkdir(queries_path)

    for gene in os.listdir(exon_pull_path):

        query_path = os.path.join(save_dir, "homology_search", gene)
        os.mkdir(query_path)

        gene_path = os.path.join(exon_pull_path, gene)

        gene_info = {"Gene": gene, "None": [], "Less": [], "Good": []}
        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            available_species = os.listdir(taxon_path)
            if len(available_species) == 0:
                gene_info["None"].append(taxon)
            elif len(available_species) < LESS_CUTOFF:
                gene_info["Less"].append(taxon)
            else:
                gene_info["Good"].append(taxon)

        if len(gene_info["Less"]) + len(gene_info["Good"]) == 0:
            print("Gene: " + gene + " seems to have no ref. seq.")
        else:

            for less_taxon in gene_info["Less"]:
                query_file = open(
                    os.path.join(query_path, less_taxon + ".fas"), "w"
                )
                query_file.write(get_random_query_sequence(
                    os.path.join(gene_path, less_taxon)
                ))
                query_file.close()

            if len(gene_info["None"]) > 0:
                query_file = open(
                    os.path.join(query_path, "none.fas"), "w"
                )
                if len(gene_info["Good"]) != 0:
                    query_file.write(get_random_query_sequence(
                        os.path.join(gene_path, gene_info["Good"][0])
                    ))
                else:
                    query_file.write(get_random_query_sequence(
                        os.path.join(gene_path, gene_info["Less"][0])
                    ))
                query_file.close()

            search_requests.append(gene_info)

    new_queries = homology_search(search_requests, queries_path)
    make_refined_query_file(new_queries, original_query_file, save_dir)

    os.remove(queries_path)


if __name__ == "__main__":
    # exon_path = r"C:\Users\tonyx\Downloads\refine_test"
    exon_path = r"C:\Users\tonyx\Downloads\of_interest\NCBI_exon_pull_results (11)"
    og_query_file = r"C:\Users\tonyx\Downloads\gene_queries - Copy.txt"
    gene_description_refiner(exon_path, r"C:\Users\tonyx\Downloads", og_query_file)
