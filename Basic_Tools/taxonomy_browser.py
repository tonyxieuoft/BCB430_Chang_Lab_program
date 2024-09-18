import os
import time

from Bio import Entrez
#from selenium import webdriver
#from selenium.webdriver.common.by import By
from typing import Dict

PHYLOGENY_COLUMN = 3
SPECIES_NAME_COLUMN = 1
CODE_COLUMN = 0

'''
def get_taxonomy_lineage_driver(species_string: str) -> Dict:
    """
    Given a string of species, returns complete lineages for each species as
    values in a dictionary, where the names of the species are the keys.

    :param species_string: a string of species names, separated by newlines
    :return: a dictionary, where keys are species names and values are lists of
    lineages.
    """
    # assemble a new driver with the "headless" option turned on
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    driver = webdriver.Chrome(chrome_options)

    # access the gene search page, which will have the refseq data
    # rate limiting factor
    driver.get("https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi")

    # text box for entering species' names
    driver.find_element(By.NAME, "tax").send_keys(species_string)

    # option for complete phylogeny
    driver.find_element(By.NAME, "lng").click()

    driver.find_element(By.XPATH, "//input[@value='Show on screen']").click()
    # table of complete phylogenies
    table = driver.find_element(By.XPATH, "//h2")
    # each row is for a separate species
    table_elements = table.find_elements(By.XPATH, ".//td")
    phylo_dict = {}

    i = 4
    while i < len(table_elements):
        phylo_dict[table_elements[i+SPECIES_NAME_COLUMN].text] = \
            table_elements[i+PHYLOGENY_COLUMN].text.split(" ")
        i += 4

    return phylo_dict
'''

def get_taxa_taxids(general_dir):

    # TODO remove the entrez email
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    genes = os.listdir(general_dir)
    taxa = os.listdir(os.path.join(general_dir, genes[0]))
    taxids = {}

    for taxon in taxa:
        handle = Entrez.esearch(db="taxonomy", term=taxon)
        taxid = Entrez.read(handle)['IdList'][0]
        taxids[taxid] = taxon

    return taxids


def get_taxonomy_lineage(species_string: str):
    """
    Given a string of species, returns complete lineages for each species as
    values in a dictionary, where the names of the species are the keys.

    :param species_string: a string of species names, separated by newlines
    :return: a dictionary, where keys are species names and values are lists of
    lineages.
    """
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    ids = ""
    all_species = species_string.split("\n")
    for species in all_species:
        handle = Entrez.esearch(db="taxonomy", term=species)
        time.sleep(0.02)
        taxid = Entrez.read(handle)['IdList'][0]
        ids += taxid + ","

    phylo_dict = {}

    handle2 = Entrez.efetch(db="taxonomy", id=ids)
    results = Entrez.read(handle2)
    for org in results:
        name = org["ScientificName"]
        phylo_dict[name] = [org['TaxId']] # add the organism's own first
        lineage = org["LineageEx"]
        for i in range(len(lineage)-1, -1, -1):
            phylo_dict[name].append(lineage[i]["TaxId"])

    return phylo_dict


if __name__ == "__main__":
    #general_dir = r"C:\Users\tonyx\Downloads\NCBI_exons_bat"
    #print(get_taxa_taxids(general_dir))
    species_string = "Tursiops truncatus\nOrcinus orca\n" \
                     "Balaenoptera acutorostrata\nBalaenoptera musculus\n" \
                     "Balaenoptera ricei\nDelphinapterus leucas\n" \
                     "Delphinus delphis\nEubalaena glacialis\n" \
                     "Globicephala melas\nKogia breviceps\n" \
                     "Lagenorhynchus albirostris\nLagenorhynchus obliquidens\n" \
                     "Mesoplodon densirostris\nMonodon monoceros\n" \
                     "Neophocaena asiaeorientalis asiaeorientalis\n" \
                     "Phocoena phocoena\nPhocoena sinus\n" \
                     "Physeter catodon"
    print(get_taxonomy_lineage(species_string))
