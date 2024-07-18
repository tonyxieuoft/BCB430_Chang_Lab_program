import os
import time
from typing import List

from Basic_Tools.lists_and_files import list_to_string
from Basic_Tools.taxonomy_browser import get_taxa_taxids, get_taxonomy_lineage
from NCBI_Genome_Blaster import driver_genome_blaster as v1

from selenium import webdriver
from selenium.common import WebDriverException
from selenium.webdriver import ActionChains, Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from NCBI_Genome_Blaster.assemble_blast_result_sequences import BlastXMLParser
from Basic_Tools.driver_tools import get_element, get_elements, try_click, \
    try_get

MAX_NUM_PROCESSES = 10


def juggle_blast_tabs(driver, num_tabs, sleep_time) -> int:
    """
    Switches between num_tabs tabs until the feature specified by the parameters
    is clickable. For this program, the only tabs that will be "juggled" are
    BLAST tabs while the results haven't finished loading. It is assumed the
    tabs to be juggled are in positions 1 - num_tabs (the NCBI datasets home
    page is position 0.

    :param driver: driver
    :param by_what: type of the description
    :param description: description of the element to click
    :param num_tabs: number of tabs to juggle between
    :param sleep_time:
    :return: the position of the tab that has loaded
    """

    switch_counter = 1
    driver.switch_to.window(driver.window_handles[1])

    while not try_click(driver, By.XPATH, "//a[text()[. = 'Download All']]") \
            and not try_get(driver, By.CLASS_NAME, "error"):
        switch_counter = (switch_counter % num_tabs) + 1
        driver.switch_to.window(driver.window_handles[switch_counter])

        time.sleep(sleep_time)

    return switch_counter



def driver_genome_blaster_v2(save_path: str, queries_path: str,
                             taxa_blast_order: List[str],
                             complete_reference_species: List[str],
                             expect_value: str, exon_pull_dir: str):

    # chromedriver initiation
    chrome_options = webdriver.ChromeOptions()

    chrome_options.add_argument('--headless=new')
    #chrome_options.add_argument('--disable-gpu')
    driver = webdriver.Chrome(chrome_options)
    driver.execute_cdp_cmd("Page.setDownloadBehavior", {
        "behavior": "allow",
        "downloadPath": save_path
    })


    # blast order from list to dictionary
    blast_order_dict = {}
    for i in range(len(taxa_blast_order)):
        blast_order_dict[taxa_blast_order[i]] = i

    driver.get("https://www.ncbi.nlm.nih.gov/")

    # get taxids of assigned_taxa in query file source directory
    taxids_to_taxa = get_taxa_taxids(exon_pull_dir)
    all_info = []

    for taxon_taxid in taxids_to_taxa:

        print("getting genome accessions for '" + taxids_to_taxa[taxon_taxid] + "' ...")
        driver.execute_script("window.open('https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=" + taxon_taxid + "&reference_only=true')")
        driver.switch_to.window(driver.window_handles[-1])

        loaded_identifier = get_element(driver, By.XPATH, "//div[text()[. = 'Rows per page']]", 40)
        table_size_text = get_element(driver, By.XPATH, "//span[@data-testid='table-size']", 40).text
        table_size = int(table_size_text.split()[0])

        if table_size > 20:
            v1.clicking_wrapper(driver, driver, By.CSS_SELECTOR,
                             ".MuiSelect-select.MuiSelect-outlined.MuiInputBase-"
                             "input.MuiOutlinedInput-input.css-zcubqt", 40)
            # selects the maximum number of results that can be displayed (100)
            v1.clicking_wrapper(driver, driver, By.XPATH, "//li[@data-value='100']", 40)

            loaded_identifier = get_element(driver, By.XPATH, "//div[text()[. = 'Rows per page']]", 40)

        species_names = get_elements(driver, By.XPATH, ".//a[@data-ga-label='taxon_name']", 40)
        accessions = get_elements(driver, By.XPATH, "//td[contains(text(), 'GCA_')]", 40)

        for i in range(len(species_names)):
            species = {"name": species_names[i].text,
                       "acc": accessions[i].text,
                       "taxon": taxids_to_taxa[taxon_taxid],
                       "del": None}
            all_info.append(species)

        num_checked = 100
        while num_checked < table_size:
            status = get_element(driver, By.XPATH, "//div[@data-testid='paging-panel__status']", 40).text
            v1.clicking_wrapper(driver, driver, By.XPATH, "//button[@data-ga-label='next_page']", 40)

            while get_element(driver, By.XPATH, "//div[@data-testid='paging-panel__status']", 40).text == status:
                time.sleep(0.01)

            species_names = get_elements(driver, By.XPATH, ".//a[@data-ga-label='taxon_name']", 40)
            accessions = get_elements(driver, By.XPATH, "//td[contains(text(), 'GCA_')]", 40)
            for i in range(len(species_names)):
                species = {"name": species_names[i].text,
                           "acc": accessions[i].text,
                           "taxon": taxids_to_taxa[taxon_taxid],
                           "del": None}
                all_info.append(species)

            num_checked += 100

        driver.close()

    species_so_far = {}
    for species in complete_reference_species:
        species_so_far[species] = True

    novel_species_info = []
    for i in range(len(all_info)):
        species = all_info[i]
        if species["name"] not in species_so_far:
            novel_species_info.append(species)
            species_so_far[species["name"]] = True

    available_species_string = ""
    if len(novel_species_info) > 0:
        available_species_string = novel_species_info[0]['name']
    for i in range(1, len(novel_species_info)):
        available_species_string += "\n" + novel_species_info[i]['name']

    print("contacting NCBI regarding taxonomy info...")
    all_lineage_info = get_taxonomy_lineage(available_species_string)

    # this actually almost takes instant time

    for species in novel_species_info:
        lineage = all_lineage_info[species["name"]]
        curr_del = None
        curr_rank = len(taxa_blast_order)
        for i in range(len(lineage)):
            if lineage[i] in blast_order_dict:
                if blast_order_dict[lineage[i]] < curr_rank:
                    curr_del = lineage[i]
                    curr_rank = blast_order_dict[lineage[i]]

        species["del"] = curr_del

    filter_out_no_del = []
    for species in novel_species_info:
        if species['del'] is not None:
            filter_out_no_del.append(species)

    novel_species_info = filter_out_no_del

    print("blasting processes have begun...")

    finished_jobs = 0
    num_processes = 0
    s_num = 0
    species_open = []
    while finished_jobs < len(novel_species_info):

        if num_processes < MAX_NUM_PROCESSES and not s_num == len(novel_species_info):

            curr_species = novel_species_info[s_num]

            driver.switch_to.window(driver.window_handles[-1])
            driver.execute_script("window.open('https://blast.ncbi.nlm.nih.gov/" +
                                  "Blast.cgi?PAGE_TYPE=BlastSearch&" +
                                  "PROG_DEF=blastn&BLAST_SPEC=GDH_" +
                                  curr_species['acc'] + "')")
            driver.switch_to.window(driver.window_handles[-1])

            reference_filepath = os.path.join(queries_path, curr_species['del'] + ".fas")
            v1.clicking_wrapper(driver, driver, By.XPATH,
                             "//*[text()[. = 'Somewhat similar sequences (blastn)']]",
                             40)
            file_input = get_element(driver, By.ID, "upl", 40)
            file_input.send_keys(reference_filepath)

            if expect_value != 0:
                v1.configure_expect_threshold(driver, expect_value)

            v1.clicking_wrapper(driver, driver, By.CLASS_NAME,
                                "blastbutton", 40)

            species_open.append(curr_species)
            num_processes += 1
            s_num += 1

        if num_processes == MAX_NUM_PROCESSES or s_num == len(novel_species_info):

            tab_no = juggle_blast_tabs(driver, num_processes, 1)

            curr_species = species_open[tab_no-1]
            if not try_get(driver, By.CLASS_NAME, "error"):

                # download the xml file for the blast results, get its path
                v1.xml_download_clicker(driver)
                file_to_analyze = v1.get_downloaded_xml_file(save_path)

                # parse the xml file, creating files that contain the
                # results in fasta format
                while True:
                    try:
                        parser = BlastXMLParser(file_to_analyze, save_path, curr_species['taxon'],
                                                curr_species['name'])
                        parser.parse_blast_xml()
                        break
                    except PermissionError:
                        # print("problem with permissions")
                        time.sleep(1)

                os.remove(file_to_analyze)

                print("(" + str(finished_jobs + 1) + "/" + str(len(novel_species_info))
                      + ")" + " finished blasting " + curr_species['name'])
            else:
                print("(" + str(finished_jobs + 1) + "/" + str(len(novel_species_info))
                      + ")" + " unexpected error " + curr_species['name'])

            driver.close()
            num_processes -= 1
            finished_jobs += 1
            species_open.pop(tab_no - 1)

# class = "error"


if __name__ == "__main__":
    save_path = r"C:\Users\tonyx\Downloads"
    queries_path = ""
    taxa_blast_order = []
    complete_reference_species = []
    expect_value = 0
    exon_pull_dir = r"C:\Users\tonyx\Downloads\NCBI_exons_bat"
    driver_genome_blaster_v2(save_path, queries_path, taxa_blast_order,
                             complete_reference_species, expect_value,
                             exon_pull_dir)


    # xiaohan.xie@mail.utoronto.ca
    # C:\Users\tonyx\Downloads\main_test
    # C:\Users\tonyx\Downloads\main_test\NCBI_exon_pull_results



