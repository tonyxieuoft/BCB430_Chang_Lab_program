import os
import time
from typing import Dict, List

from selenium.common import WebDriverException
from selenium import webdriver

from selenium.webdriver.common.by import By
from Basic_Tools.driver_tools import DriverTimeoutException, WAIT_CONSTANT
from Basic_Tools.lists_and_files import list_to_string

MAX_PROCESSES = 15


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


def enter_blast_process(driver, gene, taxon, queries_path, deficiency_type):

    driver.switch_to.window(driver.window_handles[-1])
    driver.execute_script("window.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome')")
    driver.switch_to.window(driver.window_handles[-1])

    driver.find_element(By.XPATH, "//*[text()[. = 'Somewhat similar sequences (blastn)']]").click()

    if deficiency_type == "Less":
        driver.find_element(By.NAME, "QUERYFILE"). \
            send_keys(os.path.join(queries_path, gene, taxon + ".fas"))
    else:
        driver.find_element(By.NAME, "QUERYFILE"). \
            send_keys(os.path.join(queries_path, gene, "none.fas"))

    driver.find_element(By.NAME, "EQ_TEXT").clear()
    driver.find_element(By.NAME, "EQ_TEXT"). \
        send_keys(taxon + "[orgn]")
    driver.find_element(By.XPATH, "//input[@value='BLAST']"). \
        click()


def scrape_description_tables(driver: webdriver, gene_blast_tab_order: List[str],
                              gene_efetch_order: List[str], entries: List[str],
                              total_processes: int):

    for tab_no in range(1, total_processes + 1):

        driver.switch_to.window(driver.window_handles[1])

        desc_table, success = get_desc_table(driver, 40)

        if success:

            # first "c9" is actually the E-value heading, need to remove
            e_values = desc_table.find_elements(By.CLASS_NAME, "c9")[1:]
            accessions = desc_table.find_elements(By.CSS_SELECTOR, ".c12.l.lim")

            for j in range(len(e_values)):

                if gene_blast_tab_order[tab_no][1] == "Less":
                    cutoff = 0.0
                else:
                    cutoff = 1e-50

                if float(e_values[j].text) <= cutoff:
                    entries.append(accessions[j].text)
                    gene_efetch_order.append(gene_blast_tab_order[tab_no][0])

        driver.close()

    return entries


def homology_search(search_requests: List[Dict],
                    queries_path: str) -> [str, List[str]]:

    chrome_options = webdriver.ChromeOptions()
    #chrome_options.add_argument('--headless')
    #chrome_options.add_argument('--disable-gpu')
    driver = webdriver.Chrome(chrome_options)

    driver.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome')

    gene_blast_tab_order = ["blank"]

    entries = []
    gene_efetch_order = []
    # access the gene search page, which will have the refseq data
    # rate limiting factor
    for request in search_requests:
        for key in request.keys():
            if key == "None" or key == "Less":
                for taxon in request[key]:

                    enter_blast_process(driver, request["Gene"], taxon, queries_path, key)
                    gene_blast_tab_order.append([request["Gene"], key])

                    print("opening tab")

                    if len(driver.window_handles) == MAX_PROCESSES + 1:
                        scrape_description_tables(driver, gene_blast_tab_order,
                                                  gene_efetch_order, entries,
                                                  MAX_PROCESSES)

                        gene_blast_tab_order = ["blank"]

    scrape_description_tables(driver, gene_blast_tab_order, gene_efetch_order,
                              entries, len(driver.window_handles) - 1)

    return list_to_string(entries, ","), gene_efetch_order
