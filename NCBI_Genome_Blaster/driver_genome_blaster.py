import os.path
import time
import os
from typing import List

from selenium import webdriver
from selenium.common import WebDriverException
from selenium.webdriver import ActionChains, Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from NCBI_Genome_Blaster.assemble_blast_result_sequences import parse_blast_xml
from Basic_Tools.driver_tools import get_element, get_elements, try_click

WAIT_CONSTANT = 0.1
MAX_NUM_TABS = 2


class DriverTimeoutException(Exception):
    pass


def clicking_wrapper(overall_driver, specific_driver, by_what, description,
                     timer_limit) -> None:
    """
    Clicks an element on the NCBI datasets page. Keeps attempting until
    timer_limit is reached or success occurs. On guard for a survey pop-up,
    which the program immediately closes upon detection. Note that this
    functions very similarly to selenium's built in WebDriverWait class, but
    has additional functionality tailored to NCBI

    :param overall_driver: general driver used
    :param specific_driver: a more specific driver that limits the scope of
    elements considered
    :param by_what: the type of the description (ex. By.ID, By.XPATH)
    :param description: descriptor of the element to click, as a string
    :param timer_limit: stop attempting once this much time has passed.
    :return:
    """

    counter = 0
    while counter < timer_limit:
        try:
            # finds the element and attempts to click
            specific_driver.find_element(by_what, description).click()
            break
        except WebDriverException:
            try:
                # checks if survey has popped up, not allowing clicking
                overall_driver.find_element(By.XPATH,
                                            "//button[text()[. = 'No Thanks']]").click()
            except WebDriverException:
                pass

            # pauses for a period of time
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    if counter >= timer_limit:
        # throw an exception if time limit is reached
        raise DriverTimeoutException(
            "Timeout error for clicking wrapper, description:" +
            description)


def search_button_clicker(driver, timer_limit) -> None:
    """
    The button to click after entering a taxon to search is fickle sometimes
    in that the click goes through, yet no change to the website occurs. This
    function repeatedly clicks the search button until the next page is loaded.

    :param driver: driver
    :param timer_limit: time to keep attempting until an exception is thrown
    :return: none
    """

    counter = 0
    while counter < timer_limit:
        try:
            driver.find_element(By.ID, "panel-header")
            break
        except WebDriverException:
            try_click(driver, By.XPATH, "//button[text()[. = 'Search']]")
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    if counter >= timer_limit:
        raise DriverTimeoutException("search button clicking timeout")


def blast_menu_button_clicker(driver, timer_limit) -> bool:
    """
    Clicker function specifically for the menu button "BLAST against this
    genome. The reason why it requires its own function, rather than just the
    general clicker, is because sometimes the button will not appear as an
    option. TO verify that the menu has loaded but the button is not present,
    checks to see if the "view Details" button has loaded

    :param driver: driver
    :param timer_limit: time to keep attempting until an exception is thrown
    :return: True iff "Blast against this genome" is successfully clicked
    """

    counter = 0
    while counter < timer_limit:
        try:
            # tries clicking the "BLAST against this genome" button
            clicking_wrapper(driver, driver, By.XPATH,
                             "//li[text()[. = 'BLAST against this genome']]", 4)
            return True
        except DriverTimeoutException:
            # if the View Details button is present, means that the "BLAST
            # against this genome" button is not available
            if EC.element_to_be_clickable(
                    (By.XPATH, "//li[text()[. = 'View details']]")):
                return False

        time.sleep(WAIT_CONSTANT)
        counter += WAIT_CONSTANT

    # throws exception if timer runs out
    raise DriverTimeoutException("blast menu button clicking timeout")


def xml_download_clicker(driver) -> None:
    """
    Clicker specifically to download the results of a BLAST file. After clicking
    "Download All", the "XML" button can be finicky to click, so the menu is
    looked over and over again until success is reached.

    :param driver: driver
    :return: None
    """

    time.sleep(1)
    while not try_click(driver, By.XPATH, "//a[text()[. = 'XML']]"):
        time.sleep(1)
        try_click(driver, By.XPATH,
                  "//a[text()[. = 'Download All']]")
        time.sleep(1)


def get_downloaded_xml_file(path: str) -> str:
    """
    Waits until the XML file is downloaded then returns its file path
    :param path: the path of the directory that will store the blast results
    (the XML file will be temporary and discorded soon after)
    :return: the file path of the XML file
    """

    downloaded = False
    while not downloaded:
        for filepath in os.listdir(path):
            if not os.path.isdir(filepath) and os.path.splitext(filepath)[
                1] == ".xml":
                return os.path.join(path, filepath)


def juggle_blast_tabs(driver, by_what, description, num_tabs,
                      sleep_time) -> int:
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

    while not try_click(driver, by_what, description):
        switch_counter = (switch_counter % num_tabs) + 1
        driver.switch_to.window(driver.window_handles[switch_counter])

        time.sleep(sleep_time)

    return switch_counter


def enter_taxa(driver: webdriver, taxa: str) -> None:
    """
    Enter taxa into the NCBI search bar.

    :param driver: driver
    :param taxa: taxa to enter into the search bar
    :return: None
    """
    # get the search bar element, click on it, and delete anything in it
    search_bar = get_element(driver, By.ID, "taxonomy_autocomplete", 40)
    clicking_wrapper(driver, driver, By.ID, "taxonomy_autocomplete", 40)
    search_bar.send_keys(Keys.BACKSPACE)
    # sometimes the page changes, need to get the element again
    search_bar = get_element(driver, By.ID, "taxonomy_autocomplete", 40)
    # search with the subject_taxa (name of the query file)
    search_bar.send_keys(taxa)
    search_bar.send_keys(Keys.RETURN)
    # click the search button to commence search
    search_button_clicker(driver, 40)


def configure_genomes_display(driver: webdriver) -> None:
    """
    Configure the NCBI datasets genome display to display only reference
    genomes, and display a maxmimum of 100 results

    :param driver: driver
    :return: None
    """
    # clicks the filter menu
    clicking_wrapper(driver, driver, By.ID, "panel-header", 40)
    # selects for only reference genomes to be displayed
    clicking_wrapper(driver, driver, By.XPATH, "//*[@value='reference_only']",
                     40)
    # waits for the reference genome to be displayed
    time.sleep(2)

    # clicks the "number of results displayed" menu button
    clicking_wrapper(driver, driver, By.CSS_SELECTOR,
                     ".MuiSelect-select.MuiSelect-outlined.MuiInputBase-"
                     "input.MuiOutlinedInput-input.css-zcubqt", 40)
    # selects the maximum number of results that can be displayed (100)
    clicking_wrapper(driver, driver, By.XPATH, "//li[@data-value='100']", 40)
    time.sleep(2)


def configure_expect_threshold(driver: webdriver, expect: str) -> None:
    """
    On the BLAST settings page, configure expect threshold

    :param driver: driver
    :param expect: user-entered expect threshold
    :return: None
    """

    clicking_wrapper(driver, driver, By.ID, "btnDescrOver", 40)
    expect_button = get_element(driver, By.ID, "expect", 40)
    for i in range(5):
        expect_button.send_keys(Keys.BACKSPACE)
    expect_button.send_keys(str(expect))


def driver_genome_blaster(save_path: str, queries_path: str,
                          taxa_blast_order: List[str],
                          complete_reference_species: List[str], expect_value):
    """
    Automated BLAST using a selenium driver. Pop-up chrome tabs will appear
    and files will be downloaded, but only the complete blast results will
    remain at the end execution of the program.

    :param save_path: the path that any files downloaded by the driver will be
    saved to.
    :param queries_path: directory containing query files generated by the
    Prepare_For_Blast package.
    :param taxa_blast_order: the order in which the genomes of taxa will be
    searched up in the NCBI Datasets database and blasted.
    :return:
    """

    # set the webdriver's download path, allow automatic downloads and safe
    # browsing
    chrome_options = webdriver.ChromeOptions()
    prefs = {'download.default_directory': save_path,
             'profile.default_content_setting_values.automatic_downloads': 1,
             "safebrowsing.enabled": True}
    chrome_options.add_experimental_option('prefs', prefs)
    driver = webdriver.Chrome(chrome_options)

    # access the gene search page, which will have the refseq data
    driver.get("https://www.ncbi.nlm.nih.gov/datasets/genome/")

    # a dictionary that acts as a hash table to check if species have been
    # blasted before
    species_so_far = {}
    # these are for species we already have references for
    for species in complete_reference_species:
        species_so_far[species] = True
    # the order of the subject_taxa to be searched is NON-trivial, since once a
    # species has been BLASTED, it will not be BLASTED again, even if a
    # different query sequence is being used.
    for taxa in taxa_blast_order:

        # reference query path
        reference_filepath = os.path.join(queries_path, taxa + ".fas")
        # search and enter taxa
        enter_taxa(driver, taxa)
        # configure settings for genome display
        configure_genomes_display(driver)
        # gets the total number of rows to iterate over
        genome_table_rows = get_elements(driver, By.CSS_SELECTOR,
                                         ".MuiTableRow-root.css-1mc2v0b", 40)

        tabs_open = 0
        # tracks which tabs correspond to which species (not directly evident on
        # the blast page
        species_open = []

        for i in range(1, len(genome_table_rows)):

            # need to get the genome table rows element
            genome_table_rows = get_elements(driver, By.CSS_SELECTOR,
                                             ".MuiTableRow-root.css-1mc2v0b",
                                             40)
            # the current species in the genome display we are on
            current_row = genome_table_rows[i]

            taxon_name = get_element(current_row, By.XPATH,
                                     ".//a[@data-ga-label='taxon_name']",
                                     40).text

            if taxon_name in species_so_far.keys():
                pass
                # DO NOT BLAST AGAIN if a species has already been blasted
            else:
                # open the menu
                print("on species: " + taxon_name)
                species_so_far[taxon_name] = True

                clicking_wrapper(driver, current_row, By.XPATH,
                                 ".//button[@data-ga-action='click_open_menu']",
                                 40)

                # if the species is blastable
                if blast_menu_button_clicker(driver, 40):

                    # a new BLAST tab opens -> enter necessary information
                    tabs_open += 1
                    species_open.append(taxon_name)
                    driver.switch_to.window(driver.window_handles[tabs_open])

                    # blast page options
                    clicking_wrapper(driver, driver, By.XPATH,
                                     "//*[text()[. = 'Somewhat similar sequences (blastn)']]",
                                     40)
                    file_input = get_element(driver, By.ID, "upl", 40)
                    file_input.send_keys(reference_filepath)

                    if expect_value != 0:
                        configure_expect_threshold(driver, expect_value)

                    clicking_wrapper(driver, driver, By.CLASS_NAME,
                                     "blastbutton", 40)

                # if the maximum number of blast entries has been reached
                if tabs_open == MAX_NUM_TABS:

                    # wait until the results of a tab show up
                    tab_no = juggle_blast_tabs(driver, By.XPATH,
                                               "//a[text()[. = 'Download All']]",
                                               MAX_NUM_TABS, 2)
                    print("tab no: " + str(tab_no))

                    # download the xml file for the blast results, get its path
                    xml_download_clicker(driver)
                    file_to_analyze = get_downloaded_xml_file(save_path)

                    # parse the xml file, creating files that contain the
                    # results in fasta format
                    parse_blast_xml(file_to_analyze, save_path, taxa,
                                    species_open[tab_no - 1])
                    os.remove(file_to_analyze)

                    driver.close()
                    tabs_open -= 1
                    species_open.pop(tab_no - 1)
                    print(species_open)

                driver.switch_to.window(driver.window_handles[0])
                ActionChains(driver).move_by_offset(10, 10).click().perform()

        # after blast searches have been run on every available species in the
        # subject_taxa, winds down and collects results from the remaining tabs that
        # are running
        for i in range(tabs_open, 0, -1):

            tab_no = juggle_blast_tabs(driver, By.XPATH,
                                       "//a[text()[. = 'Download All']]",
                                       i, 2)

            print("tab no: " + str(tab_no))

            xml_download_clicker(driver)
            file_to_analyze = get_downloaded_xml_file(save_path)

            parse_blast_xml(file_to_analyze, save_path, taxa,
                            species_open[tab_no - 1])
            os.remove(file_to_analyze)

            species_open.pop(tab_no - 1)

            print(species_open)
            driver.close()

        driver.switch_to.window(driver.window_handles[0])
        time.sleep(0.5)


if __name__ == "__main__":

    blast_results_path = r'C:\Users\tonyx\Downloads\blast_test_again7'
    os.mkdir(blast_results_path)
    queries_path = r'C:\Users\tonyx\Downloads\query_files (4)'
    taxa_blast_order = ['9738', '9753', '9741', '9732', '9740', '40150', '9756',
                        '119500', '90247', '27609', '9750', '9729', '30558',
                        '9726', '9722', '2746895', '9771', '9721', '378069',
                        '30483', '259919', '117887', '1158979', '117861',
                        '117851', '170819', '30503', '119203', '7778']
    reference_species = ['Balaenoptera acutorostrata', 'Balaenoptera musculus',
                         'Balaenoptera ricei', 'Delphinapterus leucas',
                         'Delphinus delphis', 'Eubalaena glacialis',
                         'Globicephala melas', 'Kogia breviceps',
                         'Lagenorhynchus albirostris',
                         'Lagenorhynchus obliquidens', 'Lipotes vexillifer',
                         'Mesoplodon densirostris', 'Monodon monoceros',
                         'Neophocaena asiaeorientalis asiaeorientalis',
                         'Orcinus orca', 'Phocoena sinus', 'Physeter catodon',
                         'Tursiops truncatus', 'Amblyraja radiata',
                         'Carcharodon carcharias', 'Chiloscyllium plagiosum',
                         'Hemiscyllium ocellatum', 'Hypanus sabinus',
                         'Leucoraja erinacea', 'Mobula hypostoma',
                         'Pristis pectinata', 'Rhincodon typus',
                         'Scyliorhinus canicula', 'Stegostoma tigrinum']

    driver_genome_blaster(blast_results_path, queries_path, taxa_blast_order,
                          reference_species, 0.1)
    pass
