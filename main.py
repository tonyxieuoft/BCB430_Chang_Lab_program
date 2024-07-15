import os
import time

from Bio import Entrez
from Bio.Blast import NCBIWWW

from After_BLAST.concatenate_gene_results import concatenate_gene_results
from Basic_Tools.lists_and_files import make_unique_directory
from Basic_Tools.numeric_user_input import numeric_user_input

# Press the green button in the gutter to run the script.
from NCBI_Exon_Puller.handle_ncbi_exon_puller import handle_ncbi_exon_puller
from NCBI_Genome_Blaster.local_genome_blaster import local_genome_blaster
from Prepare_For_BLAST.prepare_query_files import prepare_query_files
from User_Interaction.expect_threshold_user_input import \
    expect_threshold_user_input
from User_Interaction.user_exon_pulling import enter_gene_filepath, \
    get_generic_directory, get_generic_filepath

if __name__ == '__main__':

    print("Welcome to the pipeline...")

    email_confirm = 2
    email = ""
    while email_confirm == 2:

        print("NCBI Entrez and BLAST require an email to contact in case any "
              "issues arise. Please enter your email. ")
        email = input()

        print("Inputted email: \"" + email + "\"")
        email_confirm_prompt = "Enter 1 to confirm, or 2 to re-enter your email"
        email_confirm = numeric_user_input(1, 2, email_confirm_prompt)

    print("=======================================================")

    Entrez.email = email

    valid_directory = False
    download_dir = ""
    while not valid_directory:

        print("Please specify a valid directory to which files created by the "
              "pipeline can be downloaded to:")
        download_dir = input()

        valid_directory = os.path.isdir(download_dir)
        if not valid_directory:
            print("Invalid directory: " + download_dir)

    print("=======================================================")

    print("Choose from the following options by entering the corresponding "
          "number: \n"
          "(1) BLAST immediately\n"
          "(2) Pull exons from NCBI for reference species before starting BLAST")

    skip_pull_step = numeric_user_input(1, 2, "")
    print("=======================================================")

    if skip_pull_step == 1:
        print("skipping exon pulling...")
        print("Enter a path to a directory containing reference sequences "
              "previously pulled out using the NCBI exon puller.")
        exon_pull_path = get_generic_directory()
    else:
        # print statements for gene in the function itself
        gene_query_filepath = enter_gene_filepath()
        print("Enter a valid file path containing taxa of interest to pull for.")
        taxa_filepath = get_generic_filepath()

        exon_pull_path = make_unique_directory(download_dir,
                                               "NCBI_exon_pull_results")
        print("Ok, we are set to pull. The results will be at the following path: " +
              exon_pull_path)
        print("Enter any button to continue")
        input()

        print("=======================================================")

        print("Starting NCBI Exon Puller...")
        handle_ncbi_exon_puller(exon_pull_path, gene_query_filepath,
                                taxa_filepath)
        print("Finished pulling exons!")

        termination_message = "Enter: \n" \
                              "(1) to terminate the program here\n" \
                              "(2) to prepare for blast"
        if numeric_user_input(1, 2, termination_message) == 1:
            print("exiting...")
            exit(0)
    print("============================================")
    print("Now, we'll prepare query files for blasting.")
    print("============================================")
    assign_message = "We assign reference species to sub-taxa of the " \
                     "overarching taxon to be blasted to make sure query " \
                     "sequences are as similar to the subject genome as " \
                     "possible. Enter: \n" \
                     "(1) for automatic assignment\n" \
                     "(2) for manual assignment\n"
    auto_assign = numeric_user_input(1, 2, assign_message)
    print("configuring settings...")
    print("============================================")
    fill_in_message = "Sometimes, a reference species will not have sequences " \
                      "for certain genes. In cases where this occurs, the " \
                      "sequences from the closest available species will be " \
                      "used instead. Enter: \n" \
                      "(1) for this fill-in process to be automatic. \n" \
                      "(2) for manual fill-in."
    auto_fill = numeric_user_input(1, 2, fill_in_message)
    print("configuring settings...")
    queries_path = make_unique_directory(download_dir, "query_files")
    print("Great. Query files will be at the path: " + queries_path)
    time.sleep(0.5)
    print("=============================================")
    print("Making query files...")

    taxa_blast_order, complete_reference_species = \
        prepare_query_files(auto_assign, auto_fill, exon_pull_path, queries_path)

    print("Done making query files!")
    print("=============================================")
    print("Configuring BLAST...")

    print("This program can access BLAST through two different methods. First, "
          "it can emulate a web user using Selenium to access NCBI's server to "
          "BLAST remotely. Secondly, it can run BLAST locally; during this "
          "process, the program sequentially downloads genomes then deletes "
          "them immediately after. Enter:\n"
          "(1) for remote BLAST via Selenium\n"
          "(2) for local BLAST")

    remote_or_local = numeric_user_input(1, 2, "", "Incorrect. Enter 1 or 2.")

    if remote_or_local == 1:
        print("Remote BLAST selected. Please note that a pop-up window will "
              "appear during the blasting process. Do not be alarmed, the "
              "program wil be accessing BLAST as if it is a web user.")
    else:
        print("Local BLAST selected. Please note that up to 5 GB of space may "
              "be required for genome download.")

    print("=============================================")

    print("Default expect threshold is 0.05. Enter:\n"
          "(1) to proceed\n"
          "(2) to enter a custom expect threshold")
    expect_choice = numeric_user_input(1, 2, "", "Incorrect. Enter 1 or 2.")

    expect_value = 0
    if expect_choice == 2:
        print("Enter an expect threshold value. It must be greater than 0 and "
              "less than or equal to 1:")
        expect_value = expect_threshold_user_input()

    print("=============================================")

    blast_results_path = make_unique_directory(download_dir, "blast_results")
    print("Blast results will be at the path: " + blast_results_path)
    print("Enter any key to begin the blasting process")
    input()

    if remote_or_local == 1:
        #genome_blaster = driver_genome_blaster
        NCBIWWW.email = email
    else:
        #genome_blaster = local_genome_blaster
        pass

    local_genome_blaster(blast_results_path, queries_path, taxa_blast_order,
                         complete_reference_species, expect_value)

    print("BLAST complete.")
    print("=============================================")
    alignments_path = make_unique_directory(download_dir, "alignments")
    print("Organizing BLAST results into alignments, will be at the path: " +
          alignments_path)
    print("...")

    concatenate_gene_results([exon_pull_path, blast_results_path], alignments_path)

    print("=============================================")
    print("Program finished.")










