import shutil

from Bio import Entrez
import os

from Basic_Tools.lists_and_files import make_unique_directory, unique_filepath
from Basic_Tools.numeric_user_input import numeric_user_input
from Gene_Description_Refiner.gene_description_refiner import \
    gene_description_refiner
from NCBI_Exon_Puller.handle_ncbi_exon_puller import handle_ncbi_exon_puller
from Quality_Checking.get_longest_transcript import \
    optimize_transcripts_by_length
from User_Interaction.user_exon_pulling import enter_gene_filepath, \
    get_generic_directory, get_generic_filepath


class UI:

    def __init__(self):
        self.email = ""
        self.download_dir = ""

        self.gene_query_filepath = ""
        self.taxa_filepath = ""

        self.refined_gene_query_filepath = ""

        self.exon_pull_dir = ""

    def get_email(self):
        return self.email

    def email_view(self):

        email_confirm = 2
        email = ""
        while email_confirm == 2:

            print("NCBI Entrez and BLAST require an email to contact in case any "
                  "issues arise. Please enter your email. ")
            email = input()

            print("Inputted email: \"" + email + "\"")
            email_confirm_prompt = "Enter 1 to confirm, or 2 to re-enter your email"
            email_confirm = numeric_user_input(1, 2, email_confirm_prompt)

        self.email = email
        Entrez.email = email

    def base_directory_view(self):

        valid_directory = False
        download_dir = ""
        while not valid_directory:

            print("Please specify a valid directory to which files created by the "
                  "pipeline can be downloaded to:")
            download_dir = input()

            valid_directory = os.path.isdir(download_dir)
            if not valid_directory:
                print("Invalid directory: " + download_dir)

        self.download_dir = download_dir

    def main_view(self):

        print("========================= MAIN MENU =============================")
        print("(1) Pull existing exon data from the NCBI GENE database.")
        print("(2) Refine gene names and descriptions used to query the NCBI "
              "GENE database.")
        print("(3) Run NCBI BLAST to pull exons from whole genomes.")
        print("(4) Concatenate gene sequences into alignment files.")
        print("(5) Quit the program")
        print("Enter a number from 1-5 to select from one of the above options:")

        return numeric_user_input(1, 5, "")

    def exon_puller_view(self):

        print("======================== NCBI EXON PULLER =======================")

        discard_gene_query = 1
        if self.refined_gene_query_filepath != "":
            print("Path of last gene query file used: " + self.refined_gene_query_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_gene_query = numeric_user_input(1, 2, "")
            if discard_gene_query == 1:
                self.gene_query_filepath = self.refined_gene_query_filepath
        elif self.gene_query_filepath != "":
            print("Path of last gene query file used: " + self.gene_query_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_gene_query = numeric_user_input(1, 2, "")
        if self.gene_query_filepath == "" or discard_gene_query == 2:
            self.gene_query_filepath = enter_gene_filepath()

        discard_taxa_file = 1
        if self.taxa_filepath != "":
            print("Path of last taxa file used: " + self.taxa_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_taxa_file = numeric_user_input(1, 2, "")
        if self.taxa_filepath == "" or discard_taxa_file == 2:
            print("Enter a valid file path containing taxa of interest to pull for.")
            self.taxa_filepath = get_generic_filepath()

        self.exon_pull_dir = make_unique_directory(self.download_dir, "NCBI_exon_pull_results")
        print("The pulled exon results will be at the following path: " +
              self.exon_pull_dir)
        print("Enter any button to continue")
        input()

        print("-----------------------------------------------------------")

        print("Starting NCBI Exon Puller...")
        handle_ncbi_exon_puller(self.exon_pull_dir, self.gene_query_filepath,
                                self.taxa_filepath)
        print("Finished!")
        print("-----------------------------------------------------------")
        print("For each given gene and species, the exon puller may have pulled out multiple transcript variants. Enter:")
        print("(1) to keep all transcript variants")
        print("(2) to automatically select for 'optimal' variants by sequence length")
        length_optimization_choice = numeric_user_input(1, 2, "")

        if length_optimization_choice == 2:
            print("Beginning selection...")
            optimize_transcripts_by_length(self.exon_pull_dir)
            print("Optimization complete.")
        else:
            print("Keeping all pulled out sequences...")

        print("Enter any button to return to the main menu.")
        input()

    def gene_description_refiner_view(self):

        print("=================== GENE DESCRIPTION REFINER ====================")

        print("The Gene Description Refiner identifies alternative "
              "descriptions for genes of interest within the NCBI database and "
              "adds them to a gene query file that can be used as input for a "
              "subsequent iteration of the NCBI exon "
              "puller. "
              "Via sequence homology searches within the NCBI 'nuccore' "
              "database, it identifies descriptions of relevant sequences not "
              "picked up by a previous run of the NCBI exon puller. A "
              "subsequent run of the NCBI exon puller will utilize these newly "
              "identified descriptions to pull out more exons.")
        print("-----------------------------------------------------------------")
        print("TO OPERATE, the refiner REQUIRES a previous iteration of the "
              "NCBI exon puller to have ran. Previously pulled out exons act as "
              "queries for sequence homology searches. A previous gene query "
              "file must also exist for the refiner to function.")
        if self.exon_pull_dir != "":
            print("Last accessed directory created by the NCBI exon puller: " + self.exon_pull_dir)
            print("Enter:")
            print("(1) to use this directory")
            print("(2) to return to the main menu")

            directory_choice = numeric_user_input(1, 2, "")

            if directory_choice == 2:
                return None
        else:
            print("No previous iteration of the NCBI exon puller has been ran "
                  "during this program session")
            print("Enter:")
            print("(1) to enter a path to a directory previously created by the "
                  "NCBI exon puller.")
            print("(2) to return to the main menu")

            directory_choice = numeric_user_input(1, 2, "")

            if directory_choice == 1:
                self.exon_pull_dir = get_generic_directory()
            else:
                return None

        if self.gene_query_filepath != "":
            print("Last accessed gene query file: " + self.gene_query_filepath)
            print("Enter:")
            print("(1) to use this file")
            print("(2) to return to the main menu")

            filepath_choice = numeric_user_input(1, 2, "")

            if filepath_choice == 2:
                return None
        else:
            print("No previous gene query file has been accessed during this "
                  "program session.")
            print("Enter:")
            print("(1) to enter a path to a gene query file.")
            print("(2) to return to the main menu")

            filepath_choice = numeric_user_input(1, 2, "")

            if filepath_choice == 1:
                self.gene_query_filepath = get_generic_filepath()
            else:
                return None

        self.refined_gene_query_filepath = unique_filepath(os.path.join(self.download_dir, "refined_query_file.txt"))
        print("The refined query file will be at the following path: " + self.refined_gene_query_filepath)

        temp_homology_directory = make_unique_directory(self.download_dir, "homology_search")
        print("Temporary query files will be stored at the directory: " + temp_homology_directory)
        print("Enter any key to begin the refining process.")
        input()

        gene_description_refiner(self.exon_pull_dir, temp_homology_directory,
                                 self.gene_query_filepath, self.refined_gene_query_filepath)

        shutil.rmtree(temp_homology_directory)

        print("Finished!")
        print("Please double check the refined query file before using it. "
              "Sometimes, the refiner is too liberal with the addition of new "
              "gene names.")

        print("Enter any key to return to the main menu.")
        input()














