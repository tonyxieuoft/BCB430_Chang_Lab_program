import shutil

from Bio import Entrez
import os

from Bio.Blast import NCBIWWW

from After_BLAST.concatenate_gene_results import concatenate_gene_results
from Basic_Tools.lists_and_files import make_unique_directory, unique_filepath
from Basic_Tools.numeric_user_input import numeric_user_input
from Gene_Description_Refiner.gene_description_refiner import \
    gene_description_refiner
from NCBI_Exon_Puller.handle_ncbi_exon_puller import handle_ncbi_exon_puller
from NCBI_Genome_Blaster.driver_genome_blaster import driver_genome_blaster
from NCBI_Genome_Blaster.driver_genome_blaster_v2 import \
    DriverExonGenomeBlasterV2, DriverFullGenomeBlasterV2, DriverGenomeBlasterV2
from NCBI_Genome_Blaster.local_genome_blaster import local_genome_blaster
from Prepare_For_BLAST.NEPR_convert import convert_NEPR_directory
from Prepare_For_BLAST.prepare_query_files import ExonBlastPreparer, \
    ExonBlastPreparer, FullBlastPreparer, GeMoMaPreparer
from Quality_Checking.get_longest_transcript import \
    optimize_transcripts_by_length
from Quality_Checking.quality_analysis import QualityAnalyser
from Server_Genome_Blaster.automate_server_blast import ServerExonGenomeBlaster, \
    ServerFullGenomeBlaster, ServerImprovedExonGenomeBlaster, ServerImprovedFullGenomeBlaster, GemomaRunner
from User_Interaction.expect_threshold_user_input import \
    expect_threshold_user_input
from User_Interaction.file_acceptors import enter_gene_filepath, \
    get_generic_directory, get_generic_filepath, enter_taxa_filepath


class UI:

    def __init__(self):
        self.email = ""
        self.download_dir = ""

        self.taxa_filepath = ""
        self.gene_query_filepath = ""
        self.refined_gene_query_filepath = ""

        self.exon_pull_dir = ""

        self.blast_reference_dir = ""

        self.blast_results_path = ""

        self.alignments_path = ""

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

            print("Please specify a valid directory (full path) to which the pipeline can "
                  "download files:")
            download_dir = input()

            valid_directory = os.path.isdir(download_dir)
            if not valid_directory:
                print("Invalid directory: " + download_dir)

        self.download_dir = download_dir

    def main_view(self):
        print()
        print("========================= MAIN MENU =============================")
        print("(1) Pull existing sequence data from the NCBI GENE database.")
        print("(2) Create a new BR directory containing the contents of an NEPR directory")
        print("(3) Run NCBI BLAST to pull exons from whole GENOMES.")
        print("(4) Concatenate gene sequences into alignment files.")
        print("(5) Refine gene names and descriptions used to query the NCBI "
              "GENE database.")
        print("(6) Quit the program")
        print("Enter a number from 1-6 to select from one of the above options:")

        return numeric_user_input(1, 6, "")

    def exon_puller_view(self):

        print()
        print("======================== NCBI EXON PULLER =======================")

        discard_gene_query = 1
        if self.refined_gene_query_filepath != "":
            print()
            print("Path of last gene query file used: " + self.refined_gene_query_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_gene_query = numeric_user_input(1, 2, "")
            if discard_gene_query == 1:
                self.gene_query_filepath = self.refined_gene_query_filepath
        elif self.gene_query_filepath != "":
            print()
            print("Path of last gene query file used: " + self.gene_query_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_gene_query = numeric_user_input(1, 2, "")
        if self.gene_query_filepath == "" or discard_gene_query == 2:
            self.gene_query_filepath = enter_gene_filepath()

        discard_taxa_file = 1
        if self.taxa_filepath != "":
            print()
            print("Path of last assigned_taxa file used: " + self.taxa_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_taxa_file = numeric_user_input(1, 2, "")
        if self.taxa_filepath == "" or discard_taxa_file == 2:
            print()
            self.taxa_filepath = enter_taxa_filepath()

        print()
        print("The program can pull the sequences either as separate exons or "
              "full sequences. Enter:")
        print("(1) for exons")
        print("(2) for full sequences")

        exons_or_full_choice = numeric_user_input(1, 2, "")
        if exons_or_full_choice == 1:
            exons_or_full = "exons"
        else:
            exons_or_full = "full"

        self.exon_pull_dir = make_unique_directory(self.download_dir, "NCBI_exon_pull_results")
        print()
        print("The pulled exon results will be at the following path: " +
              self.exon_pull_dir)
        print("Enter any key to continue")
        input()

        print("-----------------------------------------------------------")

        print("Starting NCBI Exon Puller...")
        # ask if by exons or complete cds
        handle_ncbi_exon_puller(self.exon_pull_dir, self.gene_query_filepath,
                                self.taxa_filepath, exons_or_full)
        print("-----------------------------------------------------------")
        print("Finished!")
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

        print()
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

        # TODO Selenium disabled
        gene_description_refiner(self.exon_pull_dir, temp_homology_directory,
                                 self.gene_query_filepath, self.refined_gene_query_filepath)

        shutil.rmtree(temp_homology_directory)

        print("Finished!")
        print("Please DOUBLE CHECK the refined query file before using it. "
              "Sometimes, the refiner is TOO LIBERAL with the addition of new "
              "gene names.")

        print("Enter any key to return to the main menu.")
        input()

    def genome_blaster_view(self):

        print()
        print("====================== PREPARING FOR BLAST ======================")

        if self.taxa_filepath != "":
            print("Please enter a file path to taxa of interest you will like to pull BLAST results for.")
            print("Last taxa file path entered: " + self.taxa_filepath)
            print("Enter:\n"
                  "(1) to use this file\n"
                  "(2) to enter a different file\n"
                  "(3) to return to the main menu")

            taxapath_choice = numeric_user_input(1, 3, "")
            if taxapath_choice == 3:
                return None
            elif taxapath_choice == 2:
                self.taxa_filepath = enter_taxa_filepath()

        else:
            print("Please enter a file path to taxa of interest you will like to pull results for.")
            self.taxa_filepath = enter_taxa_filepath()
            print()

        print()
        print("-------------------------------------------------------------------")

        print("NCBI BLAST requires query sequences to identify "
              "homologous regions in genomic subject sequences. For this "
              "program, the directory containing query sequences must be "
              "in the Blast Reference (BR) format (see readme for details).")
        print()

        if self.blast_reference_dir != "":
            print("Last BR directory created during this session: " + self.blast_reference_dir)
            print("Enter:")
            print("(1) to use this directory")
            print("(2) to use a different directory (must also be in BR format!)")
            print("(3) to exit and return to the main menu")

            blast_reference_choice = numeric_user_input(1, 3, "")
            if blast_reference_choice == 3:
                return None
            elif blast_reference_choice == 2:
                self.blast_reference_dir = get_generic_directory()

        elif self.exon_pull_dir != "":
            print("Last accessed directory created by the NCBI exon puller: " + self.exon_pull_dir)
            print("Enter:")
            print("(1) to use this directory (which will be copied and automatically converted to BR format)")
            print("(2) to use a different directory (must be in BR format!)")
            print("(3) to exit and return to the main menu")

            exon_dir_choice = numeric_user_input(1, 3, "")

            if exon_dir_choice == 3:
                return None
            elif exon_dir_choice == 2:
                self.blast_reference_dir = get_generic_directory()
            elif exon_dir_choice == 1:
                self.blast_reference_dir = convert_NEPR_directory(self.exon_pull_dir, self.download_dir)
                print()
                print("New Blast Reference directory created at: " + self.blast_reference_dir)
                print()
        else:
            print("No previous BR directory creation or NCBI exon puller iteration has occurred"
                  "during this program session.")
            print("Enter:")
            print("(1) to enter a path to a directory containing reference sequences\n"
                  "    (must be in Blast Reference (BR) format! Select option 2 in the main menu"
                  "     if you wish to use results from the program's exon pulling feature)")
            print("(2) to exit and return to the main menu")

            directory_choice = numeric_user_input(1, 2, "")

            if directory_choice == 1:
                self.blast_reference_dir = get_generic_directory()
            else:
                return None

        print()
        print("In order for query sequences to be as similar as possible to "
              "the subject genomes they are blasted against, we must assign "
              "reference species to sub-assigned_taxa that they are closest to within "
              "the overarching taxon of "
              "interest. Details on how the program "
              "automatically assigns references species to sub-assigned_taxa are "
              "available in the README on GitHub. Please enter:")
        print("(1) for automatic assignment (RECOMMENDED)")
        print("(2) for manual assignment")

        auto_assign = numeric_user_input(1, 2, "")
        if auto_assign == 1:
            print("Automatic assignment selected.")
        else:
            print("Manual assignment selected. A prompt to input a path to a "
                  "manual taxon assignment file will appear later.")
        auto_fill = 1  # I don't see how this should ever be manual
        print()
        print("Keep queries as exons, or concatenate into full-length cds "
              "sequences? Enter:")
        print("(1) to keep queries as exons")
        print("(2) to concatenate them into full-length sequences")

        # incorporate this into creating query files
        exon_or_full_query_choice = numeric_user_input(1, 2, "")
        if exon_or_full_query_choice == 1:
            print("Queries will remain as exons.")
        else:
            print("Queries will be concatenated into full-length sequences.")
        print()

        print("Press 1 for normal, or 2 for GeMoMa (all above commands will be overwritten)")
        pick_gemoma = numeric_user_input(1, 2, "")

        queries_path = make_unique_directory(self.download_dir, "query_files")
        print("Great. Query files will be at the path: " + queries_path)
        print("Enter any key to continue.")
        input()

        print("-----------------------------------------------------------------")
        print("Making query files...")

        gffs_path = ""
        if pick_gemoma == 2:

            gffs_path = make_unique_directory(self.download_dir, "gffs")
            print("MAKING A GFF PATH TOO AT: " + gffs_path)
            blast_preparer = GeMoMaPreparer(self.blast_reference_dir, self.taxa_filepath, queries_path, gffs_path)

        else:
            if exon_or_full_query_choice == 1:
                blast_preparer = ExonBlastPreparer(self.blast_reference_dir, self.taxa_filepath, queries_path)
            else:
                blast_preparer = FullBlastPreparer(self.blast_reference_dir, self.taxa_filepath, queries_path)

        blast_preparer.prepare_query_files(auto_assign)
        #print(blast_preparer.get_queries_to_genes_to_exons())


        print("-----------------------------------------------------------------")
        print("Done making query files!")
        print("Enter any key to start configuring settings for BLAST.")
        input()
        print()
        print("========================= GENOME BLASTER =========================")

        print("This program can access BLAST through two different methods. First, "
              "it can emulate a web user using Selenium to access NCBI's server to "
              "BLAST remotely. Secondly (CURRENTLY UNAVAILABLE), it can run BLAST locally; during this "
              "process, the program sequentially downloads genomes then deletes "
              "them immediately after. PLEASE NOTE that local BLAST requires "
              "the NCBI BLAST+ command line application to have been "
              "previously downloaded. Additionally, it must "
              "be run in a linux-style command line. ")
        print("Enter:")
        print("(1) for remote BLAST via Selenium")
        print("(2) for local BLAST")
        print("(3) for server BLAST")
        print("(4) if you have already selected GeMoMa")

        remote_or_local_or_server = numeric_user_input(1, 4, "")

        # throwaway variable for server blast only, might change later
        genome_storage_path = ""
        gemoma_download_dir = ""
        print()
        if remote_or_local_or_server == 1:
            print("Remote BLAST selected.")
        elif remote_or_local_or_server == 2:
            print("Local BLAST selected. Please note that up to 5 GB of space may "
                  "be required for genome download.")
        elif remote_or_local_or_server == 3:
            print("Server BLAST selected.")
            print()
            print("Please enter the directory to which the server's genomes are"
                  " currently stored (at the level containing the species data"
                  " file).")
            genome_storage_path = get_generic_directory()
        else:
            print("GeMoMa confirmed")
            print()
            print("Please enter the directory to which the server's genomes (NON BLAST DB) are "
                  "currently stored, at the level containing the species data file")
            genome_storage_path = get_generic_directory()
            gemoma_download_dir = make_unique_directory(self.download_dir, "gemoma_downloads")
            print("Great, note that all GeMoMa results will be stored here: " + gemoma_download_dir)


        # if we haven't picked gemoma
        expect_value = -1
        if pick_gemoma != 2:
            print("Default expect threshold is 0.05. Enter:\n"
                  "(1) to proceed\n"
                  "(2) to enter a custom expect threshold")
            expect_choice = numeric_user_input(1, 2, "", "Incorrect. Enter 1 or 2.")
            if expect_choice == 1:
                expect_value = 0.05

            if expect_choice == 2:
                expect_value = 0
                print("Enter an expect threshold value (0 < e-value <= 1):")
                expect_value = expect_threshold_user_input()
            print()


        self.blast_results_path = make_unique_directory(self.download_dir, "blast_results")
        print("Blast results will be at the path: " + self.blast_results_path)
        print("Enter any key to begin the blasting process")
        input()
        print("-----------------------------------------------------------------")
        print("Starting BLAST...")

        if remote_or_local_or_server == 1:
            NCBIWWW.email = self.email
            if exon_or_full_query_choice == 1:
                pass
                # TODO Selenium disabled
                genome_blaster = DriverExonGenomeBlasterV2(self.blast_results_path, queries_path, blast_preparer.taxa_blast_order,
                                                           blast_preparer.complete_reference_species)
            else:
                pass
                # TODO Selenium disabled
                genome_blaster = DriverFullGenomeBlasterV2(self.blast_results_path, queries_path, blast_preparer.taxa_blast_order,
                                                           blast_preparer.complete_reference_species, blast_preparer.queries_to_genes_to_exons)
                print(blast_preparer.queries_to_genes_to_exons)
            genome_blaster.blast_genomes(expect_value, self.exon_pull_dir)

        elif remote_or_local_or_server == 2:
            local_genome_blaster(self.blast_results_path, queries_path, blast_preparer.taxa_blast_order,
                                 blast_preparer.complete_reference_species, expect_value)

        elif remote_or_local_or_server == 3:
            if exon_or_full_query_choice == 1:
                #genome_blaster = ServerExonGenomeBlaster(self.blast_results_path, queries_path, blast_preparer.taxa_blast_order,
                #                                         blast_preparer.complete_reference_species, genome_storage_path, self.exon_pull_dir)
                genome_blaster = ServerImprovedExonGenomeBlaster(self.blast_results_path, queries_path,
                                                                 blast_preparer.taxa_blast_order, blast_preparer.complete_reference_species,
                                                                 genome_storage_path, blast_preparer.taxa_to_codes)
            else:
                genome_blaster = ServerImprovedFullGenomeBlaster(self.blast_results_path, queries_path,
                                                                 blast_preparer.taxa_blast_order,
                                                                 blast_preparer.complete_reference_species,
                                                                 genome_storage_path, blast_preparer.taxa_to_codes)
                #genome_blaster = ServerFullGenomeBlaster(self.blast_results_path, queries_path, blast_preparer.taxa_blast_order,
                #                                         blast_preparer.complete_reference_species, genome_storage_path,
                #                                         blast_preparer.taxa_to_codes, blast_preparer.queries_to_genes_to_exons)
            genome_blaster.download_new_genomes()
            genome_blaster.blast_genomes(expect_value)

        else:
            genome_blaster = GemomaRunner(self.blast_results_path, queries_path,
                                          blast_preparer.taxa_blast_order, blast_preparer.complete_reference_species,
                                          genome_storage_path, blast_preparer.taxa_to_codes,
                                          gffs_path, gemoma_download_dir)# last two are the unique ones
            genome_blaster.download_new_genomes()
            genome_blaster.blast_genomes(expect_value)


        print("-----------------------------------------------------------------")
        print("Finished!")
        print("Enter any key to return to the main menu.")
        input()

    def alignment_creator_view(self):

        print("======================== ALIGNMENT CREATOR ======================")
        print("In the future, this program will automatically align sequences. "
              "For now though, its only functionality is to combine "
              "sequences from different species for the same gene into one "
              "fasta file. ")
        print()

        dirs_to_align = []
        exist_choice = 0
        none_yet_choice = 0

        if self.exon_pull_dir != "" or self.blast_results_path != "":

            print("Detected paths to directories containing sequences for alignment:")
            if self.exon_pull_dir != "":
                print(self.exon_pull_dir)
            if self.blast_results_path != "":
                print(self.blast_results_path)
            print("Enter:")
            print("(1) to select these directories")
            print("(2) to select other, different directories")
            print("(3) to exit and return to the main menu.")

            exist_choice = numeric_user_input(1, 3, "")

            if exist_choice == 1:
                if self.exon_pull_dir != "":
                    dirs_to_align.append(self.exon_pull_dir)
                if self.blast_results_path != "":
                    dirs_to_align.append(self.blast_results_path)
        else:
            print("No directories containing sequences for alignment have been "
                  "created during this program session.")
            print("Enter:")
            print("(1) to manually select directories")
            print("(2) to exit and return to the main menu")

            none_yet_choice = numeric_user_input(1, 2, "")

        if exist_choice == 3 or none_yet_choice == 2:
            return None
        elif exist_choice == 2 or none_yet_choice == 1:

            add_more_dirs = True
            while add_more_dirs:

                dirs_to_align.append(get_generic_directory())
                print("Add more directories for alignment? Enter:")
                print("(1) for yes")
                print("(2) for no")
                print("(3) to exit and return to the main menu")

                add_more_choice = numeric_user_input(1, 3, "")

                if add_more_choice == 3:
                    return None
                else:
                    add_more_dirs = (add_more_choice == 1)
        print()
        self.alignments_path = make_unique_directory(self.download_dir, "alignments")
        print("Directories selected! Gene alignments will be at the path: " + self.alignments_path)
        print("Enter any key to continue")
        input()
        print("-----------------------------------------------------------------")
        print("Organizing sequences into alignments...")

        concatenate_gene_results(dirs_to_align, self.alignments_path)

        print("-----------------------------------------------------------------")
        print("Finished!")
        print()
        print("Enter any key to return to the main menu.")
        input()

    def convert_NEPR_to_BR_view(self):

        print("Please enter a NEPR directory, whose contents will be copied and stored in "
              "in a BR directory. (Note that NEPR directories created during this program's "
              "instance will automatically be converted.")
        print()
        to_convert = get_generic_directory()
        self.blast_reference_dir = convert_NEPR_directory(to_convert, self.download_dir)

        print()
        print("Finished. The new directory is at: " + self.blast_reference_dir)

        print()

        print("Enter any key to return to the main menu.")
        input()
















