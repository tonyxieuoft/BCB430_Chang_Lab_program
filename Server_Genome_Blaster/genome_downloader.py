"""
This file contains scripts automating the process of downloading genomes onto a server. This is necessary for a
server-based genome blaster to run.
"""


import os
import subprocess
from typing import Dict

from Basic_Tools.basic_dictionaries import json_to_dict

from Basic_Tools.lists_and_files import list_to_string
from Basic_Tools.taxonomy_browser import get_taxonomy_lineage


######### species_data file constants #########
SPECIES_DATA_FILENAME = "species_data.csv"

NAME_COL = 0
ACC_COL = 1
LINEAGE_COL = 2


def read_species_data(genome_storage_path):
    """
    The species data file contains information about all genomes on a server. Each row represents a genome, and there
    are three columns per row: one for the species name, one for the accession of the genome, and one for the lineage
    of the species.

    We read the species data file to know what genomes we currently have on the server, to allow us to avoid
    downloading duplicates
    """
    # Get the path to the species data file.
    species_data_path = os.path.join(genome_storage_path, SPECIES_DATA_FILENAME)

    # if the species data file does not exist (i.e., the specified path isn't a file), then return nothing (there is
    # nothing to read)
    if not os.path.isfile(species_data_path):
        return None

    # alteratively...
    species_data = []

    # open the species data file for reading
    f = open(species_data_path, "r")

    line = f.readline()

    while line != "":
        # while there are still lines in the species data file to read...

        # split the row by commas (this is how the columns are separated)
        split_line = line.strip().split(",")

        # store the information about the particular genome
        genome = {"name": split_line[NAME_COL].strip(), "acc": split_line[ACC_COL].strip(),
                  "lineage": split_line[LINEAGE_COL].split()}

        # add the extracted data to a list of known genomes on the server
        species_data.append(genome)

        # keep reading
        line = f.readline()

    return species_data


class ServerGenomeDownloader:
    """
    Here, we define a class for downloading genomes. This is directly called by the scripts in the
    "automate_server_blast.py" file.
    """
    def __init__(self, save_path, taxa_list, genome_storage_path):

        """
        Instantiate the class. It requires three things as input:

        1) save_path: A path to generally download files (but not the genomes)
        2) taxa_list: A list of taxa to potentially download genomes for
        3) genome_storage_path: A path to a directory currently storing a number of genomes in blast db format as
        well as a species data file.
        """
        # instantiate the path to save files
        self.save_path = save_path

        # instantiate a list of taxa to download genomes for
        self.taxa_list = taxa_list

        # recall the path to which genomes are currently downloaded on the server
        self.genome_storage_path = genome_storage_path

        # get all genomes currently downloaded on the server (so we don't duplicate download)
        self.existing_accessions = {}
        self.get_existing_accessions()

        self.accessions_to_download = []

    def get_existing_accessions(self):
        """
        Get all genomes currently downloaded onto the server. This requires accessing the species data file.
        """
        # read the species data file to get all downloaded genomes on the server
        species_data = read_species_data(self.genome_storage_path)

        # if no species data file exists
        if species_data is None:
            return None

        # otherwise, record the genomes already downloaded on the server
        for species in species_data:
            self.existing_accessions[species["name"]] = True

    def get_accessions_to_download(self):
        """
        Return all the accessions of interest that are not yet on the server
        """
        return self.accessions_to_download

    def set_accessions_to_download(self):
        """
        Set the accessions needed for download onto the server
        """
        species_so_far = {}

        # for all the species in our taxa of interest...
        for taxa in self.taxa_list:

            # create an temporary output path to redirect NCBI Datasets Command Line output
            temp_summary_path = os.path.join(self.save_path, "temp_summary.txt")

            # run the NCBI Datasets command which returns all available genomes for a particular taxa, and redirect it
            # into temp_summary_path
            os.system("datasets summary genome taxon \"" + taxa + "\" > '" + temp_summary_path + "'")

            # convert the NCBI Datasets summary to a dictionary to parse: we are looking for reference sequences only
            genome_summary_dict = json_to_dict(temp_summary_path)

            # remove the temporary output file
            os.remove(temp_summary_path)

            # for all of the available genomes for the particular taxa returned by NCBI Datasets...
            for genome_record in genome_summary_dict["reports"]:

                # if the genome is a ref_seq genome, then we choose to download it for our BLAST analyses. Otherwise, no
                if "refseq_category" in genome_record["assembly_info"] and \
                    genome_record["organism"]["organism_name"] not in species_so_far:

                    # create a dictionary storing this information
                    blast_organism_dct = {"species": genome_record["organism"]["organism_name"],
                                          "accession": genome_record["accession"]}

                    # add this species to the list that we have encountered so far (so that we do not download the
                    # same genomes, again
                    species_so_far[blast_organism_dct["species"]] = True

                    if blast_organism_dct["species"] not in self.existing_accessions:
                        self.accessions_to_download.append(blast_organism_dct)

    def write_to_species_data_file(self, org):

        """
        After downloading genomes onto the server, it is important to update the species_data_file so that future
        calls to these functions will see that the genomes are downloaded.

        Here is a function that updates the species data file after a genome is downloaded.
        """

        # open (or create a new species data file)
        species_data_path = os.path.join(self.genome_storage_path, "species_data.csv")
        f = open(species_data_path, "a")

        # get the entire taxonomic lineage of the species
        name_to_lineage = get_taxonomy_lineage(org["species"])
        lineage = list_to_string(name_to_lineage[org["species"]], " ")

        # write the species, genomic accession, and taxonomic lineage to the species data file
        f.write(org["species"] + "," + org["accession"] + "," +
                    lineage + "\n")
        f.close()

    def download_genomes(self):

        """
        Download genomes onto the server. Most of the code in this function represent command line calls
        """

        # for each of our genomes to download...
        for org in self.accessions_to_download:

            # download the genome zip file via the NCBI Datasets command line
            os.system(r"datasets download genome accession " +
                      org["accession"])

            # unzip it into a generic, temporary directory called ncbi_dataset
            os.system("unzip ncbi_dataset.zip")

            # remove zip file
            os.system("rm ncbi_dataset.zip")
            # remove the readme too; these files aren't necessary
            os.system("rm README.md")

            #working_path = subprocess.check_output(["pwd"], shell=True). \
            #    decode("utf-8").strip()
            #os.system("datasets rehydrate --directory " + working_path)

            # get genome fasta file path nested within the temp directory
            genome_file_directory = os.path. \
                join("ncbi_dataset", "data", org["accession"])

            if os.path.isdir(genome_file_directory):

                # get the name of the genome fasta file downloaded by the NCBI Datasets command line using the "ls"
                # command here
                genome_filename = subprocess. \
                    check_output(["ls", genome_file_directory]). \
                    decode("utf-8").strip()
                genome_filepath = os.path. \
                    join(genome_file_directory, genome_filename)

                # create local blast database using the genome file
                blast_db_path = os.path.join(self.genome_storage_path, "blast_db")
                if not os.path.isdir(blast_db_path):
                    os.system("mkdir " + blast_db_path)

                # format a good name for the database, so that other users can also read and use it
                species_db = ""
                first = True
                for word in org["species"].split():
                    if first:
                        species_db = word
                        first = False
                    else:
                        species_db += "_" + word

                new_blast_db_name = os.path.join(blast_db_path, species_db)

                # we make a blast database out of the genomic fasta file that we downloaded
                os.system("makeblastdb -dbtype nucl -in " + genome_filepath +
                          " -out " + new_blast_db_name)

                self.write_to_species_data_file(org)

            # remove original genome file (we only keep the blast database
            os.system("rm -r ncbi_dataset")

            # remove md5sum file as well
            os.system("rm md5sum.txt")



if __name__ == "__main__":


    save_path = "/crun/tony.xie/Downloads"
    taxa_list = ["Pelodiscus sinensis",
                 "Pelusios castaneus",
                 "Platysternon megacephalum",
                 "Podocnemis expansa",
                 "Rafetus swinhoei",
                 "Sternotherus odoratus",
                 "Terrapene triunguis",
                 "Testudo graeca",
                 "Trachemys scripta elegans"]
    genome_storage = "/crun/tony.xie/GenomeStorage"

    downloader = ServerGenomeDownloader(save_path, taxa_list, genome_storage)
    downloader.set_accessions_to_download()
    #downloader.write_to_species_data_file()
    downloader.download_genomes()

