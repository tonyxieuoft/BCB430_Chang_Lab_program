import os
import subprocess

from Basic_Tools.basic_dictionaries import json_to_dict

from Basic_Tools.lists_and_files import list_to_string
from Basic_Tools.taxonomy_browser import get_taxonomy_lineage


######### species_data file constants #########
SPECIES_DATA_FILENAME = "species_data.csv"

NAME_COL = 0
ACC_COL = 1
LINEAGE_COL = 2




def read_species_data(genome_storage_path):

    species_data_path = os.path.join(genome_storage_path, SPECIES_DATA_FILENAME)

    if not os.path.isfile(species_data_path):
        return None

    species_data = []
    f = open(species_data_path, "r")

    line = f.readline()

    while line != "":

        split_line = line.strip().split(",")
        genome = {"name": split_line[NAME_COL], "acc": split_line[ACC_COL],
                  "lineage": split_line[LINEAGE_COL].split()}

        species_data.append(genome)
        line = f.readline()

    return species_data


class ServerGenomeDownloader:

    def __init__(self, save_path, taxa_list, genome_storage_path):

        self.save_path = save_path
        self.taxa_list = taxa_list
        self.genome_storage_path = genome_storage_path

        self.existing_accessions = {}
        self.get_existing_accessions()

        self.accessions_to_download = []

    def get_existing_accessions(self):

        species_data = read_species_data(self.genome_storage_path)
        if species_data is None:
            return None

        for species in species_data:
            self.existing_accessions[species["name"]] = True

    def get_accessions_to_download(self):
        return self.accessions_to_download

    def set_accessions_to_download(self):

        # need this for local, multiple species actually
        species_so_far = {}
        for taxa in self.taxa_list:

            temp_summary_path = os.path.join(self.save_path, "temp_summary.txt")
            os.system("datasets summary genome taxon \"" + taxa + "\" > '" + temp_summary_path + "'")
            genome_summary_dict = json_to_dict(temp_summary_path)
            os.remove(temp_summary_path)

            for genome_record in genome_summary_dict["reports"]:

                if "refseq_category" in genome_record["assembly_info"] and \
                    genome_record["organism"]["organism_name"] not in species_so_far:

                    blast_organism_dct = {"species": genome_record["organism"]["organism_name"],
                                          "accession": genome_record["accession"]}

                    species_so_far[blast_organism_dct["species"]] = True

                    if blast_organism_dct["species"] not in self.existing_accessions:
                        self.accessions_to_download.append(blast_organism_dct)

    def write_to_species_data_file(self):

        if not self.accessions_to_download:
            return None

        species_data_path = os.path.join(self.genome_storage_path, "species_data.csv")
        f = open(species_data_path, "a")

        species_string = ""
        first = True
        for species in self.accessions_to_download:
            if first:
                species_string = species["species"]
                first = False
            else:
                species_string += "\n" + species["species"]

        name_to_lineage = get_taxonomy_lineage(species_string)

        for species in self.accessions_to_download:
            lineage = list_to_string(name_to_lineage[species["species"]], " ")
            f.write(species["species"] + "," + species["accession"] + "," +
                    lineage + "\n")

        f.close()

    def download_genomes(self):

        for org in self.accessions_to_download:

            # download the genome zip file
            os.system(r"datasets download genome accession " +
                      org["accession"])

            # unzip it into a generic, temporary directory called ncbi_dataset
            os.system("unzip ncbi_dataset.zip")

            # remove zip file
            os.system("rm ncbi_dataset.zip")
            os.system("rm README.md")

            #working_path = subprocess.check_output(["pwd"], shell=True). \
            #    decode("utf-8").strip()
            #os.system("datasets rehydrate --directory " + working_path)

            # get genome fasta file path nested within the temp directory
            genome_file_directory = os.path. \
                join("ncbi_dataset", "data", org["accession"])
            genome_filename = subprocess. \
                check_output(["ls", genome_file_directory]). \
                decode("utf-8").strip()
            genome_filepath = os.path. \
                join(genome_file_directory, genome_filename)

            # create local blast database using the genome file
            blast_db_path = os.path.join(self.genome_storage_path, "blast_db")
            if not os.path.isdir(blast_db_path):
                os.system("mkdir " + blast_db_path)

            species_db = ""
            first = True
            for word in org["species"].split():
                if first:
                    species_db = word
                    first = False
                else:
                    species_db += "_" + word

            new_blast_db_name = os.path.join(blast_db_path, species_db)

            os.system("makeblastdb -dbtype nucl -in " + genome_filepath +
                      " -out " + new_blast_db_name)

            # remove original genome file
            os.system("rm -r ncbi_dataset")


if __name__ == "__main__":


    save_path = "/crun/tony.xie/Downloads"
    taxa_list = ["Tursiops aduncus"]
    genome_storage = "/crun/tony.xie/GenomeStorage"

    downloader = ServerGenomeDownloader(save_path, taxa_list, genome_storage)
    downloader.set_accessions_to_download()
    downloader.download_genomes()

