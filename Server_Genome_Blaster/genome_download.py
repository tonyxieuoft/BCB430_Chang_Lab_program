import os

from Basic_Tools.basic_dictionaries import json_to_dict

######### species_data file constants #########
from Basic_Tools.lists_and_files import list_to_string
from Basic_Tools.taxonomy_browser import get_taxonomy_lineage

NAME_COL = 0
ACC_COL = 1
LINEAGE_COL = 2


def read_species_data(genome_storage_path):

    species_data_path = os.path.join(genome_storage_path, "species_data.csv")

    if not os.path.isfile(species_data_path):
        return None

    species_data = []
    f = open(species_data_path, "r")

    line = f.readline()

    while line != "":

        split_line = line.strip().split(",")
        genome = {"name": split_line[NAME_COL], "acc": split_line[ACC_COL],
                  "lineage": split_line[LINEAGE_COL]}

        species_data.append(genome)
        line = f.readline()

    return species_data


class ServerGenomeDownloader:

    def __init__(self, save_path, taxa_list, genome_storage_path):

        self.save_path = save_path
        self.taxa_list = taxa_list
        self.genome_storage_path = genome_storage_path

        self.existing_accessions = {}
        self.accessions_to_download = []

    def get_existing_accessions(self):

        species_data = read_species_data(self.genome_storage_path)
        if species_data is None:
            return None

        for species in species_data:
            self.existing_accessions[species["acc"]] = True

    def get_accessions_to_download(self):

        # need this for local, multiple species actually
        species_so_far = {}
        for taxa in self.taxa_list:

            temp_summary_path = os.path.join(self.save_path, "temp_summary.txt")
            os.system("datasets summary genome taxon " + taxa + " > '" + temp_summary_path + "'")
            genome_summary_dict = json_to_dict(temp_summary_path)
            os.remove(temp_summary_path)

            for genome_record in genome_summary_dict["reports"]:

                if "refseq_category" in genome_record["assembly_info"] and \
                    genome_record["organism"]["organism_name"] not in species_so_far:

                    blast_organism_dct = {"species": genome_record["organism"]["organism_name"],
                                          "accession": genome_record["accession"]}

                    species_so_far[blast_organism_dct["species"]] = True

                    if blast_organism_dct["accession"] not in self.existing_accessions:
                        self.accessions_to_download.append(blast_organism_dct)

    def write_to_species_data_file(self):

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


if __name__ == "__main__":

    save_path = "/crun/tony.xie/Downloads"
    taxa_list = ["batoidea", "elasmobranchii"]
    genome_storage = "/crun/tony.xie/GenomeStorage"

    downloader = ServerGenomeDownloader(save_path, taxa_list, genome_storage)
    downloader.get_accessions_to_download()
    downloader.write_to_species_data_file()

    print(downloader.accessions_to_download)
    print(len(downloader.accessions_to_download))

