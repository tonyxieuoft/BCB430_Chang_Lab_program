import os

from Basic_Tools.basic_dictionaries import json_to_dict

class ServerGenomeDownloader:

    def __init__(self, save_path, taxa_list):

        self.save_path = save_path
        self.taxa_list = taxa_list

    def get_genome_accessions(self):

        species_so_far = {}

        for taxa in self.taxa_list:

            temp_summary_path = os.path.join(self.save_path, "temp_summary.txt")
            os.system("datasets summary genome taxon " + taxa + " > '" + temp_summary_path + "'")
            genome_summary_dict = json_to_dict(temp_summary_path)
            os.remove(temp_summary_path)

            blast_organisms_list = []

            for genome_record in genome_summary_dict["reports"]:

                if "refseq_category" in genome_record["assembly_info"] and \
                    genome_record["organism"]["organism_name"] not in species_so_far:

                    blast_organism_dct = {}
                    print(genome_record["organism"]["organism_name"])
                    blast_organism_dct["species"] = \
                        genome_record["organism"]["organism_name"]

                    print(genome_record["accession"])
                    blast_organism_dct["accession"] = \
                        genome_record["accession"]

                    species_so_far[blast_organism_dct["species"]] = True
                    blast_organisms_list.append(blast_organism_dct)

            print(len(blast_organisms_list))


if __name__ == "__main__":

    save_path = "/crun/tony.xie/Downloads"
    taxa_list = ["elasmobranchii"]

    downloader = ServerGenomeDownloader(save_path, taxa_list)
    downloader.get_genome_accessions()

