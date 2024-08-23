import os

from Basic_Tools.basic_dictionaries import json_to_dict

if __name__ == "__main__":

    save_path = "/crun/tony.xie/Downloads"

    species_so_far = {}

    for taxa in ["elasmobranchii"]:

        temp_summary_path = os.path.join(save_path, "temp_summary.txt")
        os.system("datasets summary genome taxon " + taxa + " > '" + temp_summary_path + "'")
        genome_summary_dict = json_to_dict(temp_summary_path)
        os.remove(temp_summary_path)

        blast_organisms_list = []

        for genome_record in genome_summary_dict["reports"]:

            if "refseq_category" in genome_record["assembly_info"]:# and \
                    #genome_record["organism"]["organism_name"] not in species_so_far:

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
