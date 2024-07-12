import os

from Basic_Tools.basic_dictionaries import json_to_dict

if __name__ == "__main__":

    os.system("datasets summary genome taxon 7797 > 'temp.txt'")
    genome_summary_dict = json_to_dict("temp.txt")
    os.remove("temp.txt")

    accession = ""
    for genome_record in genome_summary_dict["reports"]:

        print(genome_record["organism"]["organism_name"])
        print(genome_record["accession"])
