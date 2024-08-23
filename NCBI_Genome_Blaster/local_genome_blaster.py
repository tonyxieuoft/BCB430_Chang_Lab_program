import subprocess

from typing import List
import os
from Bio import Entrez

from Basic_Tools.basic_dictionaries import json_to_dict
from NCBI_Genome_Blaster.assemble_blast_result_sequences import BlastXMLParser


def local_genome_blaster(save_path: str, queries_path: str,
                         taxa_blast_order: List[str],
                         complete_reference_species: List[str],
                         expect_value):

    working_path = subprocess.check_output(["pwd"], shell=True). \
        decode("utf-8").strip()
    prev_blastdb_path = subprocess.check_output(["echo", "$BLASTDB"], shell=True). \
        decode("utf-8").strip()
    os.environ["BLASTDB"] = prev_blastdb_path + working_path + "/local_blast_db"

    species_so_far = {}

    for species in complete_reference_species:
        species_so_far[species] = True

    for taxa in taxa_blast_order:

        temp_summary_path = os.path.join(save_path, "temp_summary.txt")
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

        for org in blast_organisms_list:

            # download the genome zip file
            os.system(r"datasets download genome accession " +
                      org["accession"] + " --dehydrated")

            # unzip it into a generic, temporary directory called ncbi_dataset
            os.system("unzip ncbi_dataset.zip")

            # remove zip file
            os.system("rm ncbi_dataset.zip")
            os.system("rm README.md")

            os.system("datasets rehydrate --directory " + working_path)

            # get genome fasta file path nested within the temp directory
            genome_file_directory = os.path.\
                join("ncbi_dataset", "data", org["accession"])
            genome_filename = subprocess.\
                check_output(["ls", genome_file_directory]).\
                decode("utf-8").strip()
            genome_filepath = os.path.\
                join(genome_file_directory, genome_filename)

            # create local blast database using the genome file
            os.system("mkdir local_blast_db")
            os.system("makeblastdb -dbtype nucl -in " + genome_filepath +
                      " -out local_blast_db/local")

            # remove original genome file
            os.system("rm -r ncbi_dataset")

            # run blastn on the local database
            reference_filepath = os.path.join(queries_path, taxa + ".fas")
            xml_out_path = os.path.join(save_path, "temp.xml")
            os.system("blastn -db local -outfmt 5 -evalue 0.1 -word_size 11 "
                      "-gapopen 5 -gapextend 2 -reward 2 -penalty -3 -query " +
                      reference_filepath + " > " + xml_out_path)

            # delete the local blast database
            os.system("rm -r local_blast_db")

            # parse the blast results
            parser = BlastXMLParser(xml_out_path, save_path, org)
            parser.parse_blast_xml()
            os.remove(xml_out_path)


if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    #blast_results_path = '/mnt/c/Users/tonyx/Downloads/blast_tests_phototransduction'
    #os.mkdir(blast_results_path)
    #queries_path = "/mnt/c/Users/tonyx/Downloads/'query_files (7)'"
    #taxa_blast_order = ['9738', '9753', '9741', '9732', '9740', '40150', '9756', '119500', '90247', '27609', '9750', '9729', '30558', '9726', '9722', '2746895', '9771', '9721', '378069', '30483', '259919', '117887', '1158979', '117861', '117851', '170819', '30503', '119203', '7778']
    #reference_species = ['Balaenoptera acutorostrata', 'Balaenoptera musculus', 'Balaenoptera ricei', 'Delphinapterus leucas', 'Delphinus delphis', 'Eubalaena glacialis', 'Globicephala melas', 'Kogia breviceps', 'Lagenorhynchus albirostris', 'Lagenorhynchus obliquidens', 'Lipotes vexillifer', 'Mesoplodon densirostris', 'Monodon monoceros', 'Neophocaena asiaeorientalis asiaeorientalis', 'Orcinus orca', 'Phocoena sinus', 'Physeter catodon', 'Tursiops truncatus', 'Amblyraja radiata', 'Carcharodon carcharias', 'Chiloscyllium plagiosum', 'Hemiscyllium ocellatum', 'Hypanus sabinus', 'Leucoraja erinacea', 'Mobula hypostoma', 'Pristis pectinata', 'Rhincodon typus', 'Scyliorhinus canicula', 'Stegostoma tigrinum']

    #local_genome_blaster(blast_results_path, queries_path, taxa_blast_order,
    #                     reference_species, 0.1)

    save_path = r"C:\Users\tonyx\Downloads"

    temp_summary_path = os.path.join(save_path, "temp_summary_2.txt")
    #os.system("datasets summary genome taxon elasmobranchii > '" + temp_summary_path + "'")
    genome_summary_dict = json_to_dict(temp_summary_path)
    os.remove(temp_summary_path)

    blast_organisms_list = []

    for genome_record in genome_summary_dict["reports"]:

        if "refseq_category" in genome_record["assembly_info"]:

            blast_organism_dct = {}
            print(genome_record["organism"]["organism_name"])
            blast_organism_dct["species"] = \
                genome_record["organism"]["organism_name"]

            print(genome_record["accession"])
            blast_organism_dct["accession"] = \
                genome_record["accession"]
