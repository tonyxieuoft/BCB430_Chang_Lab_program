import os
from Basic_Tools.basic_dictionaries import json_to_dict
import subprocess


if __name__ == "__main__":
    # tempjson.txt location will instead be the folder the user specifies
    temp_file_name = r'../tests/tempjson.txt'
    os.system("datasets summary genome taxon 'orcinus orca' > " + temp_file_name)
    dct = json_to_dict(temp_file_name)

    accession = ""
    species = ""
    for genome_record in dct["reports"]:
        if "refseq_category" in genome_record["assembly_info"]:
            print(genome_record["organism"]["organism_name"])
            species = genome_record["organism"]["organism_name"]

            print(genome_record["accession"])
            accession = genome_record["accession"]

    os.system(r"datasets download genome accession " + accession)
    os.system("unzip ncbi_dataset.zip")
    os.system("rm ncbi_dataset.zip")
    os.system("rm README.md")
    genome_file_folder = "ncbi_dataset/data/" + accession
    genome_file = subprocess.check_output(["ls", genome_file_folder]).\
        decode("utf-8").strip()
    print("This is the genome file: " + genome_file)
    os.system("mkdir orca")
    os.system("makeblastdb -dbtype nucl -in " + genome_file_folder +
              "/" + genome_file + " -out orca/orca")

    # this is not working for some reason
    curr_path = subprocess.check_output(["pwd"], shell=True). \
        decode("utf-8").strip()
    prev_blastdb_path = subprocess.check_output(["echo", "$BLASTDB"], shell=True).\
        decode("utf-8").strip()
    os.environ["BLASTDB"] = prev_blastdb_path + curr_path + "/orca"
    os.system("rm -r ncbi_dataset")
    query = "/mnt/c/Users/tonyx/Downloads/'query_files (2)'/9721.fas"
    os.system("blastn -db orca -outfmt 5 -query " + query + " > test.xml")


