import json
import os

from Bio import Entrez

from Basic_Tools.lists_and_files import list_to_string


def get_longest_transcript(species_dir: str):
    """

    :param species_dir: "species" level directory within an NCBI exon
    puller-generated general directory
    :return:
    """

    mx = 0
    mx_name = ""
    for t_name in os.listdir(species_dir):
        t_length = t_name.split("_")[0]

        if t_length.isnumeric():
            if int(t_length) > mx:
                mx = int(t_length)
                mx_name = t_name

        else: # it's a blast transcript
            mx_name = t_name
            f = open(os.path.join(species_dir, t_name), "r")
            line = f.readline()
            while line != "":
                if line and line[0] == ">":
                    mx = int(line.split()[-1].split("-")[1])
                line = f.readline()

    return {"length": mx, "file_name": mx_name}


def filter_closest_length_transcript(species_dir: str, target: int):

    t_list = os.listdir(species_dir)

    min_diff = float('inf')*-1
    closest_name = ""

    for t_name in t_list:
        t_length = int(t_name.split("_")[0])

        diff = t_length - target

        if diff >= 0 and (min_diff < 0 or diff < min_diff):
            min_diff = diff
            closest_name = t_name
        elif diff < 0 and min_diff < 0 and abs(diff) < abs(min_diff):
            min_diff = diff
            closest_name = t_name

    for t_name in t_list:
        if t_name != closest_name:
            os.remove(os.path.join(species_dir, t_name))


def remove_non_three_multiples(species_dir):

    transcripts = os.listdir(species_dir)
    for t_name in transcripts:
        t_length = int(t_name.split("_")[0])
        if t_length % 3 != 0:
            os.remove(os.path.join(species_dir, t_name))

# add the ey filter detection here
def optimize_transcripts_by_length(general_dir):
    """

    :param general_dir: an NCBI exon puller-generated directory
    :return:
    """

    for gene in os.listdir(general_dir):
        gene_path = os.path.join(general_dir, gene)

        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            lt_lengths = []

            species_list = os.listdir(taxon_path)
            for species in species_list:
                species_path = os.path.join(taxon_path, species)

                remove_non_three_multiples(species_path)
                if len(os.listdir(species_path)) == 0:
                    os.rmdir(species_path)
                else:
                    lt_lengths.append(get_longest_transcript(species_path)["length"])

            median_length = -1
            lt_lengths.sort()
            if len(lt_lengths) > 0:
                if len(lt_lengths) % 2 == 0:
                    median_length = lt_lengths[len(lt_lengths)//2-1]
                else:
                    median_length = lt_lengths[len(lt_lengths) // 2]

            for species in os.listdir(taxon_path):
                species_path = os.path.join(taxon_path, species)

                filter_closest_length_transcript(species_path, median_length)


def get_transcript_tissue_origin(exon_pull_dir):

    for gene in os.listdir(exon_pull_dir):

        acc_list = []
        gene_path = os.path.join(exon_pull_dir, gene)
        for taxon in os.listdir(gene_path):

            taxon_path = os.path.join(gene_path, taxon)
            for species in os.listdir(taxon_path):

                species_path = os.path.join(taxon_path, species)
                for transcript_name in os.listdir(species_path):
                    transcript_name = os.path.splitext(transcript_name)[0]
                    acc = list_to_string(transcript_name.split("_")[1:], "_")
                    acc_list.append(acc)

                    #print(gene + " " + species)

        acc_string = list_to_string(acc_list, ",")
        print(acc_string)
        print(gene)
        print_tissue_origin(acc_string)


def print_tissue_origin(acc):

    # XM_048540187.2,XM_059967774.1
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    tissues = {}

    tissue_typing_handle = Entrez.esummary(db="nuccore", id=acc, version="2.0")
    records = Entrez.read(tissue_typing_handle)
    for record in records["DocumentSummarySet"]:
        i = 0
        found = False
        acc = ""
        while i < len(record) and not found:
            if i == 0:
                acc = record[i]
            else:
                if isinstance(record[i], str):
                    line = record[i].split("|")
                    for j in range(len(line)):
                        if line[j] == "tissue_type":
                            tissue_type = record[i+1].split("|")[j]
                            #print(tissue_type)
                            if tissue_type not in tissues:
                                tissues[tissue_type] = ""
                            found = True
                            break
            i += 1
        #if not found:
        #    print("not found")
        #print("acc: " + acc)

    print(tissues)


if __name__ == "__main__":
    #exons_path = r"C:\Users\tonyx\Downloads\NCBI_length - Copy"
    #optimize_transcripts_by_length(exons_path)

    exon_pull_dir = r"/crun/tony.xie/Downloads/NCBI_exon_pull_results2"
    get_transcript_tissue_origin(exon_pull_dir)





