import os

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


if __name__ == "__main__":
    exons_path = r"C:\Users\tonyx\Downloads\NCBI_length - Copy"
    optimize_transcripts_by_length(exons_path)








