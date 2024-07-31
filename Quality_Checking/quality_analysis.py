import os

from Basic_Tools.lists_and_files import file_to_list, list_to_string


def get_name(line: str):

    split_line = line.split()
    i = 1
    end = False
    while not end:
        if (len(split_line[i]) >= 4 and split_line[i][:4] == "mRNA") or \
                (len(split_line[i]) >= 14 and split_line[i][:14] == "reference_mrna"):
            end = True
        else:
            i += 1

    return list_to_string(split_line[1:i], " ")


class QualityAnalyser:

    """
    Analyzes alignment quality.
    """
    def __init__(self, filepath):

        self.sequences = {}
        self.species_no = 0
        self.read_fasta_alignment(filepath)

        self.STOP_CODONS = {"TAA": "stop", "TGA": "stop", "TAG": "stop"}

    def read_fasta_alignment(self, filepath):

        f = open(filepath, "r")
        line = f.readline()
        while line != "":
            if line.strip() and line[0] == ">":
                self.species_no += 1
                species_name = get_name(line)
                sequence = f.readline().strip()
                self.sequences[species_name] = sequence
            line = f.readline()

    def get_gaps_mean(self):

        averages_summed = 0
        for species in self.sequences:
            count = 0
            for ch in self.sequences[species]:
                if ch == "-":
                    count += 1

            averages_summed += (count / len(self.sequences[species]))

        return averages_summed / self.species_no

    def get_gaps_median(self):

        gaps_list = []
        for species in self.sequences:
            species_gaps = 0
            for ch in self.sequences[species]:
                if ch == "-":
                    species_gaps += 1
            gaps_list.append(species_gaps / len(self.sequences[species]))

        gaps_list.sort()

        #gaps_list[0], \
        #gaps_list[len(gaps_list) // 2], \
        #gaps_list[len(gaps_list) - 1]

        return gaps_list

    def detect_gaps(self):

        defects = set()
        for species in self.sequences:
            for ch in self.sequences[species]:
                if ch == "-":
                    defects.add(species)
                break

        return defects

    def detect_non_meth_starts(self):

        defects = set()
        for species in self.sequences:
            if len(self.sequences[species]) > 3 and \
                    "-" not in self.sequences[species][:3] and \
                    self.sequences[species][:3] != "ATG":
                defects.add(species)

        return defects

    def detect_premature_stops(self):

        defects = set()
        for species in self.sequences:
            i = 0
            while i < (len(self.sequences[species]) // 3 - 1)*3:
                codon = self.sequences[species][i:i+3]
                if "-" not in codon and codon in self.STOP_CODONS:
                    defects.add(species)
                    break
                i += 3

        return defects

    def detect_non_stop_ends(self):

        defects = set()
        for species in self.sequences:
            if len(self.sequences[species]) % 3 == 0 and \
                    self.sequences[species][-3:] not in self.STOP_CODONS:
                defects.add(species)

        return defects

    def detect_non_three_multiples(self):

        defects = set()
        for species in self.sequences:
            if len(self.sequences[species]) % 3 != 0:
                defects.add(species)

        return defects


if __name__ == "__main__":

    csv_filepath = r"C:\Users\tonyx\Downloads\gaps_all.csv"
    f = open(csv_filepath, "w")

    csv_dct = {}

    alignment_dir1 = r"C:\Users\tonyx\Downloads\main_test4\alignments1"
    alignment_dir2 = r"C:\Users\tonyx\Downloads\main_test4\alignments2"

    for gene in os.listdir(alignment_dir1):
        gene_name = os.path.splitext(gene)[0].upper()

        gene_filepath = os.path.join(alignment_dir1, gene)
        q = QualityAnalyser(gene_filepath)

        csv_dct[gene_name] = q.get_gaps_median()

    for gene_name in csv_dct:
        for value in csv_dct[gene_name]:
            f.write(gene_name + "," + str(value) + "," + "Exon BLAST" + "\n")

    csv_dct = {}

    for gene in os.listdir(alignment_dir2):
        gene_name = os.path.splitext(gene)[0].upper()

        gene_filepath = os.path.join(alignment_dir2, gene)
        q = QualityAnalyser(gene_filepath)

        csv_dct[gene_name] = q.get_gaps_median()

    for gene_name in csv_dct:
        for value in csv_dct[gene_name]:
            f.write(gene_name + "," + str(value) + "," + "Full CDS BLAST" + "\n")

    #file = r"C:\Users\tonyx\Downloads\main_test3\alignments_1e-5_all_hits\GNAT1.fas"
    #q = QualityAnalyser(file)

    #print("mean: " + str(q.get_gaps_mean()))
    #print("quartiles: " + str(q.get_gaps_median()))
    #print("non meth starts: " + str(q.detect_non_meth_starts()))
    #print("premature stops: " + str(q.detect_premature_stops()))
    #print("non stop ends: " + str(q.detect_non_stop_ends()))
    #print("non three multiples: " + str(q.detect_non_three_multiples()))

    #print(q.detect_non_meth_starts() | q.detect_gaps() |
    #      q.detect_non_stop_ends() | q.detect_non_three_multiples() |
    #      q.detect_premature_stops())
