import os

from After_BLAST.concatenate_exons import concatenate_exons
from Basic_Tools.lists_and_files import file_to_list, make_unique_directory

batoidea = ["Leucoraja erinaceus",
           "Leucoraja ocellata",
           "Amblyraja radiata",
           "Raja brachyura",
           "Malacoraja senta",
           "Okamejei kenojei",
           "Narcine bancroftii",
           "Torpedo suessi",
           "Rhina ancyclostomus",
           "Rhynchobatus australiae",
           "Pristis pectinata",
           "Hypanus sabinus",
           "Hypanus berthalutzae",
           "Mobula hypostoma",
           "Mobula birostris"]

selachii = ["Squalus suckleyi",
            "Squalus acanthias",
            "Squatina squatina",
            "Heptranchias perlo",
            "Heterodontus francisci",
            "Stegostoma tigrinum",
            "Rhincodon typus",
            "Ginglymostoma cirratum",
            "Hemiscyllium ocellatum",
            "Chiloscyllium punctatum",
            "Chiloscyllium plagiosum",
            "Carcharodon carcharias",
            "Cetorhinus maximus",
            "Isurus oxyrinchus",
            "Mustelus asterias",
            "Scyliorhinus canicula",
            "Scyliorhinus torazame",
            "Sphyrnna mokarran",
            "Negaprion brevirostris",
            "Carcharinus longimanus",
            "Prionace glauca",
            "Pristiophorus japonicus"]

def _global_aligner(seq1, seq2):
    go = -8
    ge = -1
    ma = 2
    mi = -2

    i_range = len(seq1) + 1
    j_range = len(seq2) + 1
    dp = [None] * i_range
    for i in range(i_range):
        dp[i] = [None] * j_range
        for j in range(j_range):
            dp[i][j] = {"sub": 0,
                        "gap_i": 0,
                        "gap_ip": None,
                        "gap_j": 0,
                        "gap_jp": None}

    for i in range(1, i_range):
        dp[i][0]["sub"] = -float("inf")
        dp[i][0]["gap_i"] = -float("inf")
        dp[i][0]["gap_j"] = go + ge * (i - 1)
        dp[i][0]["gap_jp"] = "gap_j"

    for j in range(1, j_range):
        dp[0][j]["sub"] = -float("inf")
        dp[0][j]["gap_i"] = go + ge * (j - 1)
        dp[0][j]["gap_ip"] = "gap_i"
        dp[0][j]["gap_j"] = -float("inf")

    for i in range(1, i_range):
        for j in range(1, j_range):

            # sub
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j]["sub"] = max(dp[i - 1][j - 1]["sub"],
                                      dp[i - 1][j - 1]["gap_i"],
                                      dp[i - 1][j - 1]["gap_j"]) + ma
            else:
                dp[i][j]["sub"] = max(dp[i - 1][j - 1]["sub"],
                                      dp[i - 1][j - 1]["gap_i"],
                                      dp[i - 1][j - 1]["gap_j"]) + mi

            # gap_i
            dp[i][j]["gap_i"] = max(dp[i][j - 1]["sub"] + go,
                                    dp[i][j - 1]["gap_i"] + ge,
                                    dp[i][j - 1]["gap_j"] + go)

            if dp[i][j]["gap_i"] == dp[i][j - 1]["sub"] + go:
                dp[i][j]["gap_ip"] = "sub"
            elif dp[i][j]["gap_i"] == dp[i][j - 1]["gap_i"] + ge:
                dp[i][j]["gap_ip"] = "gap_i"
            else:
                dp[i][j]["gap_ip"] = "gap_j"

            # gap_j
            dp[i][j]["gap_j"] = max(dp[i - 1][j]["sub"] + go,
                                    dp[i - 1][j]["gap_i"] + go,
                                    dp[i - 1][j]["gap_j"] + ge)

            if dp[i][j]["gap_j"] == dp[i-1][j]["sub"] + go:
                dp[i][j]["gap_jp"] = "sub"
            elif dp[i][j]["gap_j"] == dp[i-1][j]["gap_i"] + go:
                dp[i][j]["gap_jp"] = "gap_i"
            else:
                dp[i][j]["gap_jp"] = "gap_j"

    seq1_aligned = ""
    seq2_aligned = ""

    i, j = i_range - 1, j_range - 1
    gap_recall = None
    while i != 0 or j != 0:

        if gap_recall is not None:
            if gap_recall == "sub":
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = seq2[j - 1] + seq2_aligned
                i -= 1
                j -= 1

                gap_recall = None

            elif gap_recall == "gap_i":

                seq1_aligned = "-" + seq1_aligned
                seq2_aligned = seq2[j - 1] + seq2_aligned

                gap_recall = dp[i][j]["gap_ip"]
                j -= 1

            else:

                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = "-" + seq2_aligned

                gap_recall = dp[i][j]["gap_jp"]

                i -= 1

        else:

            box_max = max(dp[i][j]["sub"], dp[i][j]["gap_i"], dp[i][j]["gap_j"])
            if box_max == dp[i][j]["sub"]:
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = seq2[j - 1] + seq2_aligned
                i -= 1
                j -= 1
            elif box_max == dp[i][j]["gap_i"]:
                seq1_aligned = "-" + seq1_aligned
                seq2_aligned = seq2[j - 1] + seq2_aligned

                gap_recall = dp[i][j]["gap_ip"]

                j -= 1
            else:
                seq1_aligned = seq1[i - 1] + seq1_aligned
                seq2_aligned = "-" + seq2_aligned

                gap_recall = dp[i][j]["gap_jp"]

                i -= 1

    print(seq1_aligned)
    print(seq2_aligned)
    print(dp[i_range - 1][j_range - 1])

class Analyser:

    def __init__(self, alignment_dir, wd):

        self.alignment_dir = alignment_dir
        self.wd = wd

        self.gene_species_sequence_dict = {}
        for gene_file in os.listdir(self.alignment_dir):

            gene_name = os.path.splitext(gene_file)[0]
            gene_path = os.path.join(self.alignment_dir, gene_file)

            species_sequence_dict = self._read_alignment(gene_path)
            self.gene_species_sequence_dict[gene_name] = species_sequence_dict

    def _read_alignment(self, gene_path):

        species_sequence_dict = {}

        alignment_arr = file_to_list(gene_path)

        curr_species = ""
        for line in alignment_arr:
            if line[0] == ">":
                #curr_species = line[1:]
                curr_species = line.split()[1] + " " + line.split()[2]
            else:
                seq = ""
                for ch in line.strip():
                    if ch != "-":
                        seq += ch
                species_sequence_dict[curr_species] = seq

        return species_sequence_dict

    def _mafft_pairwise_aligner(self,seq1, seq2):

        new_file = os.path.join(self.wd, "tmp.fas")
        mafft_output = os.path.join(self.wd, "tmp_mafft.fas")

        f = open(new_file, "w")
        f.write(">seq1\n")
        f.write(seq1 + "\n")
        f.write(">seq2\n")
        f.write(seq2 + "\n")
        f.close()

        os.system("mafft --auto --quiet --thread -4 " + new_file + " > " + mafft_output)

        f = open(mafft_output, "r")
        line = f.readline().strip()
        seqs = []
        seq = ""
        while line != "":
            if line[0] == ">":
                if seq != "":
                    seqs.append(seq)
                    seq = ""
            else:
                seq += line

            line = f.readline().strip()

        seqs.append(seq)

        os.system("rm " + new_file)
        os.system("rm " + mafft_output)

        # print(file_friendly_species)
        # print("length1: " + str(len(seqs[0])) + " length2: " + str(len(seqs[1])))

        gap_seq1 = 0
        gap_seq2 = 0
        match = 0
        mismatch = 0

        for i in range(len(seqs[0])):
            if seqs[0][i] == seqs[1][i]:
                match += 1
            elif seqs[0][i] == "-":
                gap_seq1 += 1
            elif seqs[1][i] == "-":
                gap_seq2 += 1
            else:
                mismatch += 1
        # print("gap_full: " + str(gap_full) +
        #      "gap_exon: " + str(gap_exon) +
        #      "match: " + str(match) +
        #      "mismatch: " + str(mismatch))

        return seqs[0], seqs[1], gap_seq1, gap_seq2, match, mismatch

    def _align_against(self, other, csv_for_r):

        csv_f = open(csv_for_r, "w")
        csv_f.write("gene,gap_full,gap_exon,match,mismatch\n")

        for gene in self.gene_species_sequence_dict:
            print("=====")
            print(gene)

            if gene in other.gene_species_sequence_dict:
                ss_dict1 = self.gene_species_sequence_dict[gene]
                ss_dict2 = other.gene_species_sequence_dict[gene]

                gene_gap_full = 0
                gene_gap_exon = 0
                gene_match = 0
                gene_mismatch = 0

                species_considered = 0

                for species in ss_dict1:

                    if species in ss_dict2:

                        a_seq1, a_seq2, gap_full, gap_exon, match, mismatch = (
                            self._mafft_pairwise_aligner(ss_dict1[species], ss_dict2[species]))

                        gene_gap_full += gap_full / len(a_seq1)
                        gene_gap_exon += gap_exon / len(a_seq1)
                        gene_match += match / len(a_seq1)
                        gene_mismatch += mismatch / len(a_seq1)

                        species_considered += 1

                formatted_gap_full = "{:.5f}".format(gene_gap_full / species_considered)
                formatted_gap_exon = "{:.5f}".format(gene_gap_exon / species_considered)
                formatted_match = "{:.5f}".format(gene_match / species_considered)
                formatted_mismatch = "{:.5f}".format(gene_mismatch / species_considered)

                print(formatted_gap_full)
                print(formatted_gap_exon)
                print(formatted_match)
                print(formatted_mismatch)

                csv_f.write(gene + "," +
                            formatted_gap_full + "," +
                            formatted_gap_exon + "," +
                            formatted_match + "," +
                            formatted_mismatch + "\n")


                        #print("matching species: " + species)
                        #_global_aligner(ss_dict1[species], ss_dict2[species])

    def _align_against_reference(self, converted_NEPR, csv_for_r):

        csv_f = open(csv_for_r, "w")
        csv_f.write("gene,gap_full,gap_exon,match,mismatch\n")

        for gene in self.gene_species_sequence_dict:

            print("====")
            print(gene)

            shark_ref = None
            batoid_ref = None

            gene_gap_full = 0
            gene_gap_exon = 0
            gene_match = 0
            gene_mismatch = 0

            species_considered = 0

            NEPR_gene_dir = os.path.join(converted_NEPR, gene)
            for file in os.listdir(NEPR_gene_dir):
                file_path = os.path.join(NEPR_gene_dir, file)
                ref_species = os.path.splitext(file)[0]

                if ref_species == "Carcharodon carcharias":
                    shark_ref = concatenate_exons(file_path).split("\n")[1]
                else:
                    batoid_ref = concatenate_exons(file_path).split("\n")[1]

            ss_dict = self.gene_species_sequence_dict[gene]
            for species in ss_dict:
                seq1 = ss_dict[species]
                if species in selachii:
                    if shark_ref is not None:
                        a_seq1, a_seq2, gap_full, gap_exon, match, mismatch = (
                            self._mafft_pairwise_aligner(seq1, shark_ref))
                    else:
                        a_seq1, a_seq2, gap_full, gap_exon, match, mismatch = (
                            self._mafft_pairwise_aligner(seq1, batoid_ref))
                else:
                    if batoid_ref is not None:
                        a_seq1, a_seq2, gap_full, gap_exon, match, mismatch = (
                            self._mafft_pairwise_aligner(seq1, batoid_ref))
                    else:
                        a_seq1, a_seq2, gap_full, gap_exon, match, mismatch = (
                            self._mafft_pairwise_aligner(seq1, shark_ref))

                gene_gap_full += gap_full / len(a_seq1)
                gene_gap_exon += gap_exon / len(a_seq1)
                gene_match += match / len(a_seq1)
                gene_mismatch += mismatch / len(a_seq1)

                species_considered += 1

            formatted_gap_full = "{:.5f}".format(gene_gap_full / species_considered)
            formatted_gap_exon = "{:.5f}".format(gene_gap_exon / species_considered)
            formatted_match = "{:.5f}".format(gene_match / species_considered)
            formatted_mismatch = "{:.5f}".format(gene_mismatch / species_considered)

            print(formatted_gap_full)
            print(formatted_gap_exon)
            print(formatted_match)
            print(formatted_mismatch)

            csv_f.write(gene + "," +
                        formatted_gap_full + "," +
                        formatted_gap_exon + "," +
                        formatted_match + "," +
                        formatted_mismatch + "\n")

    def phylo_tree_generation(self, converted_NEPR):

        save_dir = make_unique_directory(self.wd, "raw_adjusted")

        for gene_file in os.listdir(self.alignment_dir):

            gene_path = os.path.join(self.alignment_dir, gene_file)

            #print("gene path: " + gene_path)
            gene_name = os.path.splitext(gene_file)[0]


            shark_ref = ""
            batoid_ref = ""
            NEPR_gene_dir = os.path.join(converted_NEPR, gene_name)
            for file in os.listdir(NEPR_gene_dir):
                file_path = os.path.join(NEPR_gene_dir, file)
                ref_species = os.path.splitext(file)[0]

                if ref_species == "Carcharodon carcharias":
                    shark_ref = (">Reference " + ref_species + "\n" +
                                 concatenate_exons(file_path).split("\n")[1] + "\n")
                else:
                    batoid_ref = (">Reference " + ref_species + "\n" +
                                 concatenate_exons(file_path).split("\n")[1] + "\n")

            in_f = open(gene_path, "r")
            out_f = open(os.path.join(save_dir, gene_file), "w")

            if shark_ref != "":
                out_f.write(shark_ref)
            if batoid_ref != "":
                out_f.write(batoid_ref)

            line = in_f.readline()
            while line != "":
                if line[0] == ">":
                    species_name = line.split()[1] + " "+ line.split()[2]
                    out_f.write(">" + species_name + "\n")
                else:
                    out_f.write(line)

                line = in_f.readline()

            in_f.close()
            out_f.close()

        mafft_folder = make_unique_directory(self.wd, "mafft")
        for file in os.listdir(save_dir):

            #print("raw_file_name: " + file)

            in_path = os.path.join(save_dir, file)
            out_path = os.path.join(mafft_folder, file)
            os.system("mafft --auto " + in_path + " > " + out_path)


        phylo_folder = make_unique_directory(self.wd, "phylo")
        for file in os.listdir(mafft_folder):

            in_path = os.path.join(mafft_folder, file)

            os.system("nice -2 /usr/local/bin/iqtree2 "
                      "-s " + in_path + " -pre " + phylo_folder + "/" + file + " -T 16")









if __name__ == "__main__":
    #alignment_path1 = "/Users/tonyx/Documents/chang_lab/1-1-4-1-e1-full-raw-score-filtered-fixed"
    #alignment_path2 = "/Users/tonyx/Documents/chang_lab/1-1-4-1-e0.5-exon_ws11_nf"
    #alignment_path3 = "/Users/tonyx/Documents/chang_lab/1-1-4-1-e0.5-exon-raw-score-no_force"

    alignment_path1 = "/Users/tonyx/Documents/chang_lab/BCB430_final_results/GeMoMa_vision_results"
    alignment_path2 = "/Users/tonyx/Documents/chang_lab/BCB430_final_results/exon_ws11_f16"

    # just play around in
    wd = "/Users/tonyx/Documents/chang_lab/mafft_test"
    full = Analyser(alignment_path1, wd)
    exon1 = Analyser(alignment_path2, wd)
    #exon2 = Analyser(alignment_path3, wd)

    ref_dir = "/Users/tonyx/Documents/chang_lab/converted_NEPR"

    csv_for_r = "/Users/tonyx/Documents/chang_lab/GeMoMa_vs_ref.csv"

    #full._align_against(exon1, csv_for_r)
    full._align_against_reference(ref_dir, csv_for_r)
    #_global_aligner("bruhasdfaskjfewj", "asdfasifew")









