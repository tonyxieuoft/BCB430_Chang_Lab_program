import os

from Bio import Entrez

from Basic_Tools.xml_extraction import file_xml_to_dictionary, get_xml_list
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from NCBI_Genome_Blaster.assemble_blast_result_sequences import \
    ExonBlastXMLParser, SEQUENCE_INDICES_FROM_MRNA_TAG

MAX_INTRON_LENGTH = 10000000
MAX_CONTIG_GAP = 20000

def extract_query_title(title_str):

    query_title = title_str.split(" ")
    gene_name = query_title[0]

    mrna_section_no = 1
    while (len(query_title[mrna_section_no]) < len("mRNA") or
           query_title[mrna_section_no][:len("mRNA")] != "mRNA") and \
            (len(query_title[mrna_section_no]) < len("reference") or
             query_title[mrna_section_no][:len("reference")] != "reference"):
        mrna_section_no += 1

    acc = query_title[mrna_section_no].split(":")[1]
    q_range = query_title[mrna_section_no + SEQUENCE_INDICES_FROM_MRNA_TAG]

    return gene_name, acc, q_range


class ImprovedExonParser(ExonBlastXMLParser):

    def parse_blast_xml(self):

        past_gene_name = ""
        dp_gene_tracker = None

        results_dict = file_xml_to_dictionary(self.xml_filepath)
        exon_iterations = results_dict['BlastOutput']['BlastOutput_iterations']

        for exon_iteration in exon_iterations:

            gene_name, ref_transcript_var, ref_sequence_range = \
                extract_query_title(exon_iteration['Iteration_query-def'])

            if gene_name != past_gene_name:
                if dp_gene_tracker is not None:
                    self._identify_best_hsps(dp_gene_tracker)
                dp_gene_tracker = []
                past_gene_name = gene_name

            dp_hits_tracker = []
            hits_list = get_xml_list(exon_iteration['Iteration_hits'])

            for hit in hits_list:

                top_hsp = get_xml_list(hit["Hit_hsps"])[0]
                hsp_attributes = {}

                if self.on_server: # if blast is server mode (downloaded NCBI)
                    hsp_attributes["contig_acc"] = hit['Hit_def'].split(" ")[0]
                else:
                    hsp_attributes["contig_acc"] = hit['Hit_accession']

                qseq_gaps = top_hsp["Hsp_qseq"].count("-")
                hseq_gaps = top_hsp["Hsp_hseq"].count("-")

                if qseq_gaps == hseq_gaps and qseq_gaps > 0:
                    hsp_attributes["qseq"] = top_hsp["Hsp_qseq"].replace("-", "")
                    hsp_attributes["hseq"] = top_hsp["Hsp_hseq"].replace("-", "")
                else:
                    hsp_attributes["qseq"] = top_hsp["Hsp_qseq"]
                    hsp_attributes["hseq"] = top_hsp["Hsp_hseq"]

                hsp_attributes["q_start"] = int(top_hsp["Hsp_query-from"])
                hsp_attributes["q_end"] = int(top_hsp["Hsp_query-to"])

                hsp_attributes["seq_range"] = abs(hsp_attributes["q_start"] - hsp_attributes["q_end"]) + 1

                hsp_attributes["h_start"] = int(top_hsp["Hsp_hit-from"])
                hsp_attributes["h_end"] = int(top_hsp["Hsp_hit-to"])

                hsp_attributes["contig_len"] = int(hit['Hit_len'])
                hsp_attributes["query_len"] = int(exon_iteration['Iteration_query-len'])

                hsp_attributes["ref_acc"] = ref_transcript_var
                hsp_attributes["ref_range"] = ref_sequence_range
                hsp_attributes["gene_name"] = gene_name

                if hsp_attributes["h_start"] < hsp_attributes["h_end"]:
                    hsp_attributes["strand"] = "Plus"
                else:
                    hsp_attributes["strand"] = "Minus"

                dp_hits_tracker.append(hsp_attributes)

            dp_gene_tracker.append(dp_hits_tracker)


        self._identify_best_hsps(dp_gene_tracker) # do it for the last one


    def _identify_best_hsps(self, exons):

        num_exons = len(exons)
        dp_table = []
        for k in range(num_exons):

            hits = exons[k]
            dp_table.append([(0, 0, 0)]*len(hits))

            # instantiate base cases
            if k == 0:
                for i in range(len(hits)):
                    dp_table[k][i] = (hits[i]["seq_range"], -1, -1)

            else:
                for i in range(len(hits)):

                    max_coverage = 0
                    prev_exon = -1
                    prev_hit = -1

                    for a in range(k-1, -1, -1): # goes down the exons

                        past_hits = exons[a]

                        for b in range(len(past_hits)):

                            # print(str(k) + " " + str(i) + " " + str(a) + " " + str(b))

                            if dp_table[a][b][0] > max_coverage and \
                                    self._is_compatible(exons[a][b], exons[k][i], k-a):
                                max_coverage = dp_table[a][b][0]
                                prev_exon = a
                                prev_hit = b

                    dp_table[k][i] = (max_coverage + hits[i]["seq_range"],
                                      prev_exon, prev_hit)

        print(dp_table)
        best_last_hsp = (-1, -1)
        best_last_score = 0
        for i in range(len(dp_table)):
            for j in range(len(dp_table[i])):
                if dp_table[i][j][0] > best_last_score:
                    best_last_score = dp_table[i][j][0]
                    best_last_hsp = (i,j)

        exon_no, hsp_no = best_last_hsp
        picked = []
        while exon_no != -1 and hsp_no != -1:

            to_insert = (exon_no, hsp_no)
            picked.insert(0, to_insert)

            exon_no, hsp_no = dp_table[exon_no][hsp_no][1], dp_table[exon_no][hsp_no][2]

        print(picked)
        #f = open(r"/crun/tony.xie/Downloads/intron_lengths2.txt", "a")
        f = open(r"C:\Users\tonyx\Downloads\intron_lengths2.txt", "a")
        f.write(self.species_name + "\n")

        for i in range(1, len(picked)):
            before_hsp = exons[picked[i-1][0]][picked[i-1][1]]
            latter_hsp = exons[picked[i][0]][picked[i][1]]

            before_exon = picked[i-1][0]
            after_exon = picked[i][0]
            exon_diff = after_exon - before_exon

            if (latter_hsp["ref_acc"] == before_hsp["ref_acc"]):
                f.write(before_hsp["gene_name"] + "_" + str(before_exon) + "_" + str(after_exon) + "\n")
                f.write(str((abs(latter_hsp["h_start"]-before_hsp["h_end"])-1) // exon_diff) + "\n")

        f.close()

        for coord in picked:
            self._length_force_and_print(exons[coord[0]][coord[1]])



    def _is_compatible(self, x, y, exon_diff):
        """

        :param x: hsp with earlier exon prio
        :param y: hsp with later exon prio
        :param exon_diff: number of exons in between
        :return:
        """

        if x["contig_acc"] == y["contig_acc"]:

            if x["strand"] != y["strand"]:
                return False

            if x["h_start"] < x["h_end"] and y["h_start"] <= x["h_end"]:
                return False

            if x["h_start"] > x["h_end"] and y["h_start"] >= x["h_end"]:
                return False

            if abs(y["h_start"] - x["h_end"])-1 > MAX_INTRON_LENGTH*exon_diff:
                return False

        else: # x["contig_acc"] != y["contig_acc"]


            if x["strand"] == "Plus":
                x_dist = x["contig_len"] - x["h_end"]
            else:
                x_dist = x["h_end"]

            if y["strand"] == "Plus":
                y_dist = y["h_start"]
            else:
                y_dist = y["contig_len"] - y["h_start"]

            if x_dist + y_dist > MAX_CONTIG_GAP:
                return False

        return True

    def _length_force_and_print(self, hsp):

        result_sequence = hsp["hseq"]

        query_bound1 = int(hsp["q_start"])
        query_bound2 = int(hsp["q_end"])

        missing_left = query_bound1 - 1
        missing_right = int(hsp["query_len"]) - query_bound2
        # given this, can easily fill up with "N"s or "-"s

        hit_bound1 = int(hsp['h_start'])
        hit_bound2 = int(hsp['h_end'])
        if hit_bound1 < hit_bound2:
            strand = "1"
        else:
            strand = "2"

        if missing_left > 0:

            if strand == "1":
                lower_bound = max(hit_bound1 - missing_left, 1)
                upper_bound = hit_bound1 - 1
                to_salvage = upper_bound - lower_bound + 1
            else:
                lower_bound = min(hit_bound1 + missing_left, hsp["contig_len"])
                upper_bound = hit_bound1 + 1
                to_salvage = lower_bound - upper_bound + 1

            string = ""
            if 0 < to_salvage < 10: # TODO change it back to 0-5 later
                print("getting gene sequence")
                arr = ncbi_get_gene_sequence(hsp["contig_acc"], lower_bound,
                                             upper_bound, strand)
                for ch in arr:
                    string += ch

                string = "-"*(missing_left - to_salvage) + string

            else:
                string = "-"*missing_left + string

            #print("missing left recovered: " + string + " iteration: " + exon_iteration['Iteration_iter-num'])

            result_sequence = string + result_sequence

        if missing_right > 0:

            if strand == "1":
                lower_bound = hit_bound2 + 1
                upper_bound = min(hit_bound2 + missing_right, hsp["contig_len"])
                to_salvage = upper_bound - lower_bound + 1
            else:
                lower_bound = hit_bound2 - 1
                upper_bound = max(hit_bound2 - missing_right, 1)
                to_salvage = lower_bound - upper_bound + 1

            string = ""
            if 10 > to_salvage > 0: # change it back later
                print("getting gene seqnece")
                arr = ncbi_get_gene_sequence(hsp["contig_acc"], lower_bound, upper_bound, strand)

                for ch in arr:
                    string += ch

                string = string + "-"*(missing_right - to_salvage)
            else:
                string = string + "-"*missing_right

            # print("missing right recovered: " + string + " iteration: " + exon_iteration['Iteration_iter-num'])

            result_sequence = result_sequence + string

        transcript_file = self.create_transcript_file(hsp["ref_acc"], hsp["gene_name"])

        fasta_heading = ">" + hsp["gene_name"] + " " + self.species_name + \
                        " reference_mrna:" + hsp["ref_acc"] + \
                        " genome:" + hsp["contig_acc"] + " " + \
                        hsp["ref_range"]
        transcript_file.write(fasta_heading + "\n")
        transcript_file.write(result_sequence + "\n")
        transcript_file.close()

if __name__ == "__main__":

    #path = r"C:\Users\tonyx\Downloads\KGMVRFG9013-Alignment.xml"
    #path = r"C:\Users\tonyx\Downloads\M2GBJPF2013-Alignment.xml"
    #path = r"C:\Users\tonyx\Downloads\M2HJANYT016-Alignment.xml"
    path = r"C:\Users\tonyx\Downloads\M3DS2H6D013-Alignment.xml"
    save_dir = r"C:\Users\tonyx\Downloads"

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    obj = ImprovedExonParser(path, save_dir,
                             {"taxon": "cetacea", "name": "Orcinus orca"},
                             False)
    obj.parse_blast_xml()

