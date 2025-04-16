import os

from Bio import Entrez

from Basic_Tools.xml_extraction import file_xml_to_dictionary, get_xml_list
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from NCBI_Genome_Blaster.assemble_blast_result_sequences import \
    ExonBlastXMLParser, SEQUENCE_INDICES_FROM_MRNA_TAG

# maximum intron length for gene models --- from Hara et al. 2018, it appears to be around 1 million for elasmobranchs
MAX_INTRON_LENGTH = 1000000

def extract_query_title(title_str):

    """
    Extracts the gene name, accession and query range for an exon reference query based on an input query heading.
    @param title_str An input query header for an exon reference
    @return Returns the name of the gene corresponding to the exon reference query, the NCBI accession of the
            subject genome that was blasted, and the exon's range in relation to the full coding sequence (e.g., 35-56)
    """

    # separate the query title based on spaces
    query_title = title_str.split(" ")

    # the gene name is the first word in the title
    gene_name = query_title[0]

    # a tracker to determine which part of the title contains the transcript accession for the full coding sequence
    mrna_section_no = 1
    # iteratively check through the words in the title
    while (len(query_title[mrna_section_no]) < len("mRNA") or
           query_title[mrna_section_no][:len("mRNA")] != "mRNA") and \
            (len(query_title[mrna_section_no]) < len("reference") or
             query_title[mrna_section_no][:len("reference")] != "reference"):
        # increment until the mRNA section is found
        mrna_section_no += 1

    # the genomic accession for the subject genome that was blasted is immediately after the "mRNA section"
    acc = query_title[mrna_section_no].split(":")[1]
    # the query range is also a constant distance away from the "mRNA section"
    q_range = query_title[mrna_section_no + SEQUENCE_INDICES_FROM_MRNA_TAG]

    return gene_name, acc, q_range

class ImprovedExonParser(ExonBlastXMLParser):

    """
    A class that extends the ExonBLASTXMLParser general class. See the "assemble_blast_results_sequences.py"
    file for more details. The general improvements made here over the previous approach include the follow:

    1) Consideration of maximum intron length
    2) Length-forcing
    3) Enforcement of HSP co-linearity: if HSP 1 is before HSP 2 in the reference sequence, it must also be before
    HSP 2 in the corresponding genome
    """

    def parse_blast_xml(self):
        """
        Parse an input XML file containing all BLAST results for a given species and all of its reference sequences
        against a single subject genome.
        """

        # track the last gene name encountered (so we know at which iteration we begin at a new gene)
        past_gene_name = ""
        dp_gene_tracker = None

        # convert BLAST XML file into dictionary (the format is exactly the same as the BLAST output file downloaded
        # from NCBI)
        results_dict = file_xml_to_dictionary(self.xml_filepath)

        # Each iteration corresponds to a single exon reference query
        # Convert the iterations into a list (note that "BlastOutput" and "BlastOutput_iterations" are the same names
        # used in the BLAST output file
        exon_iterations = get_xml_list(results_dict['BlastOutput']['BlastOutput_iterations'])

        for exon_iteration in exon_iterations:

            # from the query title, extract gene name, transcript version, reference sequence range
            gene_name, ref_transcript_var, ref_sequence_range = \
                extract_query_title(exon_iteration['Iteration_query-def'])

            # if we've moved onto a new gene
            if gene_name != past_gene_name:

                # if the current gene we just finished iterating through has hsps to parse
                if dp_gene_tracker is not None:
                    self._identify_best_hsps(dp_gene_tracker)

                # reset to the new gene
                dp_gene_tracker = []
                past_gene_name = gene_name

            # go through the hits for that single exon query
            dp_hits_tracker = []
            hits_list = get_xml_list(exon_iteration['Iteration_hits'])

            for hit in hits_list:

                # get the highest scoring HSP in the hit (it is assumed that there aren't multiple relevant hits)
                top_hsp = get_xml_list(hit["Hit_hsps"])[0]
                hsp_attributes = {}

                if self.on_server: # if blast is server mode (downloaded NCBI)
                    hsp_attributes["contig_acc"] = hit['Hit_def'].split(" ")[0]
                else:
                    hsp_attributes["contig_acc"] = hit['Hit_accession']

                # Get attributes for the HSP. Here, "hsp_attributes" represents a dictionary that will contain
                # attributes for a SINGLE HSP

                # the query and subject sequence portions of the HSP
                hsp_attributes["qseq"] = top_hsp["Hsp_qseq"]
                hsp_attributes["hseq"] = top_hsp["Hsp_hseq"]

                # the start and ending indices of the query that the HSP corresponds to
                hsp_attributes["q_start"] = int(top_hsp["Hsp_query-from"])
                hsp_attributes["q_end"] = int(top_hsp["Hsp_query-to"])

                # the sequence range (a single value) of the query
                hsp_attributes["seq_range"] = abs(hsp_attributes["q_start"] - hsp_attributes["q_end"]) + 1

                # IMPORTANT --- the query's raw score, which will be used for our dynamic programming calculations
                hsp_attributes["raw_score"] = int(top_hsp["Hsp_score"])

                # the start and ending indices of the subject that the HSP corresponds to
                hsp_attributes["h_start"] = int(top_hsp["Hsp_hit-from"])
                hsp_attributes["h_end"] = int(top_hsp["Hsp_hit-to"])

                # the length of the subject
                hsp_attributes["contig_len"] = int(hit['Hit_len'])

                # the length of the query
                hsp_attributes["query_len"] = int(exon_iteration['Iteration_query-len'])

                # the accession for the reference transcript (corresponding to the query)
                hsp_attributes["ref_acc"] = ref_transcript_var

                # the range of query in reference to the reference transcript
                hsp_attributes["ref_range"] = ref_sequence_range

                # the name of the gene that the HSP corresponds to
                hsp_attributes["gene_name"] = gene_name

                # if the starting index of the HSP (subject side) is less than the ending index of the HSP, it is on
                # the plus strand
                if hsp_attributes["h_start"] < hsp_attributes["h_end"]:
                    hsp_attributes["strand"] = "Plus"
                else:
                    # otherwise, it's on the minus strand
                    hsp_attributes["strand"] = "Minus"

                # store the HSP dictionary that we have just created in an expanding list (eventually, it will contain
                # all HSPs from a single query
                dp_hits_tracker.append(hsp_attributes)

            # the 'dp_gene_tracker" list stores "hits trackers", each containing all of the hits based on a single
            # exon query
            dp_gene_tracker.append(dp_hits_tracker)

        # once we have collected all hits from all exon references for a single gene, we move onto the stage of
        # identifying the best HSPs for the gene model!
        self._identify_best_hsps(dp_gene_tracker)


    def _identify_best_hsps(self, exons):

        """
        Identifies the best HSPs to select for the gene model. As input, this function takes in a list of lists as
        follows:
            Each index of the outer list corresponds to a single exon query. This index stores an inner list. The
            indices of the inner list each contain a dictionary with information about a single HSP for the particular
            exon query.

        """

        num_exons = len(exons)
        dp_table = []
        for k in range(num_exons):

            hits = exons[k]
            dp_table.append([(0, 0, 0)]*len(hits))

            # instantiate base cases
            if k == 0:
                for i in range(len(hits)):
                    dp_table[k][i] = (hits[i]["raw_score"], -1, -1)

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

                    dp_table[k][i] = (max_coverage + hits[i]["raw_score"],
                                      prev_exon, prev_hit)

        # print(dp_table)
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

        # print(picked)
        #f = open(r"/crun/tony.xie/Downloads/intron_lengths2.txt", "a")
        #f = open(r"C:\Users\tonyx\Downloads\intron_lengths2.txt", "a")
        #f.write(self.species_name + "\n")

        for i in range(1, len(picked)):
            before_hsp = exons[picked[i-1][0]][picked[i-1][1]]
            latter_hsp = exons[picked[i][0]][picked[i][1]]

            before_exon = picked[i-1][0]
            after_exon = picked[i][0]
            exon_diff = after_exon - before_exon

            if (latter_hsp["ref_acc"] == before_hsp["ref_acc"]):
                #f.write(before_hsp["gene_name"] + "_" + str(before_exon) + "_" + str(after_exon) + "\n")
                #f.write(str((abs(latter_hsp["h_start"]-before_hsp["h_end"])-1) // exon_diff) + "\n")
                pass

        #f.close()

        splice_site_dict = {}
        for coord in picked:
            self._length_force_and_print(exons[coord[0]][coord[1]], splice_site_dict)

        # TODO disable SPLICE

        """
        for boundary_key in splice_site_dict:
            if len(splice_site_dict[boundary_key]) == 2:
                print("=====")
                for splice_site in splice_site_dict[boundary_key]:
                    print(splice_site)

                    splice_site_file.write(splice_site["gene"] + "," + splice_site["species"] + "," +
                            str(splice_site["query_junction"]) + "," + splice_site["splice_seq"] + "," +
                            splice_site["left_or_right"] + "," + str(splice_site["fill"]) + splice_site["forced_seq"] + "\n")

        """


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

            return False

            """
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
                
            """

        return True

    def _get_splice_site(self, hsp, subject_boundary, strand, side, fill):

        splice_site = {"gene": "",
                       "species": "",
                       "query_junction": 0,
                       "splice_seq": "",
                       "left_or_right": side,
                       "fill": 0}

        # side refers to what side the splice site is on in comparison to the intron
        # the boundary is the end of the conding region
        if strand == "1" and side == "right":
            arr = ncbi_get_gene_sequence(hsp["contig_acc"], subject_boundary - 2,
                                        subject_boundary - 1, strand)
        elif strand == "1" and side == "left":
            arr = ncbi_get_gene_sequence(hsp["contig_acc"], subject_boundary + 1,
                                        subject_boundary + 2, strand)
        elif strand == "2" and side == "right":
            arr = ncbi_get_gene_sequence(hsp["contig_acc"], subject_boundary + 2,
                                         subject_boundary + 1, strand)
        else:
            arr = ncbi_get_gene_sequence(hsp["contig_acc"], subject_boundary - 1,
                                         subject_boundary - 2, strand)

        for ch in arr:
            splice_site["splice_seq"] += ch

        if side == "left":
            splice_site["query_junction"] = int(hsp["ref_range"].split("-")[1])
        else:
            splice_site["query_junction"] = int(hsp["ref_range"].split("-")[0])-1

        splice_site["gene"] = hsp["gene_name"]
        splice_site["species"] = self.species_name
        splice_site["fill"] = fill

        return splice_site


    def _length_force_and_print(self, hsp, splice_site_dict):

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

        left_allow_splice = False
        left_subject_boundary = 0
        to_salvage = 0
        forced_seq_left = ""
        if missing_left > 0:

            if strand == "1":
                lower_bound = max(hit_bound1 - missing_left, 1)
                upper_bound = hit_bound1 - 1
                to_salvage = upper_bound - lower_bound + 1

                if lower_bound - 2 >= 1:
                    left_allow_splice = True
                    left_subject_boundary = lower_bound

            else:
                lower_bound = min(hit_bound1 + missing_left, hsp["contig_len"])
                upper_bound = hit_bound1 + 1
                to_salvage = lower_bound - upper_bound + 1

                if lower_bound + 2 <= hsp["contig_len"]:
                    left_allow_splice = True
                    left_subject_boundary = lower_bound

            string = ""
            if 0 < to_salvage: # TODO change it back to 0-16 later
                print("getting gene sequence")
                arr = ncbi_get_gene_sequence(hsp["contig_acc"], lower_bound,
                                             upper_bound, strand)
                for ch in arr:
                    string += ch

                forced_seq_left = string

                # TODO this is for inserting gaps
                # string = "-"*(missing_left - to_salvage) + string

            else:
                pass
                #TODO for inserting gaps
                #string = "-"*missing_left + string

            # print("missing left recovered: " + string + " iteration: " + exon_iteration['Iteration_iter-num'])

            result_sequence = string + result_sequence

        else:
            left_subject_boundary = hit_bound1
            left_allow_splice = True


        # TODO disabled SPLICE
        if False and left_allow_splice:
            splice_site = self._get_splice_site(hsp, left_subject_boundary, strand, "right", to_salvage)

            splice_site["forced_seq"] = forced_seq_left + result_sequence[:2]

            query_junction_id = splice_site["query_junction"]
            if query_junction_id in splice_site_dict:
                splice_site_dict[query_junction_id].append(splice_site)
            else:
                splice_site_dict[query_junction_id] = [splice_site]


        right_allow_splice = False
        right_subject_boundary = 0
        to_salvage = 0
        forced_seq_right = ""
        if missing_right > 0:

            if strand == "1":
                lower_bound = hit_bound2 + 1
                upper_bound = min(hit_bound2 + missing_right, hsp["contig_len"])
                to_salvage = upper_bound - lower_bound + 1

                if upper_bound + 2 <= hsp["contig_len"]:
                    right_allow_splice = True
                    right_subject_boundary = upper_bound

            else:
                lower_bound = hit_bound2 - 1
                upper_bound = max(hit_bound2 - missing_right, 1)
                to_salvage = lower_bound - upper_bound + 1

                if upper_bound - 2 >= 1:
                    right_allow_splice = True
                    right_subject_boundary = upper_bound

            string = ""
            # 10 >
            if to_salvage > 0: # change it back later # TODO to between 0 and 16, no?
                print("getting gene seqnece")
                arr = ncbi_get_gene_sequence(hsp["contig_acc"], lower_bound, upper_bound, strand)

                for ch in arr:
                    string += ch

                forced_seq_right = string

                # TODO gap insertion
                #string = string + "-"*(missing_right - to_salvage)
            else:
                pass
                # TODO gap insertion
                #string = string + "-"*missing_right

            # print("missing right recovered: " + string + " iteration: " + exon_iteration['Iteration_iter-num'])

            result_sequence = result_sequence + string

        else:
            right_allow_splice = True
            right_subject_boundary = hit_bound2

        # TODO disable SPLICE
        if False and right_allow_splice:
            splice_site = self._get_splice_site(hsp, right_subject_boundary, strand, "left", to_salvage)

            splice_site["forced_seq"] = result_sequence[-2:] + forced_seq_right

            query_junction_id = splice_site["query_junction"]
            if query_junction_id in splice_site_dict:
                splice_site_dict[query_junction_id].append(splice_site)
            else:
                splice_site_dict[query_junction_id] = [splice_site]

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
    path = r"/Users/tonyx/Downloads/WR27YB5H016-Alignment.xml"
    save_dir = r"/Users/tonyx/Downloads"

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    obj = ImprovedExonParser(path, save_dir,
                             {"taxon": "cetacea", "name": "Orcinus orca"},
                             False)
    obj.parse_blast_xml()

