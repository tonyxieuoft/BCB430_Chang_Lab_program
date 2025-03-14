import os

from Bio import Entrez

from Basic_Tools.xml_extraction import file_xml_to_dictionary, get_xml_list
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from NCBI_Genome_Blaster.assemble_blast_result_sequences import \
    ExonBlastXMLParser, SEQUENCE_INDICES_FROM_MRNA_TAG

MAX_INTRON_LENGTH = 1000000
#MAX_CONTIG_GAP = 20000
MAX_OVERLAP = 50

MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_SCORE = -4

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


class ImprovedFullParser(ExonBlastXMLParser):

    def parse_blast_xml(self):

        # convert XML file to dictionary
        results_dict = file_xml_to_dictionary(self.xml_filepath)

        # get iterations (each corresponding to a fasta entry in the query file)
        # this means each iteration corresponds to an FULL query
        full_iterations = get_xml_list(results_dict['BlastOutput']['BlastOutput_iterations'])

        for full_iteration in full_iterations:

            # from the query title, extract gene name, transcript version, reference sequence range
            gene_name, ref_transcript_var, ref_sequence_range = \
                extract_query_title(full_iteration['Iteration_query-def'])

            contig_sep_collection = []

            # go through each hit (which will produce the best possible assembly from that hit)
            hits_list = get_xml_list(full_iteration['Iteration_hits'])

            for hit in hits_list:

                hsp_collection = []

                # we consider all of the HSPs
                hsp_list = get_xml_list(hit["Hit_hsps"])
                for hsp in hsp_list:

                    hsp_attributes = {}

                    if self.on_server:  # if blast is server mode (downloaded NCBI)
                        hsp_attributes["contig_acc"] = hit['Hit_def'].split(" ")[0]
                    else:
                        hsp_attributes["contig_acc"] = hit['Hit_accession']

                    # get attributes for the hsp

                    hsp_attributes["qseq"] = hsp["Hsp_qseq"]
                    hsp_attributes["hseq"] = hsp["Hsp_hseq"]

                    hsp_attributes["q_start"] = int(hsp["Hsp_query-from"])
                    hsp_attributes["q_end"] = int(hsp["Hsp_query-to"])

                    hsp_attributes["seq_range"] = abs(hsp_attributes["q_start"] - hsp_attributes["q_end"]) + 1
                    hsp_attributes["raw_score"] = int(hsp["Hsp_score"])

                    hsp_attributes["h_start"] = int(hsp["Hsp_hit-from"])
                    hsp_attributes["h_end"] = int(hsp["Hsp_hit-to"])

                    hsp_attributes["contig_len"] = int(hit['Hit_len'])
                    hsp_attributes["query_len"] = int(full_iteration['Iteration_query-len'])

                    hsp_attributes["ref_acc"] = ref_transcript_var
                    hsp_attributes["ref_range"] = ref_sequence_range
                    hsp_attributes["gene_name"] = gene_name

                    if hsp_attributes["h_start"] < hsp_attributes["h_end"]:
                        hsp_attributes["strand"] = "Plus"
                    else:
                        hsp_attributes["strand"] = "Minus"

                    hsp_collection.append(hsp_attributes)

                contig_sep_collection.append(hsp_collection)

            self._identify_best_version(contig_sep_collection)  # do it for the last one

    def _identify_best_version(self, contig_sep_collection):
        """
        Taking in all HSP results, organized by contig
        """

        dp_gene_level_results = []

        # we observe each contig separately, generating the DP arrays
        for i in range(len(contig_sep_collection)):

            contig_results = contig_sep_collection[i]

            # sort all of the contigs by query start
            contig_results.sort(key=lambda x: x["q_start"])

            dp_pointers = [None] * len(contig_results)
            dp_best_coverage = [None] * len(contig_results)

            for i in range(len(contig_results)):

                #print(contig_results[i])
                single_coverage = contig_results[i]["raw_score"]
                dp_pointers[i] = -1
                dp_best_coverage[i] = single_coverage
                for j in range(i):
                    if self._is_compatible(contig_results[j], contig_results[i]):
                        if contig_results[i]["q_start"] <= contig_results[j]["q_end"]:
                            overlap_adjustment = self._initial_overlap_score_corrector(contig_results[j],
                                                                                       contig_results[i])
                            curr_coverage = single_coverage + dp_best_coverage[j] + overlap_adjustment
                        else:
                            curr_coverage = single_coverage + dp_best_coverage[j]

                        if curr_coverage > dp_best_coverage[i]:
                            dp_best_coverage[i] = curr_coverage
                            dp_pointers[i] = j

            #print(dp_pointers)
            #print(dp_best_coverage)

            dp_gene_level_results.append({"pointers": dp_pointers, "best_coverage": dp_best_coverage})

        # get the contig with the best gene model
        max_global_coverage = 0
        best_contig = -1
        for i in range(len(dp_gene_level_results)):
            coverage_array = dp_gene_level_results[i]["best_coverage"]
            if max(coverage_array) > max_global_coverage:
                max_global_coverage = max(dp_gene_level_results[i]["best_coverage"])
                best_contig = i

        best_coverage_array = dp_gene_level_results[best_contig]["best_coverage"]
        best_pointers = dp_gene_level_results[best_contig]["pointers"]

        # for the contig with the best gene model, get the last HSP to use
        best_ender = -1
        for i in range(len(best_coverage_array)):
            if best_coverage_array[i] == max_global_coverage:
                best_ender = i

        # generate the traceback chain for stitching together the best sequence
        k = best_ender
        traceback_chain = [k]
        while best_pointers[k] != -1:
            traceback_chain.append(best_pointers[k])
            k = best_pointers[k]

        self._stitch_best_sequence(contig_sep_collection[best_contig],
                                                traceback_chain)

    def _initial_overlap_score_corrector(self, hsp1, hsp2):

        """
        Returns the overlap adjustment (generally a negative value)
        """

        overlap = hsp1["q_end"] - hsp2["q_start"] + 1

        q_array, s1_array = self._convert_hsp_to_array(hsp1)
        s2_array = self._convert_hsp_to_array(hsp2)[1]

        s2_overlap = s2_array[:overlap]
        s1_overlap = s1_array[len(s1_array) - overlap:]
        q_overlap = q_array[len(q_array) - overlap:]

        max_score = self._sort_overlap(q_overlap, s1_overlap, s2_overlap)[1]

        overlap_loss = 0
        for i in range(overlap):
            overlap_loss += self._score_corrector_helper(s1_overlap[i], q_overlap[i])
            overlap_loss += self._score_corrector_helper(s2_overlap[i], q_overlap[i])

        return max_score - overlap_loss

    def _score_corrector_helper(self, s_ch, q_ch):
        """
        Takes in a subject character (if insertion, multiple), as well as the query character, and returns the score
        """
        if s_ch == q_ch:
            return MATCH_SCORE
        else:
            if s_ch == "-":
                return GAP_SCORE
            elif len(s_ch) > 1:
                return GAP_SCORE * (len(s_ch) - 1)
            else:
                return MISMATCH_SCORE


    def _stitch_best_sequence(self, hsp_data, traceback_chain):

        last_hsp = hsp_data[traceback_chain[0]]
        building_array = self._convert_hsp_to_array(last_hsp)[1]
        build_start = last_hsp["q_start"]

        for index in range(1, len(traceback_chain)):

            hsp_to_add = hsp_data[traceback_chain[index]]

            query_to_compare, array_to_add = self._convert_hsp_to_array(hsp_to_add)
            to_add_end = hsp_to_add["q_end"]

            #print("===")
            #print(build_start)
            #print(to_add_end)
            if to_add_end >= build_start:

                overlap = to_add_end - build_start + 1

                build_overlap = building_array[:overlap]
                to_add_overlap = array_to_add[len(array_to_add)-overlap:]
                overlap_query_section = query_to_compare[len(query_to_compare)-overlap:]

                #print(build_overlap)
                #print(to_add_overlap)
                #print(overlap_query_section)

                best_overlap = self._sort_overlap(overlap_query_section, to_add_overlap, build_overlap)[0]
                #print(best_overlap)
                building_array = array_to_add[:len(array_to_add)-overlap] + best_overlap + building_array[overlap:]

            else:
                building_array = array_to_add + building_array

            build_start = hsp_to_add["q_start"]

        self._write_to_file(hsp_data[0], building_array)

    def _sort_overlap(self, query_section, to_add_overlap, build_overlap):
        """
        Given overlapping subject sequences and a query, return the sequence of the best compromise,
        the partition index, and compromise score, in that order
        """

        if to_add_overlap == build_overlap:

            # we just need to gauge the number of matches and mismatches to get the actual raw score
            match_no = 0
            mismatch_no = 0
            for i in range(len(to_add_overlap)):
                if query_section[i] == to_add_overlap[i]:
                    match_no += 1
                else:
                    mismatch_no += 1

            return build_overlap, 0, match_no * MATCH_SCORE + mismatch_no * MISMATCH_SCORE

        else:
            # get start alignment position (we are not interested in all the places that are the same)
            i = 0
            while to_add_overlap[i] == build_overlap[i]:
                i += 1

            # get end alignment position
            j = len(to_add_overlap) - 1
            while to_add_overlap[j] == build_overlap[j]:
                j -= 1

            best_partition = i
            max_score = -float("inf")
            for k in range(i, j + 2):
                score = self._perform_global_alignment(query_section[i:j+1], to_add_overlap[i:k] + build_overlap[k:j+1])

                if score > max_score:
                    max_score = score
                    best_partition = k

            return to_add_overlap[:best_partition] + build_overlap[best_partition:], max_score

    def _perform_global_alignment(self, seq1, seq2):

        seq1 = "".join(seq1)

        temp = ""
        for ch in seq2:
            if ch != "-":
                temp += ch
        seq2 = temp

        seq2 = "".join(seq2)

        gap = GAP_SCORE
        match = MATCH_SCORE
        mismatch = MISMATCH_SCORE

        i_bounds = len(seq1) + 1
        j_bounds = len(seq2) + 1

        dp_array = [None]*i_bounds
        for i in range(i_bounds):
            dp_array[i] = [0] * j_bounds

        for i in range(1, i_bounds):
            dp_array[i][0] = i*gap

        for j in range(1, j_bounds):
            dp_array[0][j] = j*gap

        for i in range(1, i_bounds):
            for j in range(1, j_bounds):
                if seq1[i-1] == seq2[j-1]:
                    dp_array[i][j] = max(dp_array[i][j-1] + gap,
                                         dp_array[i-1][j] + gap,
                                         dp_array[i-1][j-1] + match)
                else:
                    dp_array[i][j] = max(dp_array[i][j - 1] + gap,
                                         dp_array[i - 1][j] + gap,
                                         dp_array[i - 1][j - 1] + mismatch)

        return dp_array[i_bounds - 1][j_bounds - 1]


    def _convert_hsp_to_array(self, hsp):
        """
        Conversion: The hsp_s_array has exactly the same number of indices as the range of the query
        sequence. It collapses insertions in the subject into a neighboring index
        """
        hsp_s_array = []
        hsp_q_array = []
        cache = ""
        for i in range(len(hsp["qseq"])):
            if hsp["qseq"][i] == "-":
                cache += hsp["hseq"][i]
            else:
                hsp_s_array.append(cache + hsp["hseq"][i])
                cache = ""

                hsp_q_array.append(hsp["qseq"][i])

        return hsp_q_array, hsp_s_array


    def _is_compatible(self, x, y):
        """

        :param x: hsp with earlier exon prio
        :param y: hsp with later exon prio
        :return:
        """

        if x["contig_acc"] == y["contig_acc"]:

            if x["strand"] != y["strand"]:
                return False

            # subject sequence considerations
            if x["h_start"] < x["h_end"] and y["h_start"] <= x["h_end"]:
                return False

            if x["h_start"] > x["h_end"] and y["h_start"] >= x["h_end"]:
                 return False

            # intron length
            if abs(y["h_start"] - x["h_end"]) - 1 > MAX_INTRON_LENGTH:
                return False

            # x completely overlaps y in parts of the query
            if x["q_end"] >= y["q_end"]:
                print("e")
                return False

            if x["q_end"] - y["q_start"] + 1 > MAX_OVERLAP:
                print("f")
                return False

        return True

    def _write_to_file(self, hsp, seq_array):

        # get rid of ALL gaps
        result_sequence = ""
        for ch in seq_array:
            if ch != "-":
                result_sequence += ch

        transcript_file = self.create_transcript_file(hsp["ref_acc"], hsp["gene_name"])

        fasta_heading = (">" + hsp["gene_name"] + " " + self.species_name + \
                        " reference_mrna:" + hsp["ref_acc"] + \
                        " genome:" + hsp["contig_acc"] + " ")
                        #+ \
                        #hsp["ref_range"])

        print(fasta_heading + "\n")
        print(result_sequence + "\n")
        #transcript_file.write(fasta_heading + "\n")
        #transcript_file.write(result_sequence + "\n")
        transcript_file.close()


if __name__ == "__main__":

    path = r"/Users/tonyx/Downloads/WDAPHX45013-Alignment.xml"
    save_dir = r"/Users/tonyx/Downloads"

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    obj = ImprovedFullParser(path, save_dir,
                             {"taxon": "selachii", "name": "carcharodon"},
                             False)
    obj.parse_blast_xml()