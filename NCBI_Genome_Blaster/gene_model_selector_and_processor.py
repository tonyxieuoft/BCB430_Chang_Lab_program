"""
This file contains the HSP selection and processing algorithms for the Chang Lab exon-based program. The functions
here are not directly called by the user. Instead, the main program calls on it through an interface. The flow is
outlined in the Methods section of the final report.
"""

import os

from Bio import Entrez

from Basic_Tools.xml_extraction import file_xml_to_dictionary, get_xml_list
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from NCBI_Genome_Blaster.create_basic_XML_processor import \
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
    A class that extends the ExonBLASTXMLParser general class. See the "create_basic_XML_processor.py"
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
            Each index of the outer list corresponds to a single exon query. In particular, these outer list indices
            store an inner list. The indices of the inner list each contain a dictionary with information about a
            single HSP for the particular exon query.

        The dynamic programming algorithm is contained here. We formally define each subproblem D[i][k] as the maximum
        raw score of a compatible gene model containing HSP k at exon i, as well as at most 1 HSP for each of
        exons 1,..,i-1. Notably, this subproblem depends on previous subproblems of smaller sizes, making the process
        highly efficient if we start at zero and move our way up
        """
        # the number of exons in the full coding sequence
        num_exons = len(exons)

        # create a table necessary for dynamic programming
        dp_table = []

        # from one to the number of exons...
        for k in range(num_exons):

            # here, k refers to the k'th exon

            # "hits" represents the inner list described earlier. This variable is a list containing dictionaries
            # regarding all hits for the k'th exon
            hits = exons[k]

            # note that there are three values for each index (the first corresponds to the current maximum score)
            dp_table.append([(0, 0, 0)]*len(hits))

            # Here, it is important to define the dynamic programming model as follows:
            # 1) To form the gene model, we want to select non-clashing HSPs maximizing the raw score
            # 2) Each dynamic programming sub-problem is as follows: we select non-clashing HSPs from the first to the
            #    i'th exon query, where i iterates from 1 to 2 .... to n.

            # instantiate base cases
            if k == 0:

                # for HSPs corresponding to the first exon, there is nothing to "build" off of at the start

                for i in range(len(hits)):

                    # The second and third indices correspond to pointers (the first, a pointer to the exon,
                    # the second, a pointer to the particular HSP. We point to nothing, since we did not build from
                    # anything
                    # The maximum raw score we can obtain from selecting a single HSP from the starting exon query
                    # is just the raw score of the HSP itself

                    dp_table[k][i] = (hits[i]["raw_score"], -1, -1)

            # inductive step: now, we're looking at exon queries beyond the first, and there is something to build off
            else:
                # now, we iterate through the hits for a single exon query
                for i in range(len(hits)):

                    # set initial values (we haven't built off off anything yet, so prev. exon and hit are set to -1)
                    max_coverage = 0
                    prev_exon = -1
                    prev_hit = -1

                    for a in range(k-1, -1, -1): # we iteratively examine the k exons before the current exon

                        # this contains all the HSP from a previous exon
                        past_hits = exons[a]

                        for b in range(len(past_hits)):

                            # see if the current HSP and the previous HSP are mutually compatible. if so, link them
                            # together and lengthen the chain (ie. set the previous exon to point to them)
                            if dp_table[a][b][0] > max_coverage and \
                                    self._is_compatible(exons[a][b], exons[k][i], k-a):
                                max_coverage = dp_table[a][b][0]
                                prev_exon = a
                                prev_hit = b

                    # the "maximum raw score" for the current HSP is its own raw score plus the maximum raw score of
                    # all previous, compatible HSPs that it links to
                    dp_table[k][i] = (max_coverage + hits[i]["raw_score"],
                                      prev_exon, prev_hit)

        # print(dp_table)

        # Here, identify the best HSP that will occur last in the gene model
        # Essentially what we're doing here is just identifying the HSP with the "maximum" raw score, and it will
        # link to all HSPs before it
        best_last_hsp = (-1, -1)
        best_last_score = 0
        # go through the DP table
        for i in range(len(dp_table)):
            for j in range(len(dp_table[i])):
                # observe maximum HSP scores, make changes to "best_last_score" if the maximum score of a
                # particular HSP is greater (meaning that it should be at the end of the gene model)
                if dp_table[i][j][0] > best_last_score:
                    best_last_score = dp_table[i][j][0]
                    best_last_hsp = (i,j)

        # here, we walk through the HSP-pointing chain, starting from the HSP at the end of the gene model and moving
        # ourselves up
        exon_no, hsp_no = best_last_hsp
        picked = []
        while exon_no != -1 and hsp_no != -1:

            # add every HSP that we encounter in the chain; this forms the gene model
            to_insert = (exon_no, hsp_no)
            picked.insert(0, to_insert)

            exon_no, hsp_no = dp_table[exon_no][hsp_no][1], dp_table[exon_no][hsp_no][2]

        # finally, given that we have picked all of our HSPs of interest to form the gene model, we iterate through
        # them, length force them, and print to an output file
        for i in range(len(picked)):
            self._length_force_and_print(exons[picked[i][0]][picked[i][1]])


    def _is_compatible(self, x, y, exon_diff):
        """
        Determines if two HSPs are compatible. This forms the crux of the contraints placed on the dyanmic programming
        algorithm. To reiterate, the following criteria must be met for two HSPs to be compatible
        1) They must be on the same strand
        2) They must be on the same contig
        3) They must be less than MAX_INTRON_LENGTH base pairs a part with reference to their position on the contig
        4) They must be co-linear (same order in the transcript as in the genome)

        :param x: An HSP corresponding to an exon query earlier in the coding sequence
        :param y: An HSP corresponding to an exon query later in the coding sequence
        :param exon_diff: The number of exons in between them. This informs the MAX_INTRON_LENGTH adjustment
        :return: True iff x and y are compatible
        """

        # the HSPs must be on the same contig
        if x["contig_acc"] == y["contig_acc"]:

            # if the HSPs are not on the same strand, return False
            if x["strand"] != y["strand"]:
                return False

            # breach of co-linearity (if they're on the plus strand)
            if x["h_start"] < x["h_end"] and y["h_start"] <= x["h_end"]:
                return False

            # breach of co-linearity (if they're on the minus strand)
            if x["h_start"] > x["h_end"] and y["h_start"] >= x["h_end"]:
                return False

            # if the HSPs are farther apart than the maximum intron length, return False
            if abs(y["h_start"] - x["h_end"])-1 > MAX_INTRON_LENGTH*exon_diff:
                return False

        # only if all the constraints are met, do we reach here
        return True


    def _length_force_and_print(self, hsp):

        # the current prediction, without any length forcing
        result_sequence = hsp["hseq"]

        # boundaries of the exon reference with reference to the full coding sequence
        query_bound1 = int(hsp["q_start"])
        query_bound2 = int(hsp["q_end"])

        # how much we're missing from the front of the HSP (assuming conserved intron-exon structure
        missing_left = query_bound1 - 1

        # how much we're missing from the end of the HSP (again, assuming conserved intron-exon structure)
        missing_right = int(hsp["query_len"]) - query_bound2
        # given this, can easily fill up with "N"s or "-"s

        # the position of the HSP with respect to the subject genome (we have to length force immediately adjacent to
        # these positions)
        hit_bound1 = int(hsp['h_start'])
        hit_bound2 = int(hsp['h_end'])

        # first, determine again what strand the HSP is on
        if hit_bound1 < hit_bound2:
            strand = "1"
        else:
            strand = "2"

        to_salvage = 0
        # if we're missing part of the front of the HSP with reference to the exon, we need to do length forcing
        if missing_left > 0:

            # if we're on the plus strand
            if strand == "1":
                # lower boundary of the base pair indices to force
                # max function so that we don't try to force past the limits of the contig itself
                lower_bound = max(hit_bound1 - missing_left, 1)
                # upper boundary of the base pair indices to force
                upper_bound = hit_bound1 - 1
                # here equivalent to the amount that we are length forcing by
                to_salvage = upper_bound - lower_bound + 1

            # if we're on the minus strand
            else:
                # do the reverse of the above
                # make sure we aren't forcing past the limits of the contig
                lower_bound = min(hit_bound1 + missing_left, hsp["contig_len"])
                upper_bound = hit_bound1 + 1
                to_salvage = lower_bound - upper_bound + 1

            string = ""
            # if the prospective amount to length force is between 0 and 16
            if 0 < to_salvage < 17:
                print("getting gene sequence")

                # call the NCBI API here to length force based on the specified boundaries
                arr = ncbi_get_gene_sequence(hsp["contig_acc"], lower_bound,
                                             upper_bound, strand)

                # format the length-forced segment into a string
                for ch in arr:
                    string += ch

                # TODO this is for inserting gaps
                # string = "-"*(missing_left - to_salvage) + string

            else:
                pass
                #TODO for inserting gaps
                #string = "-"*missing_left + string

            result_sequence = string + result_sequence

        to_salvage = 0
        # if we're missing the of the HSP with reference to the exon
        if missing_right > 0:

            # if we're on the plus strand
            if strand == "1":
                # again, set the boundaries for length-forcing, being wary of not pulling past the contig
                lower_bound = hit_bound2 + 1
                upper_bound = min(hit_bound2 + missing_right, hsp["contig_len"])
                to_salvage = upper_bound - lower_bound + 1

            else:
                # set the boundaries, given that we're on the negative strand
                lower_bound = hit_bound2 - 1
                upper_bound = max(hit_bound2 - missing_right, 1)
                to_salvage = lower_bound - upper_bound + 1

            string = ""
            # same "17 base pair limit" is enforced here
            if 17 > to_salvage > 0:
                print("getting gene seqnece")

                # use NCBI API to pull out the portion to length force from the relevant contig
                arr = ncbi_get_gene_sequence(hsp["contig_acc"], lower_bound, upper_bound, strand)

                # format into string
                for ch in arr:
                    string += ch

                # TODO gap insertion
                #string = string + "-"*(missing_right - to_salvage)
            else:
                pass
                # TODO gap insertion
                #string = string + "-"*missing_right

            # print("missing right recovered: " + string + " iteration: " + exon_iteration['Iteration_iter-num'])

            result_sequence = result_sequence + string


        # create a new output transcript file to store the predictions
        transcript_file = self.create_transcript_file(hsp["ref_acc"], hsp["gene_name"])

        # creates a fasta heading with a specified format: gene name first, then species name,
        # then the accession of the reference transcript, then accession of the subject genome, then the range
        # of the prediction (with respect to the coding sequence)
        fasta_heading = ">" + hsp["gene_name"] + " " + self.species_name + \
                        " reference_mrna:" + hsp["ref_acc"] + \
                        " genome:" + hsp["contig_acc"] + " " + \
                        hsp["ref_range"]

        # write the fasta heading
        transcript_file.write(fasta_heading + "\n")

        # write out the predicted sequence
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

