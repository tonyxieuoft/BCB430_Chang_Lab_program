import os
from abc import abstractmethod

from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from Bio import Entrez
from Basic_Tools.xml_extraction import file_xml_to_dictionary
from Basic_Tools.xml_extraction import get_xml_list

SEQUENCE_INDICES_FROM_MRNA_TAG = 2


def get_directory(parent_directory: str, directory_name: str) -> str:
    """
    Returns the directory specified by

    :param parent_directory:
    :param directory_name:
    :return:
    """

    directory_path = os.path.join(parent_directory, directory_name)
    if not os.path.isdir(directory_path):
        os.mkdir(directory_path)

    return directory_path


def gap_corrector(q_seq, h_seq):

    q_seq_gaps = 0
    q_no_gaps = ""
    for ch in q_seq:
        if ch == "-":
            q_seq_gaps += 1
        else:
            q_no_gaps += ch

    q_gaps_counter = q_seq_gaps
    h_seq_gaps = 0
    h_corrected_gaps = ""
    for ch in h_seq:
        if ch == "-":
            h_seq_gaps += 1
            q_gaps_counter -= 1
        if ch != "-" or q_gaps_counter < 0:
            h_corrected_gaps += ch

    return h_corrected_gaps


class BlastXMLParser:

    def __init__(self, file, save_dir, curr_species):

        self.xml_filepath = file
        self.save_dir = save_dir

        self.taxon_name = curr_species['taxon']
        self.species_name = curr_species['name']

    def create_transcript_file(self, ref_transcript_var, gene_name):

        gene_folder = get_directory(self.save_dir, gene_name)
        taxa_folder = get_directory(gene_folder, self.taxon_name)
        species_folder = get_directory(taxa_folder, self.species_name)

        transcript_filepath = os.path.join(species_folder,
                                           "_Reference_" +
                                           ref_transcript_var + ".fas")
        return open(transcript_filepath, "a")

    @abstractmethod
    def parse_blast_xml(self):
        pass


class ExonBlastXMLParser(BlastXMLParser):

    def parse_blast_xml(self):

        results_dict = file_xml_to_dictionary(self.xml_filepath)

        exon_iterations = results_dict['BlastOutput']['BlastOutput_iterations']
        for exon_iteration in exon_iterations:

            query_title = exon_iteration['Iteration_query-def'].split(" ")
            query_seq_length = int(exon_iteration['Iteration_query-len'])
            gene_name = query_title[0]

            mrna_section_no = 1
            while len(query_title[mrna_section_no]) < len("mRNA") or \
                    query_title[mrna_section_no][:len("mRNA")] != "mRNA":
                mrna_section_no += 1

            ref_transcript_var = query_title[mrna_section_no].split(":")[1]
            ref_sequence_range = query_title[mrna_section_no +
                                             SEQUENCE_INDICES_FROM_MRNA_TAG]

            hits_list = get_xml_list(exon_iteration['Iteration_hits'])
            if hits_list: # if not empty

                # get the top hit for the exon
                top_hit = hits_list[0]

                accession = top_hit['Hit_accession']
                # accession = top_hit['Hit_def'].split(" ")[0]
                hit_max = int(top_hit['Hit_len'])

                # if hit, there exists at least one hsp
                seq_data = get_xml_list(top_hit["Hit_hsps"])[0]

                result_sequence = seq_data["Hsp_hseq"]

                """
                temp_filter_out_gaps = ""
                for ch in result_sequence:
                    if ch != "-":
                        temp_filter_out_gaps += ch

                result_sequence = temp_filter_out_gaps
                """

                query_bound1 = int(seq_data["Hsp_query-from"])
                query_bound2 = int(seq_data["Hsp_query-to"])

                missing_left = query_bound1 - 1
                missing_right = query_seq_length - query_bound2
                # given this, can easily fill up with "N"s or "-"s

                hit_bound1 = int(seq_data['Hsp_hit-from'])
                hit_bound2 = int(seq_data['Hsp_hit-to'])
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
                        lower_bound = min(hit_bound1 + missing_left, hit_max)
                        upper_bound = hit_bound1 + 1
                        to_salvage = lower_bound - upper_bound + 1

                    string = ""
                    if 0 < to_salvage < 5:
                        arr = ncbi_get_gene_sequence(accession, lower_bound,
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
                        upper_bound = min(hit_bound2 + missing_right, hit_max)
                        to_salvage = upper_bound - lower_bound + 1
                    else:
                        lower_bound = hit_bound2 - 1
                        upper_bound = max(hit_bound2 - missing_right, 1)
                        to_salvage = lower_bound - upper_bound + 1

                    string = ""
                    if 5 > to_salvage > 0:
                        arr = ncbi_get_gene_sequence(accession, lower_bound, upper_bound, strand)

                        for ch in arr:
                            string += ch

                        string = string + "-"*(missing_right - to_salvage)
                    else:
                        string = string + "-"*missing_right

                    # print("missing right recovered: " + string + " iteration: " + exon_iteration['Iteration_iter-num'])

                    result_sequence = result_sequence + string

                transcript_file = self.create_transcript_file(ref_transcript_var, gene_name)

                fasta_heading = ">" + gene_name + " " + self.species_name + \
                                " reference_mrna:" + ref_transcript_var + \
                                " genome:" + accession + " " + \
                                ref_sequence_range
                transcript_file.write(fasta_heading + "\n")
                transcript_file.write(result_sequence + "\n")
                transcript_file.close()


class FullBlastXMLParser(BlastXMLParser):

    def __init__(self, file, save_dir, curr_species,
                 queries_to_genes_to_exons):
        super().__init__(file, save_dir, curr_species)
        self.genes_to_exons = queries_to_genes_to_exons[curr_species['del']]

    def parse_blast_xml(self):

        results_dict = file_xml_to_dictionary(self.xml_filepath)

        exon_iterations = results_dict['BlastOutput']['BlastOutput_iterations']
        for exon_iteration in exon_iterations:

            query_title = exon_iteration['Iteration_query-def'].split(" ")
            query_seq_length = int(exon_iteration['Iteration_query-len'])
            gene_name = query_title[0]

            mrna_section_no = 1
            while len(query_title[mrna_section_no]) < len("mRNA") or \
                    query_title[mrna_section_no][:len("mRNA")] != "mRNA":
                mrna_section_no += 1

            ref_transcript_var = query_title[mrna_section_no].split(":")[1]
            ref_sequence_range = query_title[mrna_section_no +
                                             SEQUENCE_INDICES_FROM_MRNA_TAG]

            hits_list = get_xml_list(exon_iteration['Iteration_hits'])
            if hits_list: # if not empty

                # get the top hit for the exon
                top_hit = hits_list[0]

                accession = top_hit['Hit_accession']
                # accession = top_hit['Hit_def'].split(" ")[0]
                hit_max = int(top_hit['Hit_len'])

                pg = PredictedGene()

                # if hit, there exists at least one hsp
                hsps = get_xml_list(top_hit["Hit_hsps"])

                for hsp in hsps:
                    pg.add_fragment(hsp)

                #result_sequence = seq_data["Hsp_hseq"]

                transcript_file = self.create_transcript_file(ref_transcript_var, gene_name)

                fasta_heading = ">" + gene_name + " " + self.species_name + \
                                " reference_mrna:" + ref_transcript_var + \
                                " genome:" + accession + " " + \
                                ref_sequence_range
                transcript_file.write(fasta_heading + "\n")
                #transcript_file.write(result_sequence + "\n")
                transcript_file.close()


class BlastFragment:

    def __init__(self, start, stop):
        self.seq = None
        self.start = int(start)
        self.stop = int(stop)

        self.next = None
        self.prev = None

    def set_seq(self, sequence):
        self.seq = sequence

    def merge_fragments(self, other):
        pass


class PredictedGene:

    def __init__(self):
        self.head = None

    def add_fragment(self, hsp):

        bf = BlastFragment(hsp["Hsp_query-from"], hsp["Hsp_query-to"])
        if self.head is None:
            self.head = bf
        else:
            curr_f = self.head






if __name__ == "__main__":
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    # parse_blast_xml(r"C:\Users\tonyx\Downloads\51JMZR2N016-Alignment.xml", r"C:\Users\tonyx\Downloads\xml_readtest", "altered2", "some weird one")





