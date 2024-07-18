import os
from abc import abstractmethod

from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from Bio import Entrez
from Basic_Tools.xml_extraction import file_xml_to_dictionary

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


class BlastXMLParser:

    def __init__(self, file, save_dir, taxon_name, species_name):

        self.xml_filepath = file
        self.save_dir = save_dir

        self.taxon_name = taxon_name
        self.species_name = species_name

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

            if isinstance(exon_iteration['Iteration_hits'], list) or \
                    isinstance(exon_iteration['Iteration_hits'], dict): # if not empty

                # get the top hit for the exon
                if isinstance(exon_iteration['Iteration_hits'], list):
                    top_hit = exon_iteration['Iteration_hits'][0]
                else:
                    top_hit = exon_iteration['Iteration_hits']['Hit']
                # get the top sequence for the top hit
                if isinstance(top_hit["Hit_hsps"], list):
                    seq_data = top_hit["Hit_hsps"][0]
                else:
                    seq_data = top_hit["Hit_hsps"]["Hsp"]

                result_sequence = seq_data["Hsp_hseq"]

                temp_filter_out_gaps = ""
                for ch in result_sequence:
                    if ch != "-":
                        temp_filter_out_gaps += ch

                result_sequence = temp_filter_out_gaps

                query_bound1 = int(seq_data["Hsp_query-from"])
                query_bound2 = int(seq_data["Hsp_query-to"])

                missing_left = query_bound1 - 1
                missing_right = query_seq_length - query_bound2
                # given this, can easily fill up with "N"s or "-"s

                # accession = top_hit['Hit_def'].split(" ")[0]
                accession = top_hit['Hit_accession']
                #print("accessing accession: " + accession)

                hit_max = int(top_hit['Hit_len'])

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

    def parse_blast_xml(self):
        pass

        # basically, for each iteration:
        # we go to the top hit, and go through all of its hit hsps


if __name__ == "__main__":
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    # parse_blast_xml(r"C:\Users\tonyx\Downloads\51JMZR2N016-Alignment.xml", r"C:\Users\tonyx\Downloads\xml_readtest", "altered2", "some weird one")





