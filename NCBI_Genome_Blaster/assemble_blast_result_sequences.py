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

        transcript_filepath = os.path.join(species_folder, "Reference_" +
                                           ref_transcript_var + ".fas")
        return open(transcript_filepath, "a")

    @abstractmethod
    def parse_blast_xml(self):
        pass


class ExonBlastXMLParser(BlastXMLParser):

    def __init__(self, file, save_dir, curr_species, on_server):
        super().__init__(file, save_dir, curr_species)
        self.on_server = on_server

    def parse_blast_xml(self):

        results_dict = file_xml_to_dictionary(self.xml_filepath)

        exon_iterations = results_dict['BlastOutput']['BlastOutput_iterations']
        for exon_iteration in exon_iterations:

            query_title = exon_iteration['Iteration_query-def'].split(" ")
            query_seq_length = int(exon_iteration['Iteration_query-len'])
            gene_name = query_title[0]

            mrna_section_no = 1
            while (len(query_title[mrna_section_no]) < len("mRNA") or \
                    query_title[mrna_section_no][:len("mRNA")] != "mRNA") and \
                    (len(query_title[mrna_section_no]) < len("reference") or
                     query_title[mrna_section_no][:len("reference")] != "reference"):
                mrna_section_no += 1

            ref_transcript_var = query_title[mrna_section_no].split(":")[1]
            ref_sequence_range = query_title[mrna_section_no +
                                             SEQUENCE_INDICES_FROM_MRNA_TAG]

            hits_list = get_xml_list(exon_iteration['Iteration_hits'])
            if hits_list: # if not empty

                # get the top hit for the exon
                top_hit = hits_list[0]

                if self.on_server:
                    # TODO this only works for NCBI, be careful
                    accession = top_hit['Hit_def'].split(" ")[0]
                else:
                    accession = top_hit['Hit_accession']

                hit_max = int(top_hit['Hit_len'])

                # if hit, there exists at least one hsp
                seq_data = get_xml_list(top_hit["Hit_hsps"])[0]

                qseq_num_gaps = 0
                for ch in seq_data["Hsp_qseq"]:
                    if ch == "-":
                        qseq_num_gaps += 1

                hseq_num_gaps = 0
                for ch in seq_data["Hsp_hseq"]:
                    if ch == "-":
                        hseq_num_gaps += 1

                if qseq_num_gaps == hseq_num_gaps and qseq_num_gaps > 0:
                    result_sequence = ""
                    for ch in seq_data["Hsp_hseq"]:
                        if ch != "-":
                            result_sequence += ch

                else:
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
                    if 0 < to_salvage < 10: # TODO change it back to 0-5 later
                        print("getting gene sequence")
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
                    if 10 > to_salvage > 0: # change it back later
                        print("getting gene seqnece")
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

        exon_iterations = get_xml_list(results_dict['BlastOutput']['BlastOutput_iterations'])
        for exon_iteration in exon_iterations:

            query_title = exon_iteration['Iteration_query-def'].split(" ")
            query_seq_length = int(exon_iteration['Iteration_query-len'])
            gene_name = query_title[0]
            print(gene_name)

            mrna_section_no = 1
            while (len(query_title[mrna_section_no]) < len("mRNA") or \
                    query_title[mrna_section_no][:len("mRNA")] != "mRNA") and \
                    (len(query_title[mrna_section_no]) < len("reference") or
                     query_title[mrna_section_no][:len("reference")] != "reference"):
                mrna_section_no += 1

            ref_transcript_var = query_title[mrna_section_no].split(":")[1]
            ref_sequence_range = query_title[mrna_section_no +
                                             SEQUENCE_INDICES_FROM_MRNA_TAG]

            hits_list = get_xml_list(exon_iteration['Iteration_hits'])
            if hits_list: # if not empty

                # get the top hit for the exon

                pg = PredictedGene(self.genes_to_exons[gene_name])
                accession = hits_list[0]['Hit_accession']
                # accession = top_hit['Hit_def'].split(" ")[0]
                """
                all_hsps = []
                for hit in hits_list:
                    hsps = get_xml_list(hit["Hit_hsps"])
                    all_hsps += hsps

                for i in range(len(all_hsps)):
                    min_index = -1
                    min_value = float('inf')
                    for j in range(i, len(all_hsps)):
                        if float(all_hsps[j]["Hsp_evalue"]) < min_value:
                            min_index = j
                            min_value = float(all_hsps[j]["Hsp_evalue"])
                    all_hsps[i], all_hsps[min_index] = all_hsps[min_index], all_hsps[i]
                    
                """
                hit_counter = 0
                while hit_counter < len(hits_list): # hit_counter < 2 and

                    hsps = get_xml_list(hits_list[hit_counter]["Hit_hsps"])
                    for hsp in hsps:
                        pg.add_fragment(hsp)
                    hit_counter += 1

                result_sequence = pg.construct_full_seq()

                #result_sequence = seq_data["Hsp_hseq"]

                transcript_file = self.create_transcript_file(ref_transcript_var, gene_name)

                fasta_heading = ">" + gene_name + " " + self.species_name + \
                                " reference_mrna:" + ref_transcript_var + \
                                " genome:" + accession + " " + \
                                ref_sequence_range
                transcript_file.write(fasta_heading + "\n")
                transcript_file.write(result_sequence + "\n")
                transcript_file.close()


class BlastFragment:

    def __init__(self, start, stop, q_seq, h_seq):

        self.start = int(start)
        self.stop = int(stop)

        self.next = None
        self.prev = None

        self.seq = []

        num_q_gaps = 0
        for ch in q_seq:
            if ch == "-":
                num_q_gaps += 1

        num_h_gaps = 0
        for ch in h_seq:
            if ch == "-":
                num_h_gaps += 1

        if num_q_gaps == num_h_gaps:
            for ch in h_seq:
                if ch != "-":
                    self.seq.append(ch)
        else:
            i = 0
            while i < len(q_seq):
                curr = h_seq[i]
                i += 1

                while i < len(q_seq) and q_seq[i] == "-":
                    curr += h_seq[i]
                    i += 1

                self.seq.append(curr)


class PredictedGene:

    def __init__(self, ref_exons):
        self.head = None
        self.ref_exons = ref_exons

    def add_fragment(self, hsp):

        bf = BlastFragment(hsp["Hsp_query-from"], hsp["Hsp_query-to"],
                           hsp["Hsp_qseq"], hsp["Hsp_hseq"])

        if self.head is None:
            self.head = bf
        else:
            curr_f = self.head
            while curr_f.next is not None and bf.stop > curr_f.stop:
                curr_f = curr_f.next

            # if we've reached the end, treat curr as behind
            if curr_f.next is None and bf.stop > curr_f.stop:

                if bf.start <= curr_f.start:
                    return None

                if (curr_f.stop - bf.start + 1) / (bf.stop - bf.start + 1) > 0.25:
                    print("curr_f.stop: " + str(curr_f.stop) + " bf.start: " + str(bf.start) + " bf.stop: " + str(bf.stop))
                    return None

                curr_f.next = bf
                bf.prev = curr_f

            # we're at the beginning, and the only curr is forward
            elif curr_f == self.head:

                if bf.start >= curr_f.start:
                    return None

                if (bf.stop - curr_f.start + 1) / (bf.stop - bf.start + 1) > 0.25:
                    print("curr_f.start: " + str(curr_f.start) + " bf.start: " + str(bf.start) + " bf.stop: " + str(bf.stop))
                    return None

                curr_f.prev = bf
                bf.next = curr_f
                self.head = bf

            # we're in the middle somewhere
            else:  # bf.stop <= curr_f.stop:
                if bf.start >= curr_f.start:  # encompassed inside
                    return None

                if bf.start <= curr_f.prev.start: # encompassed inside
                    return None

                if (curr_f.prev.stop - bf.start + 1) / (bf.stop - bf.start + 1) > 0.25:
                    print("connected to next: start:" + str(curr_f.start) + " stop: " + str(curr_f.stop))
                    print("curr_f.prev.stop: " + str(curr_f.prev.stop) + " bf.start: " + str(bf.start) + " bf.stop: " + str(bf.stop))
                    return None

                if (bf.stop - curr_f.start + 1) / (bf.stop - bf.start + 1) > 0.25:
                    print("curr_f.start (mid): " + str(curr_f.start) + " bf.start: " + str(bf.start) + " bf.stop: " + str(bf.stop))
                    return None

                bf.next = curr_f
                bf.prev = curr_f.prev
                curr_f.prev.next = bf
                curr_f.prev = bf

            # now, adjust for overlap
            if bf.next is not None and bf.stop >= bf.next.start:
                self._resolve_overlap(bf, bf.next)
            if bf.prev is not None and bf.prev.stop >= bf.start:
                self._resolve_overlap(bf.prev, bf)

    def _resolve_overlap(self, bf1, bf2):

        junction_start = self._look_for_exon_boundary(bf2.start, bf1.stop)

        if junction_start is None:
            new_seq = bf1.seq + ["-"]*(bf2.stop - bf1.stop)

            prio = ""
            for ch in bf1.seq[(bf2.start - bf1.start):]:
                if len(ch) > 1 or ch == "-" or ch == "N":
                    prio = "bf2"
                    break

            for ch in bf2.seq[:(bf1.stop - bf2.start + 1)]:
                if len(ch) > 1 or ch == "-" or ch == "N":
                    if prio == "bf2":
                        prio = ""
                    else:
                        prio = "bf1"
                    break

            rel_start = bf2.start - bf1.start
            if prio == "bf2":
                for i in range(len(bf2.seq)):
                    new_seq[i + rel_start] = bf2.seq[i]
            elif prio == "bf1":
                for i in range(bf1.stop - bf2.start + 1, len(bf2.seq)):
                    new_seq[i + rel_start] = bf2.seq[i]
            else:
                for i in range(len(bf2.seq)):
                    if new_seq[i + rel_start] != "-" and new_seq[i + rel_start] != bf2.seq[i]:
                        print(new_seq[i + rel_start] + " " + bf2.seq[i] + " index: " + str(i))
                        new_seq[i + rel_start] = "N"
                    else:
                        new_seq[i + rel_start] = bf2.seq[i]

        else:
            new_seq = bf1.seq[:(junction_start - bf1.start + 1)] + \
                      bf2.seq[(junction_start + 1 - bf2.start):]

        bf1.seq = new_seq
        bf1.stop = bf2.stop
        bf1.next = bf2.next
        if bf2.next is not None:
            bf2.next.prev = bf1

        """
        junction_start = self._look_for_exon_boundary(bf2.start, bf1.stop)
        if junction_start is not None:

            bf1_change = bf1.stop - junction_start
            bf1.stop = junction_start
            bf1.seq = bf1.seq[:-bf1_change]

            bf2_change = junction_start + 1 - bf2.start
            bf2.start = junction_start + 1
            bf2.seq = bf2.seq[bf2_change:]

        else:
            # implement exon junction checking here later
            bf1_overlap_seq = bf1.seq[(bf2.start - bf1.start):]
            bf2_overlap_seq = bf2.seq[:(bf1.stop - bf2.start + 1)]
            averaged_overlap = ""
            if len(bf1_overlap_seq) != len(bf2_overlap_seq):
                print("overlap is wrong...")
            else:
                for i in range(len(bf1_overlap_seq)):
                    if bf1_overlap_seq[i] == bf2_overlap_seq[i]:
                        averaged_overlap += bf1_overlap_seq[i]
                    else:
                        averaged_overlap += "N"

            forced_junction = (bf2.start + bf1.stop) // 2
            ao_split1 = averaged_overlap[:((len(averaged_overlap) + 1) // 2)]
            ao_split2 = averaged_overlap[((len(averaged_overlap) + 1) // 2):]

            bf1.seq = bf1.seq[:bf2.start - bf1.start] + ao_split1
            bf1.stop = forced_junction

            bf2.seq = ao_split2 + bf2.seq[(bf1.stop - bf2.start + 1):]
            bf2.start = forced_junction + 1
            
        """

    def _look_for_exon_boundary(self, overlap_start, overlap_end):

        for exon in self.ref_exons:
            if overlap_start - 1 <= exon[1] <= overlap_end:
                return exon[1]

        return None

    def construct_full_seq(self):

        result_seq = "-" * (self.head.start - 1)

        curr = self.head
        while curr.next is not None:
            print(str(curr.start) + " " + str(curr.stop))
            for ch in curr.seq:
                result_seq += ch
            result_seq += (curr.next.start - curr.stop - 1) * "-"
            curr = curr.next

        print(str(curr.start) + " " + str(curr.stop))
        for ch in curr.seq:
            result_seq += ch
        result_seq += (self.ref_exons[-1][1] - curr.stop) * "-"
        return result_seq


class PredictedGeneV2(PredictedGene):

    def __init__(self, ref_exons):
        super().__init__(ref_exons)
        self.gene_arr = ["-"]*(ref_exons[-1][1] + 1)
        self.gap_arr = [""]*(ref_exons[-1][1] + 1)

    def add_fragment(self, hsp):
        seq = gap_corrector(hsp["Hsp_qseq"], hsp["Hsp_hseq"])
        start = int(hsp["Hsp_query-from"])

        if seq is not None:
            not_occupied = 0
            for i in range(len(seq)):
                if self.gene_arr[i + start] == "-":
                    not_occupied += 1

            if not_occupied / len(seq) > 0.75:
                for i in range(len(seq)):
                    if self.gene_arr[i + start] != "-" \
                            and self.gene_arr[i + start] != seq[i]:
                        self.gene_arr[i + start] = "N"
                    else:
                        self.gene_arr[i + start] = seq[i]

    def construct_full_seq(self):
        result = ""
        for ch in self.gene_arr[1:]:
            result += ch

        return result


if __name__ == "__main__":
    pass
    # xiaohan.xie@mail.utoronto.ca
    # C:\Users\tonyx\Downloads\main_test4
    # C:\Users\tonyx\Downloads\main_test3\genes.txt
    # C:\Users\tonyx\Downloads\main_test3\taxa.txt
    # C:\Users\tonyx\Downloads\main_test4\NCBI_exon_pull_results
    # parse_blast_xml(r"C:\Users\tonyx\Downloads\51JMZR2N016-Alignment.xml", r"C:\Users\tonyx\Downloads\xml_readtest", "altered2", "some weird one")





