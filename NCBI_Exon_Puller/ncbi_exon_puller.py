import os
import time
from typing import List, Dict, Tuple
from Bio import Entrez
import re

from Basic_Tools.lists_and_files import list_to_string


NCBI_CALL_ATTEMPTS = 5


class NcbiCallFailException(Exception):
    pass


def ncbi_gene_search(search_query: str, retmax=1000) -> Dict:
    """
    Searches NCBI's gene database with the given search_query

    :param search_query: properly formatted search query including gene name(s)
    and/or gene description(s)
    :param retmax: number of search results to return
    :return: the NCBI XML output, converted to dictionary format
    """
    # builds the search query
    # Entrez API, searchs the gene database returns matching results in
    # accession format
    handle = Entrez.esearch(db="gene", retmax=retmax, term=search_query, idtype="acc")
    # parse and return the XML output returned by Entrez
    return Entrez.read(handle)


def ncbi_get_gene_page_summaries(string_ids: str) -> List[Dict]:
    """
    Given a string of ids, gets gene summary pages for each.
    :param string_ids: a string of ids
    :return: a list of dictionaries each containing summary information for an
    id
    """
    # contained within a wrapper as the Entrez esummary call sometimes fails
    success, attempts = False, 0
    while not success and attempts < 5:
        try:
            esummary_handle = Entrez.esummary(db='gene', id=string_ids)
            success = True
        except:
            time.sleep(0.5)
            attempts += 1
            print("attempting to acquire gene summaries...")

    # prases the NCBI XML output and returns the list of summaries in dict. form
    return Entrez.read(esummary_handle)["DocumentSummarySet"]["DocumentSummary"]


def read_gene_table(text: str, seq_start: int, seq_stop: int) -> Dict:
    """
    Extracts exon information from a gene table imported from NCBI using the
    Entrez "efetch" feature.

    :param text: the gene table, in .txt format
    :param seq_start: the starting range of the gene
    :param seq_stop: the end range of the gene.
    :return: a dictionary, with each mrna transcript accession as a key and
    a list of exons based on a transcript as a value. Each exon is represented
    as a tuple of its endpoints.
    """
    # split the text file by line
    table = text.split("\n")

    # output dictionary
    dct = {}

    # the "line" in the table we are currently on
    line_no = 0
    while line_no < len(table):
        # iterates through lines until reaching dashes indicating the start of
        # an exon table
        while line_no < len(table) and (table[line_no] == "" or table[line_no][0] != "-"):
            line_no += 1

        # if we haven't reached the end of the table
        if line_no != len(table):
            # invariant: the name of the transcript that the exon table is based
            # off of will always be two lines up (subject to change)
            table_columns = table[line_no - 1].split("\t\t")

            # get the coding column -> if it doesn't exist, we don't do anything
            coding_no = 0
            while coding_no < len(table_columns) and \
                    table_columns[coding_no] != "Genomic Interval Coding":
                coding_no += 1

            table_heading = table[line_no - 2].split()
            if coding_no != len(table_columns) and len(table_heading) >= 5 and \
                    list_to_string(table_heading[:4], " ") == "Exon table for mRNA":

                transcript_name = table_heading[4]

                line_no += 1
                exons = []
                # iterate through the exon table
                while table[line_no] != "":
                    # split by double tabs (what separates the columns)
                    exon_line = table[line_no].split("\t\t")
                    # coding exon endpoints stored by dashes
                    endpoints = exon_line[coding_no].split("-")

                    # only include the "exon" if its endpoints are within the gene
                    if min(seq_start, seq_stop) <= min(int(endpoints[0]),
                                                       int(endpoints[1])) and \
                            max(int(endpoints[0]),
                                int(endpoints[1])) <= max(seq_start, seq_stop):
                        exons.append(endpoints)

                    line_no += 1

                dct[transcript_name] = exons

            else:
                # this is to skip over the "------" -> we hunt for the next one
                line_no += 1
        # continues for each transcript variant (which will have an exon table
    return dct


def ncbi_get_gene_features(acc: str, seq_start: int, seq_stop: int) -> Dict:
    """
    Returns the gene table information for a particular accession.

    :param acc: accession to fetch the gene table for
    :param seq_start: an endpoint of the gene
    :param seq_stop: an endpoint of the gene
    :return: a dictionary, with each mrna transcript accession in the gene table
    as a key and a list of exons based on a transcript as a value. Each exon is
    represented as a tuple of its endpoints.
    """
    # contained within a wrapper, necessasry as it is the API call that fails
    # the most often
    table_in_text = ""
    success, attempts = False, 0
    while not success and attempts < NCBI_CALL_ATTEMPTS:
        try:
            # fetch the gene table of a gene page, the only format is .txt
            gene_table_file = Entrez.efetch(db='gene', id=acc, rettype='gene_table',
                                            retmode="text")
            table_in_text = str(gene_table_file.read())
            success = True
        except:
            time.sleep(0.5)
            attempts += 1
            print("attempting to fetch table of gene features...")

    if attempts == NCBI_CALL_ATTEMPTS:
        raise NcbiCallFailException("Failure to fetch table of gene features")

    return read_gene_table(table_in_text, seq_start, seq_stop)


def ncbi_get_gene_sequence(chromosome_accession: str, seq_start: int,
                           seq_stop: int, strand: str) -> List[str]:
    """
    Fetches the segment of a chromosome corresponding to chromosome_accession
    from seq_start to seq_stop.

    :param chromosome_accession: accession of the chromosome to fetch from
    :param seq_start: endpoint of the gene
    :param seq_stop: endpoint of the gene
    :param strand: Takes on '1' or '2' for plus or mind strand, respectively
    :return: the gene segment in array format, where each index contains a
    single nucleotide.
    """
    # wrapper for potential error
    sequence_handle = ""
    success, attempts = False, 0
    while not success and attempts < NCBI_CALL_ATTEMPTS:
        try:
            # gets the fasta file in .txt format
            sequence_handle = Entrez.efetch(db='nuccore', id=chromosome_accession, rettype='fasta',
                                            retmode='text', seq_start=seq_start,
                                            seq_stop=seq_stop, strand=strand)
            success = True
        except:
            time.sleep(0.5)
            attempts += 1
            print("attempting to fetch sequence...")

    if attempts == NCBI_CALL_ATTEMPTS:
        raise NcbiCallFailException

    raw_text = str(sequence_handle.read())
    # splits the fasta output by line and omits the header
    raw_array = raw_text.split("\n")[1:]
    dna_array = []

    # converts the string to array for efficient exon indexing
    for line in raw_array:
        line = line.strip()
        for ch in line:
            dna_array.append(ch)

    return dna_array


def get_exon_sequence(exon: List[str], seq_start: int, seq_stop: int,
                      strand: str, dna_array: List[str]) -> str:
    """
    Given exon endpoints and an assembled DNA array of the gene, returns the
    exon sequence.

    :param exon: a tuple of the endpoints of the exon
    :param seq_start: start index of the gene / dna_array in chromosome context
    :param seq_stop: stop index of the gene / dna_array in chromosome context
    :param strand: plus or minus, can be either '1' or '2'
    :param dna_array: array representing the gene nucleotide sequence
    :return: the exon in string format
    """
    # everything is ascending from beginning to end of dna_array
    if strand == "1":
        start_index = int(exon[0]) - seq_start
        end_index = int(exon[1]) - seq_start
    # for minus strand, order is flipped and descends from beginning to end of
    # dna_array
    else:
        start_index = seq_start - int(exon[0])
        end_index = seq_start - int(exon[1])

    sequence = ""
    for i in range(start_index, end_index + 1):
        sequence += dna_array[i]

    return sequence


def write_exon_to_gene_file(transcript_filename: str, fasta_heading: str,
                            exon_sequence: str) -> None:
    """
    Appends a fasta sequence to a gene file.

    :param transcript_filename: name of file to append to
    :param fasta_heading: string
    :param exon_sequence: string
    :return:
    """
    # open for appending, so that all the exons of a gene can be appended to
    # the same file
    f = open(transcript_filename, "a")
    f.write(fasta_heading)
    f.write(exon_sequence + "\n")
    f.close()


def ncbi_exon_puller(search_query: str, gene_queries: List[str],
                     description_queries: List[str], taxon_folder: str,
                     gene_name: str):

    # search NCBI to get relevant accession numbers, then obtain summaries for
    # each accession
    search_results = ncbi_gene_search(search_query)
    ids = list_to_string(search_results["IdList"], ",")
    if ids != "":
        summaries = ncbi_get_gene_page_summaries(ids)
    else:
        summaries = []

    # summary_no is the current summary we are on
    summary_no = -1
    for summary in summaries:

        summary_no += 1

        # gene name, description and scientific name of the organism are
        # obtained in the summary to verify if the accession is relevant to our
        # search
        curr_gene_name = summary['Name']
        gene_description = summary['Description']
        scientific_name = summary["Organism"]['ScientificName']

        if (curr_gene_name.upper() in gene_queries or
            gene_description.upper() in description_queries) and \
                not os.path.isdir(os.path.join(taxon_folder, scientific_name)) and \
                len(summary['LocationHist']) != 0:
            # relevance is if 1) summary gene name or description matches one
            # of our queries, and 2) we have not yet encountered this species
            # for the particular gene (to avoid redundant visits)

            # genome accession that the gene's sequence is pulled from
            # print(summary)
            genome_accession = summary["LocationHist"][0]["ChrAccVer"]

            # make the folder that will contain transcripts for the species
            species_folder = os.path.join(taxon_folder, scientific_name)
            os.mkdir(species_folder)

            print("on species " + scientific_name)

            # obtain gene start-stop boundaries
            gene_start = int(summary["LocationHist"][0]["ChrStart"]) + 1
            gene_stop = int(summary["LocationHist"][0]["ChrStop"]) + 1

            # determines + or - strand
            strand = '1'
            if gene_start > gene_stop:
                strand = '2'

            # get the gene table (containing exon information)
            gene_features = ncbi_get_gene_features(search_results["IdList"][summary_no], gene_start, gene_stop)
            # get the sequence within the gene boundaries
            dna_array = ncbi_get_gene_sequence(genome_accession,
                                               gene_start, gene_stop,
                                               strand)
            time.sleep(0.5)
            # each transcript produces a different set of exons
            for transcript in gene_features.keys():

                transcript_filename = \
                    os.path.join(species_folder, transcript + ".fas")

                # iterates through the exons for each transcript
                exon_no = 1
                curr_exon_start = 1
                for exon in gene_features[transcript]:

                    curr_exon_length = abs(int(exon[1]) - int(exon[0]))

                    fasta_heading = ">" + gene_name + " " + \
                                    scientific_name + " mRNA:" + \
                                    transcript + " genome:" + \
                                    genome_accession + " " + \
                                    str(curr_exon_start) + "-" + \
                                    str(curr_exon_start + curr_exon_length) + \
                                    " exon number: " + str(exon_no) + "\n"

                    curr_exon_start += curr_exon_length + 1
                    # get the exon from the already extracted gene sequence
                    exon_sequence = get_exon_sequence(exon, gene_start,
                                                      gene_stop, strand,
                                                      dna_array)

                    write_exon_to_gene_file(transcript_filename,
                                            fasta_heading, exon_sequence)
                    exon_no += 1

if __name__ == "__main__":
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    gene_table_file = Entrez.efetch(db='gene', id="131764012,136129011,133094716", rettype='gene_table',
                                    retmode="text")
    table_in_text = str(gene_table_file.read())
    print(table_in_text)
