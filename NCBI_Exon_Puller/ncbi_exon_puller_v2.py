import os
import time
from typing import Dict, List

from Bio import Entrez

from Basic_Tools.lists_and_files import list_to_string
from NCBI_Exon_Puller import ncbi_exon_puller as v1
from NCBI_Exon_Puller.Transcript import Transcript


def ncbi_get_gene_features_v2(ids: str, id_to_org: Dict[str, str]) -> Dict:

    table_in_text = ""
    success, attempts = False, 0
    while not success and attempts < v1.NCBI_CALL_ATTEMPTS:
        try:
            # fetch the gene table of a gene page, the only format is .txt
            gene_table_file = Entrez.efetch(db='gene', id=ids, rettype='gene_table',
                                            retmode="text")
            table_in_text = str(gene_table_file.read())
            success = True
        except Exception as e:
            print(e)
            time.sleep(0.5)
            attempts += 1
            print("attempting to fetch table of gene features...")

    if attempts == v1.NCBI_CALL_ATTEMPTS:
        raise v1.NcbiCallFailException("Failure to fetch table of gene features")

    return read_gene_table_v2(table_in_text, id_to_org)


def read_gene_table_v2(text: str, id_to_org: Dict[str, str]) -> Dict:

    table = text.split("\n")
    organisms = []
    line_no = 0

    while table[line_no].strip():
        split_line = table[line_no].split()
        if len(split_line) >= 3 and split_line[0] == "Gene" and split_line[1] == "ID:":
            curr_id = split_line[2][:-1]
            organisms.append(id_to_org[curr_id])
        line_no += 1

    orgs_to_transcripts = {}

    org_no = 0
    curr_org = ""
    while line_no < len(table):

        line_keywords = table[line_no].split()

        # new org
        if line_keywords and line_keywords[0] == "Reference":

            if curr_org in orgs_to_transcripts and not orgs_to_transcripts[curr_org]:
                orgs_to_transcripts.pop(curr_org)

            curr_org = organisms[org_no]
            orgs_to_transcripts[curr_org] = []
            org_no += 1

        # stumble on new exon table
        elif len(line_keywords) >= 5 and \
                list_to_string(line_keywords[:4], " ") == "Exon table for mRNA":

            line_no += 1
            table_columns = table[line_no].split("\t\t")

            # finds the respective columns
            exon_column_no = -1
            cds_column_no = -1
            for i in range(len(table_columns)):

                if table_columns[i] == "Genomic Interval Exon":
                    exon_column_no = i
                elif table_columns[i] == "Genomic Interval Coding":
                    cds_column_no = i

                if cds_column_no != -1 and exon_column_no != -1:
                    break

            acc_name = line_keywords[4]
            acc_quality = acc_name[:2]

            # if there exists both the exon and cds column
            if cds_column_no != -1 and exon_column_no != -1: # and acc_quality != "XM":

                transcript = Transcript(acc_name)

                line_no += 2
                table_data = table[line_no].split("\t\t")

                # keep going before hitting the empty line
                while len(table_data) > 1:

                    exon = table_data[exon_column_no].split("-")
                    cds = table_data[cds_column_no].split("-")

                    if transcript.check_valid_cds(exon, cds):
                        transcript.add_cds(exon, cds)
                    else:
                        transcript.add_to_exons_no_cds_space(exon)

                    line_no += 1
                    table_data = table[line_no].split("\t\t")

                orgs_to_transcripts[curr_org].append(transcript)

        line_no += 1

    return orgs_to_transcripts


def ncbi_get_transcript_sequences(accs: List[str]) -> Dict:

    accs_string = list_to_string(accs, ",")

    sequence_handle = ""
    success, attempts = False, 0
    while not success and attempts < v1.NCBI_CALL_ATTEMPTS:
        try:
            # gets the fasta file in .txt format
            sequence_handle = Entrez.efetch(db='nuccore', id=accs_string,
                                            rettype='fasta', retmode='text')
            success = True
        except:
            time.sleep(0.5)
            attempts += 1
            print("attempting to fetch sequence...")

    if attempts == v1.NCBI_CALL_ATTEMPTS:
        raise v1.NcbiCallFailException

    raw_arr = str(sequence_handle.read()).split("\n")

    line_no = 0
    transcript_to_sequence = {}
    while line_no < len(raw_arr):

        if raw_arr[line_no] and raw_arr[line_no][0] == ">":
            # get the name of the sequence, which is right after the ">"
            transcript_name = raw_arr[line_no].split()[0][1:]
            transcript_to_sequence[transcript_name] = []
            line_no += 1
            while raw_arr[line_no].strip() != "":
                for ch in raw_arr[line_no].strip():
                    transcript_to_sequence[transcript_name].append(ch)

                line_no += 1

        line_no += 1

    return transcript_to_sequence


def ncbi_exon_puller_v2(search_query: str, gene_queries: List[str],
                        description_queries: List[str], taxon_folder: str,
                        gene_name: str, exons_or_full):

    search_results = v1.ncbi_gene_search(search_query)
    ids_arr = search_results["IdList"]
    ids_string = list_to_string(ids_arr, ",")

    if ids_string == "":
        return None

    summaries = v1.ncbi_get_gene_page_summaries(ids_string)

    relevant_ids = []
    relevant_organisms = []
    summary_no = -1
    for summary in summaries:

        summary_no += 1

        curr_gene_name = summary['Name']
        gene_name_aliases = summary["OtherAliases"].split(",")
        gene_description = summary['Description']
        scientific_name = summary["Organism"]["ScientificName"]

        relevant = False
        if gene_name_aliases != "":

            alias_no = 0
            while alias_no < len(gene_name_aliases) and \
                    gene_name_aliases[alias_no].strip().upper() not in gene_queries:
                alias_no += 1

            if alias_no < len(gene_name_aliases):
                relevant = True

        if (curr_gene_name.upper() in gene_queries or
            gene_description.upper() in description_queries) and \
                not os.path.isdir(os.path.join(taxon_folder, scientific_name)) and \
                len(summary['LocationHist']) != 0:

            relevant = True

        if relevant:
            relevant_ids.append(ids_arr[summary_no])
            relevant_organisms.append(scientific_name)

    relevant_ids_string = list_to_string(relevant_ids, ",")

    id_to_org = {}
    for i in range(len(relevant_ids)):
        id_to_org[relevant_ids[i]] = relevant_organisms[i]

    orgs_to_transcripts = ncbi_get_gene_features_v2(relevant_ids_string,
                                                    id_to_org)

    transcript_list = []
    for org in orgs_to_transcripts:

        species_folder = os.path.join(taxon_folder, org)
        os.mkdir(species_folder)

        for transcript in orgs_to_transcripts[org]:
            transcript_list.append(transcript.accession)

    transcript_to_sequence = ncbi_get_transcript_sequences(transcript_list)

    for org in orgs_to_transcripts:
        for transcript in orgs_to_transcripts[org]:
            if transcript.cds_end <= len(transcript_to_sequence[transcript.accession]):
                transcript_filename = \
                    os.path.join(taxon_folder, org,
                                 str(transcript.get_cds_length()) + "_" +
                                 transcript.accession + ".fas")

                transcript_file = open(transcript_filename, "w")

                if exons_or_full == "exons":
                    cds_no = 1
                    for cds in transcript.cds:

                        fasta_heading = ">" + gene_name + " " + org + " mRNA:" + \
                                        transcript.accession + " genome:na " + \
                                        str(cds[0] - transcript.cds_start + 1) + \
                                        "-" + \
                                        str(cds[1] - transcript.cds_start + 1) + \
                                        " exon number: " + str(cds_no) + "\n"

                        cds_no += 1
                        cds_sequence = ""

                        for i in range(cds[0]-1, cds[1]):
                            cds_sequence += transcript_to_sequence[transcript.accession][i]

                        transcript_file.write(fasta_heading + cds_sequence + "\n")
                else:
                    fasta_heading = ">" + gene_name + " " + org + " mRNA:" + \
                                    transcript.accession + " genome:na " + \
                                    "1-" + \
                                    str(transcript.cds_end - transcript.cds_start + 1) + \
                                    "\n"

                    cds_sequence = ""
                    for i in range(transcript.cds_start-1, transcript.cds_end):
                        cds_sequence += transcript_to_sequence[transcript.accession][i]

                    transcript_file.write(fasta_heading + cds_sequence + "\n")

                transcript_file.close()


if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"
    #gene_table_file = Entrez.efetch(db='gene', id="105573713", rettype='gene_table',
    #                                retmode="text")
    #table_in_text = str(gene_table_file.read())
    #print(table_in_text)

    handle = Entrez.efetch(db="taxonomy", id='elasmobranchii')
    print(Entrez.read(handle))


    #for summary in summaries:
    #    print(summary)
    #ncbi_exon_puller_v2("cetacea[orgn] rho[gene]", ["RHO"], [], "hello", "")





