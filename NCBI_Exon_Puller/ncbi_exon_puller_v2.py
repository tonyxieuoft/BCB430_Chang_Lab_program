import os
import time
from typing import Dict, List

from Bio import Entrez

from Basic_Tools.lists_and_files import list_to_string
from NCBI_Exon_Puller import ncbi_exon_puller as v1
from NCBI_Exon_Puller.Transcript import Transcript
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence


splice_file_output = "/crun/tony.xie/Downloads/reference_splice_site_output.csv"


def ncbi_get_gene_features_v2(ids: str, id_to_org: Dict[str, str], gene_name) -> Dict:

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

    return read_gene_table_v2(table_in_text, id_to_org, gene_name)


def read_gene_table_v2(text: str, id_to_org: Dict[str, str], gene_name: str) -> Dict:

    #print(text)
    # text is the entire table
    table = text.split("\n")

    organisms = []
    line_no = 0

    # while we haven't reached an empty line
    while table[line_no].strip():
        split_line = table[line_no].split()
        if len(split_line) >= 3 and split_line[0] == "Gene" and split_line[1] == "ID:":
            curr_id = split_line[2][:-1]
            organisms.append(id_to_org[curr_id])
        line_no += 1

    orgs_to_transcripts = {}

    org_no = 0
    curr_org = ""
    genomic_id = ""
    while line_no < len(table):

        # only relevant for storing the current genomic ID

        line_keywords = table[line_no].split()

        # new org
        if line_keywords and line_keywords[0] == "Reference":

            # ======= To get splice sites ============
            for keyword in line_keywords:
                if len(keyword) >= 3 and keyword[0] == "N" and keyword[2] == "_":

                    genomic_id = keyword

                    break


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

            ### ============ Get Splice Site ========================
            relevant_transcript_ids = ['XM_059956423.1', 'XM_041173625.1', 'XM_059945786.1', 'XM_041207015.1', 'XM_059950643.1', 'XM_041195929.1', 'XM_059955473.1', 'XM_041213792.1', 'XM_059946101.1', 'XM_041176112.1', 'XM_059975050.1', 'XM_041178838.1', 'XM_059947580.1', 'XM_041210456.1', 'XM_059952777.1', 'XM_041174680.1', 'XM_059963252.1', 'XM_041199335.1', 'XM_059946625.1', 'XM_041178670.1', 'XM_059948319.1', 'XM_041216848.1', 'XM_052017529.1', 'XM_069936896.1', 'XM_059944490.1', 'XM_055647860.1', 'XM_033036884.1', 'XM_063067743.1', 'XM_041190931.1', 'XM_059990205.1', 'XM_041193610.1', 'XM_059967774.1', 'XM_041211548.1', 'XM_059985789.1', 'XM_059956044.1', 'XM_041180364.1', 'XM_059986830.1', 'XM_041203252.1', 'XM_059963352.1', 'XM_041194575.1', 'XM_059963252.1', 'XM_041199335.1', 'XM_059952428.1', 'XM_041206315.1', 'XM_059989449.1', 'XM_041173738.1', 'XM_059950157.1', 'XM_059971408.1', 'XM_041196729.1', 'XM_059993484.1', 'XM_041192083.1', 'XM_059952479.1', 'XM_033044102.1', 'XM_059961194.1', 'XM_041215304.1', 'XM_063064065.1', 'XM_041207809.1']

            # if there exists both the exon and cds column
            if cds_column_no != -1 and exon_column_no != -1 and acc_name in relevant_transcript_ids: # and acc_quality != "XM":

                # ======= For splice site ==========
                real_range = []
                # ==================================

                transcript = Transcript(acc_name)

                line_no += 2
                table_data = table[line_no].split("\t\t")

                # keep going before hitting the empty line
                while len(table_data) > 1:

                    exon = table_data[exon_column_no].split("-")
                    cds = table_data[cds_column_no].split("-")

                    if transcript.check_valid_cds(exon, cds):
                        transcript.add_cds(exon, cds)
                        real_range.append([int(cds[0]), int(cds[1])])

                    else:
                        transcript.add_to_exons_no_cds_space(exon)

                    line_no += 1
                    table_data = table[line_no].split("\t\t")

                # =============== For splice site ============
                # =============================================

                junction_position = 0
                strand = ""

                for i in range(len(real_range)):

                    get_start = True
                    get_end = True

                    if i == 0:
                        get_start = False
                        if real_range[i][0] < real_range[i][1]:
                            strand = "1"
                        else:
                            strand = "2"

                    elif i == len(real_range) - 1:
                        get_end = False

                    if strand == "1":
                        arr = ncbi_get_gene_sequence(genomic_id, real_range[i][0],
                                                     real_range[i][1], strand)
                    else:
                        arr = ncbi_get_gene_sequence(genomic_id, real_range[i][1],
                                                     real_range[i][0], strand)
                    string = ""
                    for ch in arr:
                        string += ch

                    if get_start:
                        if strand == "1":
                            arr = ncbi_get_gene_sequence(genomic_id, real_range[i][0] - 2,
                                                        real_range[i][0] - 1, strand)
                        else:
                            arr = ncbi_get_gene_sequence(genomic_id, real_range[i][0] + 2,
                                                         real_range[i][0] + 1, strand)

                        splice_seq = ""
                        for ch in arr:
                            splice_seq += ch

                        splice_f = open(splice_file_output, "a")
                        print(splice_seq + ",end," + curr_org + "," + str(junction_position) + "," + gene_name + "\n")
                        splice_f.close()

                    if get_end:

                        # very important for updating junction position
                        junction_position += abs(real_range[i][0] - real_range[i][1]) + 1

                        if strand == "1":
                            arr = ncbi_get_gene_sequence(genomic_id, real_range[i][1] + 1,
                                                         real_range[i][1] + 2, strand)
                        else:
                            arr = ncbi_get_gene_sequence(genomic_id, real_range[i][1] - 1,
                                                         real_range[i][1] - 2, strand)

                        splice_seq = ""
                        for ch in arr:
                            splice_seq += ch

                        splice_f = open(splice_file_output, "a")
                        splice_f.write(splice_seq + ",start," + curr_org + "," + str(junction_position) + "," + gene_name + "\n")
                        splice_f.close()

                # =======================

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

    # get the initial search results (based on initial query)
    search_results = v1.ncbi_gene_search(search_query)

    # then, get the ids that come out from it
    ids_arr = search_results["IdList"]
    ids_string = list_to_string(ids_arr, ",")

    if ids_string == "":
        return None

    # get summaries
    summaries = v1.ncbi_get_gene_page_summaries(ids_string)

    relevant_ids = []
    relevant_organisms = []
    summary_no = -1
    # look through the summaries
    for summary in summaries:

        summary_no += 1

        curr_gene_name = summary['Name']
        gene_name_aliases = summary["OtherAliases"].split(",")
        gene_description = summary['Description']
        scientific_name = summary["Organism"]["ScientificName"]

        # relevant is true if the ID (obtained from basic search) is relevant
        relevant = False

        # look at the aliases
        if gene_name_aliases != "":

            alias_no = 0
            # go through all of the aliases
            while alias_no < len(gene_name_aliases) and \
                    gene_name_aliases[alias_no].strip().upper() not in gene_queries:
                alias_no += 1

            # if one of the aliases matched the name of our gene of interest
            if alias_no < len(gene_name_aliases):
                relevant = True

        # if (the ID's gene name or description matches) and we haven't come across
        # it yet (reflected in the directory checking)
        if (curr_gene_name.upper() in gene_queries or
            gene_description.upper() in description_queries) and \
                not os.path.isdir(os.path.join(taxon_folder, scientific_name)) and \
                len(summary['LocationHist']) != 0:

            relevant = True

        # add it to our relevant list
        if relevant:
            relevant_ids.append(ids_arr[summary_no])
            relevant_organisms.append(scientific_name)

    # relevant IDs concatenate (to prepare for search)
    relevant_ids_string = list_to_string(relevant_ids, ",")

    id_to_org = {}
    for i in range(len(relevant_ids)):
        # basic search ID to organism mapping
        id_to_org[relevant_ids[i]] = relevant_organisms[i]

    # get the features
    orgs_to_transcripts = ncbi_get_gene_features_v2(relevant_ids_string,
                                                    id_to_org, gene_name)

    # create a list
    transcript_list = []
    for org in orgs_to_transcripts:

        # prepare the species directories
        species_folder = os.path.join(taxon_folder, org)
        os.mkdir(species_folder)

        for transcript in orgs_to_transcripts[org]:
            transcript_list.append(transcript.accession)

    # and then get the sequences themselves (all at once, using RNA ids)
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

    #handle = Entrez.efetch(db="taxonomy", id='elasmobranchii')
    #print(Entrez.read(handle))


    #for summary in summaries:
    #    print(summary)
    ncbi_exon_puller_v2("cetacea[orgn] rho[gene]",
                        ["RHO"],
                        [],
                        "/Users/tonyx/Downloads/taxon_dir",
                        "hello",
                        "exons")





