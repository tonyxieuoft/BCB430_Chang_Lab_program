import time
from os import path

from Bio import Entrez

def string_empty(str):
    return str != ""

if __name__ == "__main__":

    Entrez.email = input("Please enter your email. This is required to access "
                         "NCBI Genbank via the API.")

    accession_filepath = input("Please enter the complete filepath to a file "
                                "containing the accessions. The file should be "
                                "in .txt format, with only one accession per "
                                "line.")

    accession_dirpath = input("Please enter the complete path to a "
                               "directory you wish the file to the outputted "
                               "to. ")
    while not path.isdir(accession_dirpath):
        accession_dirpath = input("Invalid directory path, please enter again.")

    out_path_name = "cds_output"

    file_num = 1
    out_path = path.join(accession_dirpath, out_path_name + ".fas")
    while path.isfile(out_path):

        out_path = path.join(accession_dirpath, out_path_name +
                             " (" + str(file_num) + ").fas")
        file_num += 1

    output_f = open(out_path, "w")

    accessions = filter(string_empty, open(accession_filepath).read().split())

    for acc in accessions:

        # get taxid
        summary_handle = Entrez.esummary(db='nuccore', id=acc)
        summaries = Entrez.parse(summary_handle)
        taxa = "0"
        for summary in summaries:
            taxa = int(summary["TaxId"])

        # get name from taxid
        tax_fetch = Entrez.efetch(db="taxonomy", id=taxa)
        tax = Entrez.parse(tax_fetch)
        first, last = None, None
        for t in tax:
            first, last = t["ScientificName"].split()[:2]

        # get CDS range
        success, attempts = False, 0
        while not success and attempts < 5:
            try:
                # gets the fasta file in .txt format
                sequence_handle = Entrez.efetch(db='nuccore', id=acc,
                                                rettype='ft',
                                                retmode='text')
                success = True
            except:
                time.sleep(0.5)
                attempts += 1
                print("attempting to fetch sequence ")

        if attempts == 5:
            raise Exception("Failure to get feature table for accession: " + acc)

        output = str(sequence_handle.read()).split()
        start = None
        end = None
        for i in range(len(output)):
            if (i - 2 >= 0) and output[i] == "CDS" and output[i-1].isnumeric() and \
                    output[i-2].isnumeric() and int(output[i-2]) < int(output[i-1]):
                start = output[i-2]
                end = output[i-1]

        # get sequence
        if start is None:
            print("Issue: no CDS for accession " + acc)
        else:
            success, attempts = False, 0
            while not success and attempts < 5:
                try:
                    # gets the fasta file in .txt format
                    sequence_handle = Entrez.efetch(db='nuccore', id=acc,
                                                    rettype='fasta',
                                                    retmode='text',
                                                    seq_start = start,
                                                    seq_stop = end)
                    success = True
                except:
                    time.sleep(0.5)
                    attempts += 1
                    print("attempting to fetch sequence " + acc)

            if attempts == 5:
                raise Exception("Failure to retrieve sequence with accession: " + acc)

            output_f.write(">" + first[0] + last[:3] + "\n" + "\n".join(sequence_handle.read().split("\n")[1:]))
            print("success for accession: " + acc)




# testing
# 1) xiaohan.xie@mail.utoronto.ca
# 2) C:\Users\tonyx\Downloads\sample_acc.txt
# 3) C:\Users\tonyx\Downloads
