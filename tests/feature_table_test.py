import time

from Bio import Entrez

if __name__ == "__main__":
    Entrez.email = "xiaohan.xie@mail.utoronto.ca"

    success, attempts = False, 0
    while not success and attempts < 5:
        try:
            # gets the fasta file in .txt format
            sequence_handle = Entrez.efetch(db='nuccore', id="NM_174221.2",
                                            rettype='ft',
                                            retmode='text')
            success = True
        except:
            time.sleep(0.5)
            attempts += 1
            print("attempting to fetch sequence ")

    if attempts == 5:
        raise Exception("Failure to retrieve sequence with accession: ")

    output = str(sequence_handle.read()).split()
    start = None
    end = None
    for i in range(len(output)):
        if (i - 2 >= 0) and output[i] == "CDS" and output[i-1].isnumeric() and \
                output[i-2].isnumeric() and int(output[i-2]) < int(output[i-1]):
            start = output[i-2]
            end = output[i-1]

    if start is None:
        print("Error, no CDS for " + )
    else:


