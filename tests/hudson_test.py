from Bio import Entrez
from Bio.Entrez.Parser import IntegerElement

if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@mail.utoronto.ca"


    summary_handle = Entrez.esummary(db='nuccore', id="XM_057528069.1")
    summaries = Entrez.parse(summary_handle)
    for summary in summaries:
        taxa = int(summary["TaxId"])

    tax_fetch = Entrez.efetch(db="taxonomy", id=taxa)
    tax = Entrez.parse(tax_fetch)
    for t in tax:
        first, last = t["ScientificName"].split()

    print(first[0])
    print(last[:3])


    sequence_handle = Entrez.efetch(db='nuccore', id="XM_057528069.1",
                                    rettype='ft',
                                    retmode='text')
    output = str(sequence_handle.read()).split()
    print(output)
