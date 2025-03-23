import os
from Bio import Entrez

GENOMIC_ACCESSION = 0
TYPE = 2
BOUND1 = 3
BOUND2 = 4
STRAND = 6
ANNOTATIONS = 8


class GeneModel:

    # just a convenient way to store info, no real functions in itself until later
    def __init__(self, accession, score, strand, gene_bound1, gene_bound2):

        self.score = score
        self.accession = accession
        self.gene_range = (int(gene_bound1), int(gene_bound2))

        if strand == "+":
            self.strand = "1"
        else:
            self.strand = "2"

        self.cds_ranges = []

    def add_to_range(self, range1, range2):
        self.cds_ranges.append((int(range1), int(range2)))


def get_info_from_annotation(annotation):

    gene_name = ""
    score = 0

    for entry in annotation:
        entry_type = entry.split("=")[0]
        if entry_type == "ID":
            # in the form ID=genename_R1,   need to get rid of the R1
            gene_name = entry.split("=")[1].split("_")[0]
        if entry_type == "score":
            score = int(entry.split("=")[1])

    return gene_name, score

def process_gemoma_results(annotation_path, species):

    """
    The whole file needs to be read and processed first, since there are multiple

    Gene model dictionary is in the form {gene_name: GeneModel}
    """

    gene_model_dict = {}

    annot_f = open(annotation_path, "r")

    # read the first few lines
    line = annot_f.readline()
    line = annot_f.readline()
    line = annot_f.readline()

    while line.strip() != "":

        # at the beginning, the invariant is that we're at a gene line or ## genomic delimeter

        # skip two if delimiter, one for gene line
        if line[0] == "#":
            line = annot_f.readline()
        line = annot_f.readline()

        # now, this line must be the mRNA line
        gff_arr = line.split()
        if gff_arr[TYPE] != "mRNA":
            print("something is wrong")

        # get the useful info out of the annotations section for the mRNA
        annotations_arr = gff_arr[ANNOTATIONS].split(";")
        gene_name, score = get_info_from_annotation(annotations_arr)

        # if gene name is not in gene_model_dict, we haven't encountered it yet
        # if it exists in the dict, it means a dfferent version is currently in it and we compare the score
        if gene_name not in gene_model_dict or gene_model_dict[gene_name].score < score:

            # create the new model
            new_model = GeneModel(annotations_arr[GENOMIC_ACCESSION], score,
                                  annotations_arr[STRAND],
                                  annotations_arr[BOUND1], annotations_arr[BOUND2])

            gene_model_dict[gene_name] = new_model

            line = annot_f.readline()
            while line[0] != "#" and line.split()[TYPE] != "gene":
                gff_arr = line.split()
                new_model.add_to_range(gff_arr[BOUND1], gff_arr[BOUND2])

        else:
            line = annot_f.readline()
            while line.strip() != "" and line[0] != "#" and line.split()[TYPE] == "CDS":
                line = annot_f.readline()

    return gene_model_dict

def get_gene_model_sequences():
    pass

if __name__ == "__main__":
    pass












    # the gene name is encoded for here, the species is externally supplied, and we
    # can actually enforce range and even genomic accession



    pass

if __name__ == "__main__":




    pass