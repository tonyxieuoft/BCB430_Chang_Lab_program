"""
This script converts information from reference sequences pulled from NCBI into a GFF format. This is necessary for
GeMoMa to run.
"""

import os.path
import sys

from Basic_Tools.lists_and_files import file_to_list

def exon_file_to_gff(exon_file: str, gff_file: str, fasta_file: str):
    """
    Converts a file containing exon references pulled from NCBI into GFF format.

    The function requires three inputs:

    1) exon_file: the file containing exons (note that the file must be in NEPR format, though)
    2) gff_file: a path to a growing GFF containing information about exon references (for input into GeMoMa)
    3) fasta_file: a path to a growing FASTA file containing full coding sequences (for input into GeMoMa)
    """
    # if we don't have a single exon in our exon file, quit the function
    exon_file_array = file_to_list(exon_file)
    if len(exon_file_array) < 2:
        print("not enough lines in the exon file")
        sys.exit()

    # if the first line of the exon file isn't a fasta heading:
    header = exon_file_array[0]
    if header[0] != ">":
        print("not in fasta heading")
        sys.exit()

    # we need to identify the position in the exon headers that contains the identifier from which the exon comes from
    mrna_pos = -1
    for i in range(len(header.split())):
        # iterate until the mRNA related header is found
        if len(header.split()[i]) > 4 and header.split()[i][:4] == "mRNA":
            mrna_pos = i
            break

    if mrna_pos == -1:
        print("no mrna name")
        sys.exit()

    # the position of the accession in the header tells the position where the full-coding sequence range of the query
    # is stored
    cds_pos = mrna_pos + 2

    # the first 4 characters of the mRNA position is "mRNA:", to get the identifier we take the characters at position
    # 5 and later
    mrna_name = header.split()[mrna_pos][5:]
    # the gene name is at position 0 in the header
    gene_name = header.split()[0][1:]

    # if the gff file does not exist yet, make a new one
    if not os.path.isfile(gff_file):
        gff_f = open(gff_file, "w")
        # all gff files must contain this at the top
        gff_f.write("##gff-version 3\n")
    else:
        # otherwise, add to an existing GFF file
        gff_f = open(gff_file, "a")

    # just get the cds_ranges and sequence from the file
    cds_ranges = []
    seq = ""

    # IMPORTANT: here, we acquire the coding sequence ranges that correspond to the exons -> these will directly be
    # used for GFF ranges
    for i in range(len(exon_file_array)):
        if exon_file_array[i][0] == ">":
            cds = exon_file_array[i].split()[cds_pos].split("-")
            cds_ranges.append(cds)
        else:
            seq += exon_file_array[i].strip()

    # the length of the full coding sequence is the ending range for the last exon
    cds_len = cds_ranges[-1][1]

    # write the general gene (this entity is necessary to give the CDS features something to point to)
    gff_f.write(gene_name + "\t" +
                "." + "\t" +
                "mRNA" + "\t" +
                "1" + "\t" +
                cds_len + "\t" +
                "." + "\t" +
                "+" + "\t" +
                "." + "\t" +
                "ID=" + gene_name + "\n")

    # write each of the CDS features (corresponding to the ranges of the exons relative to the transcript)
    # a little bit of a hack, since the genomic range isn't used
    for cds in cds_ranges:
        gff_f.write(gene_name + "\t" +
                    "." + "\t" +
                    "CDS" + "\t" +
                    cds[0] + "\t" +
                    cds[1] + "\t" +
                    "." + "\t" +
                    "+" + "\t" +
                    "." + "\t" +
                    "Parent=" + gene_name + "\n")

    fasta_f = open(fasta_file, "a")
    fasta_f.write(">" + gene_name + " " + mrna_name + "\n")
    fasta_f.write(seq + "\n")


if __name__ == "__main__":

    exon_file = "/Users/tonyx/Documents/chang_lab/converted_NEPR5/RARA/Carcharodon carcharias.fas"
    gff_file = "/Users/tonyx/Documents/chang_lab/gffs/RARA_test.gff"
    fasta_file = "/Users/tonyx/Documents/chang_lab/gffs/RARA_test.fas"

    exon_file_to_gff(exon_file,gff_file,fasta_file)



