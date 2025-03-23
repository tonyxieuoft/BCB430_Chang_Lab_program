import os.path
import sys

from Basic_Tools.lists_and_files import file_to_list

def exon_file_to_gff(exon_file: str, gff_file: str, fasta_file: str):
    # if we don't have a single fasta entry, leave (that's not going to be the case though, is it?
    exon_file_array = file_to_list(exon_file)
    if len(exon_file_array) < 2:
        print("not enough lines in the exon file")
        sys.exit()

    header = exon_file_array[0]
    if header[0] != ">":
        print("not in fasta heading")
        sys.exit()

    # TODO add back if mRNA is needed

    mrna_pos = -1

    for i in range(len(header.split())):
        if len(header.split()[i]) > 4 and header.split()[i][:4] == "mRNA":
            mrna_pos = i
            break

    if mrna_pos == -1:
        print("no mrna name")
        sys.exit()

    # the only reason I'm getting the mRNA position is bc it's relative tot he CDS
    cds_pos = mrna_pos + 2

    mrna_name = header.split()[mrna_pos][5:]
    gene_name = header.split()[0][1:]

    if not os.path.isfile(gff_file):
        gff_f = open(gff_file, "w")
        gff_f.write("##gff-version 3\n")
    else:
        gff_f = open(gff_file, "a")

    # just get the cds_ranges and sequence from the file
    cds_ranges = []
    seq = ""

    for i in range(len(exon_file_array)):
        if exon_file_array[i][0] == ">":
            cds = exon_file_array[i].split()[cds_pos].split("-")
            cds_ranges.append(cds)
        else:
            seq += exon_file_array[i].strip()

    # the length of the coding sequence is the last one
    cds_len = cds_ranges[-1][1]

    # technically in the mRNA category (seems like this is all we need though, LOL)
    gff_f.write(gene_name + "\t" +
                "." + "\t" +
                "mRNA" + "\t" +
                "1" + "\t" +
                cds_len + "\t" +
                "." + "\t" +
                "+" + "\t" +
                "." + "\t" +
                "ID=" + gene_name + "\n")

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



