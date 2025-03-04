import os
import sys

from Basic_Tools.lists_and_files import file_to_list

if __name__ == "__main__":

    NEPR_path = "/Users/tonyx/Documents/chang_lab/good_batoid_references"
    gff_file = "/Users/tonyx/Documents/chang_lab/gffs/combined2.gff"
    fasta_file = "/Users/tonyx/Documents/chang_lab/gffs/combined2.fas"

    gff_f = open(gff_file, "w")
    gff_f.write("##gff-version 3\n")

    fasta_f = open(fasta_file, "w")

    for gene in os.listdir(NEPR_path):

        gene_path = os.path.join(NEPR_path, gene)
        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            for species in os.listdir(taxon_path):
                species_path = os.path.join(taxon_path, species)

                for exon_file in os.listdir(species_path):

                    exon_file = os.path.join(species_path, exon_file)

                    exon_file_array = file_to_list(exon_file)
                    if len(exon_file_array) < 2:
                        print("not enough lines in the exon file")
                        sys.exit()

                    header = exon_file_array[0]
                    if header[0] != ">":
                        print("not in fasta heading")
                        sys.exit()

                    mrna_pos = -1

                    for i in range(len(header.split())):
                        if len(header.split()[i]) > 4 and header.split()[i][:4] == "mRNA":
                            mrna_pos = i
                            break

                    if mrna_pos == -1:
                        print("no mrna name")
                        sys.exit()

                    cds_pos = mrna_pos + 2

                    mrna_name = header.split()[mrna_pos][5:]
                    gene_name = header.split()[0][1:]

                    cds_ranges = []

                    seq = ""
                    for i in range(len(exon_file_array)):
                        if exon_file_array[i][0] == ">":
                            cds = exon_file_array[i].split()[cds_pos].split("-")
                            cds_ranges.append(cds)
                        else:
                            seq += exon_file_array[i].strip()

                    cds_len = cds_ranges[-1][1]
                    gff_f.write(gene_name + "\t" +
                                "." + "\t" +
                                "gene" + "\t" +
                                "1" + "\t" +
                                cds_len + "\t" +
                                "." + "\t" +
                                "+" + "\t" +
                                "." + "\t" +
                                "ID=gene-" + gene_name +"\n")

                    gff_f.write(gene_name + "\t" +
                                "." + "\t" +
                                "mRNA" + "\t" +
                                "1" + "\t" +
                                cds_len + "\t" +
                                "." + "\t" +
                                "+" + "\t" +
                                "." + "\t" +
                                "ID=rna-" + mrna_name + "\n")

                    for cds in cds_ranges:
                        gff_f.write(gene_name + "\t" +
                                    "." + "\t" +
                                    "CDS" + "\t" +
                                    cds[0] + "\t" +
                                    cds[1] + "\t" +
                                    "." + "\t" +
                                    "+" + "\t" +
                                    "." + "\t" +
                                    "Parent=rna-" + mrna_name + "\n")

                    fasta_f.write(">" + gene_name + " " + mrna_name + "\n")
                    fasta_f.write(seq + "\n")

    gff_f.write("###")




