import os
from Bio import Entrez

from Basic_Tools.lists_and_files import list_to_string
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from NCBI_Genome_Blaster.create_basic_XML_processor import get_directory

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
        self.relative_cds_ranges = []

    def add_to_range(self, range1, range2):
        self.cds_ranges.append((int(range1), int(range2)))
        if self.strand == "1":
            self.relative_cds_ranges.append((int(range1) - self.gene_range[0],
                                             int(range2) - self.gene_range[0]))
        else:
            self.relative_cds_ranges.append((self.gene_range[1] - int(range2),
                                             self.gene_range[1] - int(range1)))


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


class GeMoMaProcessor:

    def __init__(self, annotation_path, taxon_and_name, save_dir):

        self.taxon_name = taxon_and_name["taxon"]
        self.species_name = taxon_and_name["name"]

        self.annotation_path = annotation_path
        self.save_dir = save_dir

        self.gene_model_dict = {}


    def process_gemoma_results(self):

        self.read_annotation_file()
        self.get_gene_model_sequences()

    def read_annotation_file(self):

        """
        The whole file needs to be read and processed first, since there are multiple

        Gene model dictionary is in the form {gene_name: GeneModel}
        """

        # redundant but makes more sense this way
        self.gene_model_dict = {}

        annot_f = open(self.annotation_path, "r")

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
            if gene_name not in self.gene_model_dict or self.gene_model_dict[gene_name].score < score:

                # create the new model
                new_model = GeneModel(gff_arr[GENOMIC_ACCESSION], score,
                                      gff_arr[STRAND],
                                      gff_arr[BOUND1], gff_arr[BOUND2])

                self.gene_model_dict[gene_name] = new_model

                line = annot_f.readline()
                while line.strip() != "" and line[0] != "#" and line.split()[TYPE] != "gene":
                    gff_arr = line.split()
                    new_model.add_to_range(gff_arr[BOUND1], gff_arr[BOUND2])

                    line = annot_f.readline()

            else:
                line = annot_f.readline()
                while line.strip() != "" and line[0] != "#" and line.split()[TYPE] == "CDS":
                    line = annot_f.readline()


    def create_transcript_file(self, gene_name):

        gene_folder = get_directory(self.save_dir, gene_name)
        taxa_folder = get_directory(gene_folder, self.taxon_name)
        species_folder = get_directory(taxa_folder, self.species_name)

        transcript_filepath = os.path.join(species_folder, "Reference_GeMoMa_result" +
                                           ".fas")
        return open(transcript_filepath, "a")

    def get_gene_model_sequences(self):


        for gene in self.gene_model_dict:

            # create (or simply locate) the existing transcript file for output
            transcript_file = self.create_transcript_file(gene)

            # create the gene model
            model = self.gene_model_dict[gene]

            # note that if strand = 2, auto correction happens, and all of the ranges
            # change
            # this is addressed in the relative_cds_ranges
            gene_region_arr = ncbi_get_gene_sequence(model.accession,
                                                     model.gene_range[0],
                                                     model.gene_range[1],
                                                     model.strand)

            completed_range = 0
            for cds in model.relative_cds_ranges:

                cds_range = abs(cds[0] - cds[1]) + 1

                fasta_heading = ">" + gene + " " + self.species_name + \
                                " reference_mrna:placeholder" + \
                                " genome:" + model.accession + " " + \
                                str(completed_range + 1) + "-" + str(completed_range + cds_range)

                transcript_file.write(fasta_heading + "\n")

                seq = ""
                for i in range(cds[0], cds[1] + 1):
                    seq += gene_region_arr[i]

                transcript_file.write(seq + "\n")

                completed_range += cds_range


if __name__ == "__main__":

    Entrez.email = "xiaohan.xie@amil.utoronto.ca"

    annotation_path = "/Users/tonyx/Documents/chang_lab/final_annotation.gff"
    taxon_and_name = {"taxon": "test", "name": "species"}
    save_dir = "/Users/tonyx/Documents/chang_lab/gemom_parse_test"

    os.mkdir("/Users/tonyx/Documents/chang_lab/gemom_parse_test")

    processor = GeMoMaProcessor(annotation_path, taxon_and_name, save_dir)
    processor.process_gemoma_results()




    # the gene name is encoded for here, the species is externally supplied, and we
    # can actually enforce range and even genomic accession




