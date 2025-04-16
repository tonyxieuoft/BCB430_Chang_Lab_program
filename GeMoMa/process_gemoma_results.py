"""
This script processes GeMoMa's results. GeMoMa's output comes in the form of genomic ranges, and we need to extract
the sequences that correspond to these genomic ranges for dN/dS analysis
"""

import os
from Bio import Entrez

from Basic_Tools.lists_and_files import list_to_string
from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from HSP_Selector_and_Processor.create_basic_XML_processor import get_directory

# define constants within GeMoMa's GFF output
GENOMIC_ACCESSION = 0
TYPE = 2
BOUND1 = 3
BOUND2 = 4
STRAND = 6
ANNOTATIONS = 8


class GeneModel:
    """
    Here, we define a class that represents a single GeMoMa gene model. We note that a single reference may produce
    multiple gene models, and we differentiate the best using their score (row BLAST score, I believe)
    """
    def __init__(self, accession, score, strand, gene_bound1, gene_bound2):
        """
        Instantiating the gene model requires several things, among which includes:

        1) accession: The genomic accession to which the Gene Model maps (necessary to extract the sequence)
        2) score: the raw alignment score associated with the Gene Model
        3) strand: the strand on which the gene model lies
        4) gene_bound1 and gene_bound2: the minimum boundaries containing the entire genomic gene sequence
        """

        # set raw score
        self.score = score

        # set genomic accession
        self.accession = accession

        # set genomic gene range
        self.gene_range = (int(gene_bound1), int(gene_bound2))

        # set strand
        if strand == "+":
            self.strand = "1"
        else:
            self.strand = "2"

        # set the absolute cds ranges (corresponding to the positions of the subject genomic sequence)
        self.cds_ranges = []

        # set the relative cds ranges (where 1 indicates the beginning of the gene region
        self.relative_cds_ranges = []

    def add_to_range(self, range1, range2):
        """
        Add a new CDS region to the gene model at hand. Auto convert to relative CDS range
        """
        self.cds_ranges.append((int(range1), int(range2)))
        if self.strand == "1":
            self.relative_cds_ranges.append((int(range1) - self.gene_range[0],
                                             int(range2) - self.gene_range[0]))
        else:
            self.relative_cds_ranges.append((self.gene_range[1] - int(range2),
                                             self.gene_range[1] - int(range1)))


def get_info_from_annotation(annotation):

    """
    GFF files contain a column with many miscellaneous annotations, separated by semicolons. Here,
    "annotation" refers to a list containing the column's annotations split by the semi colons
    """

    # we are trying to find the gene name and score assocate with the annotation
    gene_name = ""
    score = 0

    # for each of the annotation's features, check if the beginning starts with "ID": if so, then it's a gene name
    # also, check if the beginning starts with "score": if so, then it is a raw score value
    for entry in annotation:
        entry_type = entry.split("=")[0]
        if entry_type == "ID":
            # in the form ID=genename_R1,   need to get rid of the R1
            gene_name = entry.split("=")[1].split("_")[0]
        if entry_type == "score":
            score = int(entry.split("=")[1])

    return gene_name, score


class GeMoMaProcessor:

    """
    Here, we define a class that processes GeMoMa results first by creating GeneModel objects, then extracts the
    sequences corresponding to genomic ranges
    """
    def __init__(self, annotation_path, taxon_and_name, save_dir):
        """
        The processor requires three inputs for initiation:

        1) annotation_path; the genomic annotation path produced by GeMoMa
        2) taxon_and_name: the taxon and scientific name of the species that GeMoMa is making predictions for
        3) save_dir: the directory to save the processed sequence results
        """

        # set the taxon and scientific name of the species
        self.taxon_name = taxon_and_name["taxon"]
        self.species_name = taxon_and_name["name"]

        # set the annotation path
        self.annotation_path = annotation_path
        self.save_dir = save_dir

        self.gene_model_dict = {}

    def process_gemoma_results(self):
        """
        Processing GeMoMa's results comes in two stages: first, we read the genomic annotation file and organize them
        into GeneModels; second, we take the GeneModels and extract sequences from the genomes based on the information
        they contain.
        """
        self.read_annotation_file()
        self.get_gene_model_sequences()

    def read_annotation_file(self):

        """
        The whole file needs to be read and processed first, since there are multiple potentially relevant gene models.
        As a next step, maybe grab all of these gene models and see if incorporating all decreases GeMoMa's risk of
        paralogous sequence inclusion.

        Anyways, the Gene model dictionary produces key-value pairings in in the form {gene_name: GeneModel}
        """
        # redundant but makes more sense this way
        self.gene_model_dict = {}

        # open the GeMoMa genomic annotation file for reading
        annot_f = open(self.annotation_path, "r")

        # read the first few lines to skip over irrelevant parts
        line = annot_f.readline()
        line = annot_f.readline()
        line = annot_f.readline()

        while line.strip() != "":
            # now, we iterate through the file

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

        """
        Create a transcript file to store the results. This process is exactly the same as the Chang Lab algorithm's
        process.
        """

        # successfully create or get the path to gene, taxon and species folder (NEPR format)
        gene_folder = get_directory(self.save_dir, gene_name)
        taxa_folder = get_directory(gene_folder, self.taxon_name)
        species_folder = get_directory(taxa_folder, self.species_name)

        transcript_filepath = os.path.join(species_folder, "Reference_GeMoMa_result" +
                                           ".fas")
        return open(transcript_filepath, "a")

    def get_gene_model_sequences(self):
        """
        Given that we have created GeneModel objects after reading the genomic annotation file, we now extract the
        sequences that the annotations map to.
        """
        # for each of our gene models:
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

            # now, we have pulled out the entire gene region, and need to zone in specifically on the CDS regions
            completed_range = 0
            for cds in model.relative_cds_ranges:
                # get each relative cds range
                cds_range = abs(cds[0] - cds[1]) + 1

                # create a fasta heading in NEPR format
                fasta_heading = ">" + gene + " " + self.species_name + \
                                " reference_mrna:placeholder" + \
                                " genome:" + model.accession + " " + \
                                str(completed_range + 1) + "-" + str(completed_range + cds_range)

                transcript_file.write(fasta_heading + "\n")

                # get the sequence in the gene region corresponding to the coding sequence
                seq = ""
                for i in range(cds[0], cds[1] + 1):
                    # note that the above range implies cds[0], cds[0] + 1, ... , cds[1] - 1, cds[1]
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




