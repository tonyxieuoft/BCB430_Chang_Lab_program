"""
This file contains scripts automating the process of exon-by-exon reference BLAST against a subject genome. This is
the first step of the Chang Lab prediction program: using nucleotide-based, specificity-oriented BLAST to obtain a pool
of exons to construct putative gene models

"""

import os
from abc import abstractmethod

from Basic_Tools.lists_and_files import list_to_string
from Basic_Tools.taxonomy_browser import get_taxa_taxids
from GeMoMa.process_gemoma_results import GeMoMaProcessor
from HSP_Selector_and_Processor.gene_model_selector_and_processor import ImprovedExonParser
from Server_Genome_Blaster.genome_downloader import SPECIES_DATA_FILENAME, \
    ServerGenomeDownloader, \
    read_species_data
from Server_Genome_Blaster.genome_downloader_fasta import ServerFastaGenomeDownloader

class ServerGenomeBlaster:

    """
    Here, we define the main server genome blaster class. The exon-based approach merely extends the functionality
    exhibited here.
    """

    # TODO make a Genome Blaster class?
    def __init__(self, save_path, queries_path, taxa_blast_order,
                 complete_reference_species, genome_storage_path, taxa_to_codes):
        """
        Initializing the server genome blaster class requires the following as input:

        1) save_path: the path to which gene predictions made by "gene_model_selector_and_processor.py" will be
        outputted
        2) queries_path: the path containing reference sequences, organized for genome-wide species-based BLAST
        3) taxa_blast_order: the order in which taxa of interest will be blasted. Consult the ReadMe for details on
        why this matters
        4) complete_reference_species: the list of total species for which the reference sequences come from
        5) genome_storage_path: a path to a directory all genomes in blast-databases format (made through the
        makeblastdb command)
        6) taxa_to_codes: a dictionary mapping the names of taxa to their corresponding NCBI lineage codes
        """

        # output path definition
        self.save_path = save_path

        # reference sequence path definition
        self.queries_path = queries_path

        # the order in which taxa of interest will be blasted (again, consult the ReadMe for more details)
        self.blast_order_dict = {}
        self.set_blast_order_dict(taxa_blast_order)

        # define the list of total species for which the reference sequences come from
        self.species_so_far = {}
        self.set_species_so_far(complete_reference_species)

        # here, we map NCBI lineage codes to their corresponding taxonomic groups
        self.taxids_to_taxa = {}
        for taxon in taxa_to_codes:
            self.taxids_to_taxa[taxa_to_codes[taxon][0]] = taxon

        self.genome_storage_path = genome_storage_path

    def set_blast_order_dict(self, taxa_blast_order):
        """
        Define the order in which taxa will be blasted against (again, see the ReadMe for more details)
        """
        for i in range(len(taxa_blast_order)):
            self.blast_order_dict[taxa_blast_order[i]] = i

    def set_species_so_far(self, complete_reference_species):
        """
        Create a dictionary containing all of the species encountered so far during the automated server BLAST
        """
        for species in complete_reference_species:
            self.species_so_far[species] = True

    def download_new_genomes(self):
        """
        As the function title suggests, download genomes for all taxa of interest that currently do not have genomes
        on the server. This relies on creating GenomeDownloader classes, which are referred to in depth in the
        "genome_downloader.py" file
        """

        # check to see if new genomes need to be downloaded
        print("downloading new genomes...")

        taxa_list = list(self.taxids_to_taxa.keys())

        # instantiate a Genome Downloader class to check if new genomes need to be downloaded, and if so, to automate
        # the process of downloading them. See the "genome_downloader python file for more details
        downloader = ServerGenomeDownloader(self.save_path, taxa_list,
                                            self.genome_storage_path)
        downloader.set_accessions_to_download()
        if len(downloader.get_accessions_to_download()) == 0:
            print("No new accessions to download!")
        else:
            # automated process of downloading new genomes.
            downloader.download_genomes()

    def blast_genomes(self, expect_value: str):

        """
        Automate the process of nucleotide BLAST against genomes from taxa of interest.
        """

        # get all genomes downloaded onto the server (after the GenomeDownloader has done its job
        available_genome_data = read_species_data(self.genome_storage_path)
        # theoretically, available_genome_data should contain data about all of NCBI's reference genomes for the taxa
        # of interest to BLAST
        print("blast_order_dict: " + str(self.blast_order_dict))

        # cycle through the genomes
        for genome in available_genome_data:

            # here, we determine which reference species is best suited to BLAST against the particular subject genome

            # set the current reference taxon
            curr_ref_taxon = None

            # set the priority of the current reference taxon (that the higher the number, the lower the priority)
            curr_prio = float("inf")

            # the "overarching taxon" variable is used to categorize what directory to place the gene predictions for
            # the subject genome in
            overarching_taxon = None

            # walk through the subject genome's species' entire lineage, until a code corresponding to a reference
            # species' lineage is reached
            for taxon in genome["lineage"]:
                if taxon in self.taxids_to_taxa:

                    # if the overarching taxon isn't set, set it
                    if overarching_taxon is None:
                        overarching_taxon = self.taxids_to_taxa[taxon]
                    else:
                        # if we've encountered multiple overarching taxa, it means the user-specified taxa of interest
                        # containined species overlapping in different taxa
                        print(genome["name"] + " is a member of multiple inputted taxa, will be placed in " + overarching_taxon)

                # reset the reference species that will blast our subject of interest, if the reference species is of
                # other preset priority
                if taxon in self.blast_order_dict and self.blast_order_dict[taxon] < curr_prio:
                    curr_ref_taxon = taxon
                    curr_prio = self.blast_order_dict[taxon]

            # if we have identified a good reference species for the taxon:
            if curr_ref_taxon is not None:

                # identify a path to the blast database
                blast_db = os.path.join(self.genome_storage_path, "blast_db", list_to_string(genome["name"].split(), "_"))

                # identify the path to the file containing the reference sequences
                reference_filepath = os.path.join(self.queries_path, curr_ref_taxon + ".fas")

                # the path of the output XML file (this will be fed into the gene model creater and HSP processor!!)
                xml_out_path = os.path.join(self.save_path, "temp.xml")

                # print statement to let the user know BLAST is about to occur
                print("blasting " + genome["name"] + " against " + curr_ref_taxon + "...")

                # here, we run BLAST with word size 11, gap open 4, gap extend 1, match 1, mismatch 1: this is optimized
                # towards sequences with roughly 75% sequence similarity (as is roughly the case for Elasmobranchii)
                os.system("nice -5 blastn -db " + blast_db +
                          " -outfmt 5 -evalue " + str(expect_value) +
                          " -word_size 11 -gapopen 4 -gapextend 1 -reward 1 "
                          "-penalty -1 "
                          "-num_threads 16 "
                          "-query " + reference_filepath +
                          " > " + xml_out_path)

                # grab the taxon and species name of the subject genome we just made predictions for
                taxon_and_name = {"taxon": overarching_taxon, "name": genome["name"]}
                print("parsing...")

                # parse BLAST XML runs the gene model creater and processor in "gene_model_selector_and_processor.py"
                self.parse_blast_xml(xml_out_path, taxon_and_name)
                print("done parsing")
                os.remove(xml_out_path)

    @abstractmethod
    def parse_blast_xml(self, file_to_analyze, curr_species):
        """
        A skeleton for parsing and processing output BLAST XML files. Expanded on in the
        "gene_model_selector_and_processor" file
        """
        pass

class GemomaRunner(ServerGenomeBlaster):
    """
    Create a separate class specifically for running GeMoMa
    """
    def __init__(self, save_path, queries_path, taxa_blast_order,
                 complete_reference_species, genome_storage_path, taxa_to_codes,
                 gff_path, gemoma_out):
        super().__init__(save_path, queries_path, taxa_blast_order,
                         complete_reference_species, genome_storage_path, taxa_to_codes)

        # here, aside from just the basic constructor, we need to also have a gff file (normally for genome
        # annotations, but for us it's transcript-level intron-exon boundaries.
        self.gff_path = gff_path
        self.gemoma_out = gemoma_out

    def download_new_genomes(self):

        """
        Instead of downloading BLAST databases, the GeMoMa runner downloads compelte fasta files. This is because
        GeMoMa insists on creating blast databases from scratch every time it is run anew.
        """

        print("downloading new genomes...")
        taxa_list = list(self.taxids_to_taxa.keys())

        # Here, yet again, relies on a version of the server genome downloader class to download genomes
        downloader = ServerFastaGenomeDownloader(self.save_path, taxa_list,
                                                 self.genome_storage_path)
        downloader.set_accessions_to_download()
        if len(downloader.get_accessions_to_download()) == 0:
            print("No new accessions to download!")
        else:
            downloader.download_genomes()

    def blast_genomes(self, expect_value: str):

        """
        A little bit of a misnaming...

        This function actually automates an instance of GeMoMa. Here, we run GeMoMa in default settings,
        feeding it an input GFF file as well as a curated fasta file.

        """

        # get available genome data on the server
        available_genome_data = read_species_data(self.genome_storage_path)
        print("blast_order_dict: " + str(self.blast_order_dict))

        # iterating through each genome, we see...
        for genome in available_genome_data:

            # observing each genome separately
            curr_ref_taxon = None
            curr_prio = float("inf")
            overarching_taxon = None

            # go through the lineage
            for taxon in genome["lineage"]:

                # if we have a match to one of our input taxa to predict sequences for
                if taxon in self.taxids_to_taxa:
                    # the first match is the overarching taxon
                    if overarching_taxon is None:
                        overarching_taxon = self.taxids_to_taxa[taxon]
                    else:
                        print(genome["name"] + " is a member of multiple inputted taxa, will be placed in " + overarching_taxon)

                # what to blast first (and everything else gets excluded, yes)
                # find the best prio (which I feel like should just be the first regardless)
                if taxon in self.blast_order_dict and self.blast_order_dict[taxon] < curr_prio:
                    curr_ref_taxon = taxon
                    curr_prio = self.blast_order_dict[taxon]

            # if we have identified a good reference sequence match for our subject genome:
            if curr_ref_taxon is not None:

                # get the path to the genome (in fasta format)
                genome_fasta = os.path.join(self.genome_storage_path, "genome_fastas",
                                            list_to_string(genome["name"].split(), "_") + ".fna")

                # get the path to where the reference sequences are being stored
                sequence_filepath = os.path.join(self.queries_path, curr_ref_taxon + ".fas")

                # get the path to the GFF path containing the exon-boundary information necessary for GeMoMa to function
                gff_filepath = os.path.join(self.gff_path, curr_ref_taxon + ".gff")

                # a print statement notifying the user what stage the program is on
                print("running " + genome["name"] + " against " + curr_ref_taxon + "...")

                # here, we run the command line version of GeMoMa as follows:
                os.system("GeMoMa -Xms5G -Xmx50G "
                          "GeMoMaPipeline t=" + genome_fasta +
                          " a=" + gff_filepath + " g="  + sequence_filepath +
                          " AnnotationFinalizer.r=NO threads=16 outdir=" + self.gemoma_out
                          )

                # GeMoMa returns three files, notably:

                # 1) genomic locations for the predictions it makes. This is honestly
                # a little bit of a pain, since we are only interesting in nucleotide sequences and not locations
                gemoma_annotations = os.path.join(self.gemoma_out, "final_annotation.gff")

                # 2) GeMoMa's predictions in protein format. Given the redundancy of the codon table, we can't convert
                # this back to nulcoetide format, so this is useless to us for dN/dS analyses
                gemoma_protein_predictions = os.path.join(self.gemoma_out, "predicted_proteins.fasta")

                # 3) GeMoMa's protocol: irrelevant for our purposes, but the path is necessary to note down so that we
                # can delete it.
                gemoma_protocol = os.path.join(self.gemoma_out, "protocol_GeMoMaPipeline.txt")

                # the name of the subject genome that we are running GeMoMa against
                taxon_and_name = {"taxon": overarching_taxon, "name": genome["name"]}
                print("parsing...")

                # here, we process GeMoMa's results (i.e., consult the NCBI API to extract the genomic ranges of
                # interest into nucleotide sequence format
                processor = GeMoMaProcessor(gemoma_annotations, taxon_and_name, self.save_path)
                processor.process_gemoma_results()

                #self.parse_blast_xml(xml_out_path, taxon_and_name)
                print("done parsing")

                os.remove(gemoma_annotations)
                os.remove(gemoma_protein_predictions)
                os.remove(gemoma_protocol)

class ServerImprovedExonGenomeBlaster(ServerGenomeBlaster):

    def parse_blast_xml(self, file_to_analyze, curr_species):

        parser = ImprovedExonParser(file_to_analyze, self.save_path,
                                    curr_species, on_server=True)
        parser.parse_blast_xml()


