"""
This file contains minimal code, and is only meant to instantiate a basic XML parser class for the gene model selector
contained in "gene_model_selector_and_processor" to run on
"""

import os
from abc import abstractmethod

from NCBI_Exon_Puller.ncbi_exon_puller import ncbi_get_gene_sequence
from Bio import Entrez
from Basic_Tools.xml_extraction import file_xml_to_dictionary
from Basic_Tools.xml_extraction import get_xml_list

SEQUENCE_INDICES_FROM_MRNA_TAG = 2

def get_directory(parent_directory: str, directory_name: str) -> str:
    """
    Returns either a new directory, or simply a full path specified by 'parent_directory / directory_name'

    :param parent_directory: The parent directory that will contain the new directory
    :param directory_name: The name of the directory to create
    :return:
    """

    # first, concatenate the parent directory and directory name to get the full path
    directory_path = os.path.join(parent_directory, directory_name)
    if not os.path.isdir(directory_path):
        os.mkdir(directory_path)

    return directory_path


class BlastXMLParser:

    """
    Here outlines a basic BLAST XML Parser class.
    """

    def __init__(self, file, save_dir, curr_species):
        """
        The class takes in three inputs:
        1) file: the path to the BLAST XML file
        2) save_dir: the directory to store the Chang Lab's gene predictions upon parsing the file
        3) curr_species: the species of the genome blasted against
        """

        # define the XML file path
        self.xml_filepath = file

        # define the output directory to save results
        self.save_dir = save_dir

        # get the taxonomic group of the species whose genome is getting blasted against
        self.taxon_name = curr_species['taxon']
        # get the species name that the genome corresponds to
        self.species_name = curr_species['name']

    def create_transcript_file(self, ref_transcript_var, gene_name):

        """
        Creates a transcript file to store predictions made by gene_model_selector_and_processor after it reads the
        BLAST input file. Note that this function only makes sense in the context in which it is called in the
        "gene_model_selector_and_processor" file

        @param ref_transcript_var the accession of the reference transcript for each a particular exon query belongs
        to (undefined here)

        @param gene_name the name of the gene for which the program is making predictions.
        """
        # either creates a new gene directory or gets the path an existing one
        gene_folder = get_directory(self.save_dir, gene_name)
        # same thing, but on the taxon level
        taxa_folder = get_directory(gene_folder, self.taxon_name)
        # same, but on the the species level
        species_folder = get_directory(taxa_folder, self.species_name)

        transcript_filepath = os.path.join(species_folder, "Reference_" +
                                           ref_transcript_var + ".fas")
        return open(transcript_filepath, "a")

    @abstractmethod
    def parse_blast_xml(self):
        """
        Here, a skeleton is defined for which the gene model to build upon.
        """
        pass


class ExonBlastXMLParser(BlastXMLParser):

    """
    Here, we define a basic Exon Parser class. The class is honestly an artifact from the multiple different
    iterations of the pipeline, only extending functionality by storing whether we are accessing a server version of
    BLAST or not
    """

    def __init__(self, file, save_dir, curr_species, on_server):
        # use the same constructor as the basic BlastXMLParser class
        super().__init__(file, save_dir, curr_species)
        self.on_server = on_server

    def parse_blast_xml(self):
        pass




if __name__ == "__main__":
    pass
    # xiaohan.xie@mail.utoronto.ca
    # C:\Users\tonyx\Downloads\main_test4
    # C:\Users\tonyx\Downloads\main_test3\genes.txt
    # C:\Users\tonyx\Downloads\main_test3\taxa.txt
    # C:\Users\tonyx\Downloads\main_test4\NCBI_exon_pull_results
    # parse_blast_xml(r"C:\Users\tonyx\Downloads\51JMZR2N016-Alignment.xml", r"C:\Users\tonyx\Downloads\xml_readtest", "altered2", "some weird one")





