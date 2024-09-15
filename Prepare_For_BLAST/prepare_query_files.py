import os
from abc import abstractmethod
from typing import Dict, List

from Basic_Tools.taxonomy_browser import get_taxonomy_lineage
from Basic_Tools.lists_and_files import file_to_list, list_to_string
from Basic_Tools.basic_dictionaries import dict_get_values
from Prepare_For_BLAST.taxa_assigner import get_assignments
from Quality_Checking.get_longest_transcript import get_longest_transcript


class BlastPreparer:
    """

    Attributes
    ----------
    ref_seq_path: str
        path pointing to a directory of reference sequences
    taxa_to_ref_species: Dict[str,List]
    genes_to_available_species:

    """

    def __init__(self, ref_seq_path, save_path):

        self.ref_seq_path = ref_seq_path
        self.save_path = save_path

        self.taxa_to_ref_species = None
        self.genes_to_available_species = None

        self.lineage_dict = None
        self.taxid_codes = None

        self.complete_reference_species = None
        self.taxa_blast_order = None

    def get_taxa_blast_order(self):
        return self.taxa_blast_order

    def get_complete_ref_species(self):
        return self.complete_reference_species

    def set_reference_species(self) -> None:
        """
        Given a ref_seq_path pointing to a folder of reference sequences, identifies
        the reference species present and organizes them by assigned_taxa

        :param ref_seq_path: a string path pointing to a folder of reference
        sequences generated by the NCBI_Exon_Puller package.
        :return: a dictionary where keys are assigned_taxa and values are reference species
        belonging to a taxon.
        """
        # species encountered so far. hashmap data structure for quick verification
        # that a species has already been encountered.
        encountered_species = {}
        # output dictionary
        taxon_to_species = {}
        # iterate through folder structure
        # gene level
        for gene_folder in os.listdir(self.ref_seq_path):
            gene_path = os.path.join(self.ref_seq_path, gene_folder)
            # taxon level
            for taxon_folder in os.listdir(gene_path):
                taxon_path = os.path.join(gene_path, taxon_folder)
                # species level
                for species_folder in os.listdir(taxon_path):
                    # if we haven't already added the species to the output
                    if species_folder not in encountered_species:
                        encountered_species[species_folder] = True
                        # if the taxon of the species is already in the output,
                        # just append
                        if taxon_folder.upper() in taxon_to_species:
                            taxon_to_species[taxon_folder.upper()].append(
                                species_folder)
                        else:
                            # otherwise, create a new key
                            taxon_to_species[taxon_folder.upper()] = [species_folder]

        self.taxa_to_ref_species = taxon_to_species

    def set_lineage_dict_and_taxid_codes(self):
        """
        Given the reference species and associated assigned_taxa pulled out via the
        NCBI_Exon_Puller package, find lineages/taxids for each of them.

        :param taxa_to_species_dict: a dictionary, where keys are assigned_taxa and values
        are reference species within that assigned_taxa
        :return:
        """
        # get a list of all assigned_taxa present among reference species

        taxa_list = []
        for taxon in self.taxa_to_ref_species:
            taxa_list.append(taxon.upper())
        # get a list of all reference_species present
        species_list = dict_get_values(self.taxa_to_ref_species)

        # a list of assigned_taxa/species to lookup lineages for
        # ideally, only call get_taxonomy_lineage once, because it takes time to
        # contact the website
        species_string = list_to_string(taxa_list + species_list, "\n")
        print("contacting the NCBI taxonomy API...")
        lineage_dict = get_taxonomy_lineage(species_string)

        # a dictionary where keys are assigned_taxa and values are taxids corresponding to
        # the assigned_taxa

        taxid_codes = {}
        lineage_keys = list(lineage_dict.keys())
        # remove assigned_taxa from the results, we want lineage_dict to just have species
        for key in lineage_keys:
            if key.upper() in taxa_list:  # slightly inefficient, O(T), but T is pretty short
                taxid_codes[key.upper()] = lineage_dict[key][0]
                # only pop if it's not a species...
                if len(key.split()) == 1:
                    lineage_dict.pop(key, None)

        self.lineage_dict = lineage_dict
        self.taxid_codes = taxid_codes

    def prepare_query_files(self, auto_assign):
        """
        Prepare query files based on reference sequences pulled using the
        NCBI_Exon_Puller module.

        :param auto_assign: '1' for automatic delegation of sub-branches of an
        overhead taxon to reference species in the taxon. '0' for manual assignment.
        :param auto_fill_in: '1' for automatic pulling from the closest available
        species when a species is missing sequences for a gene. '0' for automatic
        filling.
        :param ref_seq_path: The directory containing reference sequences from the
        NCBI_Exon_Puller in the original folder hierarchy.
        :param save_path: A directory to save query files to.
        """
        # get a dictionary where the keys are assigned_taxa and values are
        # reference species
        # corresponding to the assigned_taxa
        self.set_reference_species()
        # get the lineage dictionary and taxids for the reference species
        self.set_lineage_dict_and_taxid_codes()

        # get assignments for the reference species. this is based on the
        # lineage
        # dictionary if automatic
        assignments = get_assignments(auto_assign, self.lineage_dict,
                                      self.taxid_codes,
                                      self.taxa_to_ref_species)

        # for each reference species -> sub-taxon assignment, make a query file
        # for it to use in BLAST. Sometimes, the reference species are missing
        # sequences for specific genes, and automatic/manual pulling from other
        # available species is offered.
        self.taxa_blast_order = []
        self.complete_reference_species = []

        print(assignments)

        for assignment in assignments:
            if self.single_query_assembly(assignment[0], assignment[1]):
                self.complete_reference_species.append(assignment[0])
            self.taxa_blast_order.append(assignment[1])

    def select_fill_in(self, species_name,
                       available_species: List[Dict]) -> Dict:
        """
        Called when building a query sequence file for ref_species_name, but is missing
        sequences for some genes. Automatically pulls from available species the
        species with the closest lineage (as stored in NCBI).

        :param species_name: name of species to build a query sequence file for
        :param available_species: a list of dictionaries, each storing a species
        that has the missing gene and the assigned_taxa the species is a part of, for species
        that have the missing gene.
        :param lineage_dict: a dictionary where the keys are species and the values
        are lists corresponding to the species' lineages.
        :return: a dictionary storing the species with the closest lineage to
        ref_species_name
        """
        # Stores the lineage of ref_species_name in hashmap format for quick retrieval.
        # Each key represents a taxon in the lineage, and its value is the
        # corresponding depth (0 for species level)
        host_lineage_hash = {}

        host_depth = 0
        for taxon in self.lineage_dict[species_name]:
            host_lineage_hash[taxon] = host_depth
            host_depth += 1

        # trace up the lineages of available species to determine how similar they
        # are to ref_species_name. store the closest one to return later
        closest_depth = float('inf')
        closest_species = ""
        for organism in available_species:

            # trace up the species lineage, stopping when encountering the first
            # assigned_taxa shared with the host
            avail_depth = 0
            curr_taxon = self.lineage_dict[organism["species"]][avail_depth]
            while curr_taxon not in host_lineage_hash:
                avail_depth += 1
                curr_taxon = self.lineage_dict[organism["species"]][avail_depth]

            # if the depth is smaller than the closest depth so far, this species
            # is more similar
            if host_lineage_hash[curr_taxon] < closest_depth:
                closest_depth = host_lineage_hash[curr_taxon]
                closest_species = organism

        return closest_species

    @abstractmethod
    def single_query_assembly(self, ref_species_name: str,
                              assigned_taxa: str) -> bool:
        pass


class ExonBlastPreparer(BlastPreparer):

    def single_query_assembly(self, ref_species_name: str,
                              assigned_taxa: str) -> bool:
        """
        Assumes that exons for reference sequences have been pulled out and are in
        the folder format: general folder -> gene -> taxon -> species -> transcript.
        Concatenate all results for one species into a single query_file, and name the
        query_file after the taxon that the query will be blasted against. If no reference
        sequence of a species exists for a given gene, asks the user to input an
        alternative species to draw the sequence from (or finds one automatically).

        :param autofill: '1' for automatically pulling from the closest available
        species when a sequence is missing, '0' for manual user input.
        :param ref_species_name: the species to concatenate results for
        :param reference_seq_path: the general folder containing reference sequences
        :param save_path: str
        :param assigned_taxa: the subject_taxa that the query will be blasted against
        :return: True iff the species has reference sequences for every
        specified gene. Also, a query_file will be created in the specified save_path
        """

        # query query_file to contain reference sequences to blast against a taxon
        query_file = open(os.path.join(self.save_path, assigned_taxa + ".fas"),
                          "a")

        complete_species = True

        # follows through the folder structure
        for gene in os.listdir(self.ref_seq_path):
            gene_path = os.path.join(self.ref_seq_path, gene)

            species_found = False
            available_species = []

            for taxa_folder in os.listdir(gene_path):
                taxa_path = os.path.join(gene_path, taxa_folder)

                for species_folder in os.listdir(taxa_path):
                    species_path = os.path.join(taxa_path, species_folder)

                    available_species.append({"taxon": taxa_folder,
                                              "species": species_folder})
                    # if the species is the one we are looking for
                    if species_folder.upper() == ref_species_name.upper():

                        self.write_to_query(species_path, query_file)
                        species_found = True

            if not species_found:

                complete_species = False
                to_select = self.select_fill_in(ref_species_name,
                                                available_species)

                species_path = os.path.join(gene_path, to_select["taxon"],
                                            to_select["species"])
                self.write_to_query(species_path, query_file)

        return complete_species

    def write_to_query(self, species_path, query_file):

        transcript_file = get_longest_transcript(species_path)['file_name']
        transcript_path = os.path.join(species_path, transcript_file)
        query_file.write(open(transcript_path, "r").read() + "\n")


class FullBlastPreparer(BlastPreparer):

    def __init__(self, ref_seq_path, save_path):

        super().__init__(ref_seq_path, save_path)
        self.queries_to_genes_to_exons = {}

    def get_queries_to_genes_to_exons(self):
        return self.queries_to_genes_to_exons

    def single_query_assembly(self, ref_species_name: str,
                              assigned_taxa: str) -> bool:
        """
        Assumes that exons for reference sequences have been pulled out and are in
        the folder format: general folder -> gene -> taxon -> species -> transcript.
        Concatenate all results for one species into a single query_file, and name the
        query_file after the taxon that the query will be blasted against. If no reference
        sequence of a species exists for a given gene, asks the user to input an
        alternative species to draw the sequence from (or finds one automatically).

        :param autofill: '1' for automatically pulling from the closest available
        species when a sequence is missing, '0' for manual user input.
        :param ref_species_name: the species to concatenate results for
        :param reference_seq_path: the general folder containing reference sequences
        :param save_path: str
        :param assigned_taxa: the subject_taxa that the query will be blasted against
        :return: True iff the species has reference sequences for every
        specified gene. Also, a query_file will be created in the specified save_path
        """

        # query query_file to contain reference sequences to blast against a taxon
        query_file = open(os.path.join(self.save_path, assigned_taxa + ".fas"),
                          "a")
        self.queries_to_genes_to_exons[assigned_taxa] = {}
        genes_to_exons = self.queries_to_genes_to_exons[assigned_taxa]

        complete_species = True

        # follows through the folder structure
        for gene in os.listdir(self.ref_seq_path):
            gene_path = os.path.join(self.ref_seq_path, gene)

            species_found = False
            available_species = []

            for taxa_folder in os.listdir(gene_path):
                taxa_path = os.path.join(gene_path, taxa_folder)

                for species_folder in os.listdir(taxa_path):
                    species_path = os.path.join(taxa_path, species_folder)

                    available_species.append({"taxon": taxa_folder,
                                              "species": species_folder})
                    # if the species is the one we are looking for
                    if species_folder.upper() == ref_species_name.upper():

                        genes_to_exons[gene] = []
                        self.write_to_query(species_path, query_file, genes_to_exons[gene])
                        species_found = True

            if not species_found:

                complete_species = False
                to_select = self.select_fill_in(ref_species_name,
                                                available_species)

                species_path = os.path.join(gene_path, to_select["taxon"],
                                            to_select["species"])
                genes_to_exons[gene] = []
                self.write_to_query(species_path, query_file, genes_to_exons[gene])

        return complete_species

    def write_to_query(self, species_path, query_file, exon_list):

        transcript_file = get_longest_transcript(species_path)['file_name']
        transcript_path = os.path.join(species_path, transcript_file)

        t_length = transcript_file.split("_")[0]

        t_file = open(transcript_path, "r")
        t_line = t_file.readline()

        exon_section = 0
        sample_header = t_line.split(" ")
        while len(sample_header[exon_section]) < 6 or \
                sample_header[exon_section][:6] != "genome":
            exon_section += 1
        exon_section += 1

        new_fasta_heading = list_to_string(sample_header[:exon_section], " ") + " 1-" + t_length + "\n"

        # just think about exon files for now
        # this wouldn't work if there are gaps
        full_sequence = ""
        while t_line != "":
            if t_line[0] == ">":
                exons = t_line.split(" ")[exon_section].split("-")
                exon_list.append([int(exons[0]), int(exons[1])])

            elif t_line.strip():
                full_sequence += t_line.strip()
            t_line = t_file.readline()

        full_sequence += "\n"
        query_file.write(new_fasta_heading + full_sequence + "\n")


if __name__ == "__main__":
    pass
    # path = r"C:\Users\tonyx\Downloads\api_pull_complete_results8"
    # save_path = r"C:\Users\tonyx\Downloads\test_taxa6"

    # if not os.path.isdir(save_path):
    #    os.mkdir(save_path)
    # print(prepare_query_files(1, 1, path, save_path))
    # print(get_reference_species(r'C:\Users\tonyx\Downloads\NCBI_exon_pull_results (2)'))
