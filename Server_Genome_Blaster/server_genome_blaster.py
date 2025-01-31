import os
from abc import abstractmethod

from Basic_Tools.lists_and_files import list_to_string
from Basic_Tools.taxonomy_browser import get_taxa_taxids
from NCBI_Genome_Blaster.assemble_blast_result_sequences import BlastXMLParser, \
    ExonBlastXMLParser, FullBlastXMLParser
from NCBI_Genome_Blaster.improved_assembler import ImprovedExonParser
from Server_Genome_Blaster.genome_downloader import SPECIES_DATA_FILENAME, \
    ServerGenomeDownloader, \
    read_species_data


class ServerGenomeBlaster:

    # TODO make a Genome Blaster class?
    def __init__(self, save_path, queries_path, taxa_blast_order,
                 complete_reference_species, genome_storage_path, taxa_to_codes):

        self.save_path = save_path

        self.queries_path = queries_path

        self.blast_order_dict = {}
        self.set_blast_order_dict(taxa_blast_order)

        self.species_so_far = {}
        self.set_species_so_far(complete_reference_species)

        self.taxids_to_taxa = {}
        for taxon in taxa_to_codes:
            self.taxids_to_taxa[taxa_to_codes[taxon][0]] = taxon

        self.genome_storage_path = genome_storage_path

    # TODO not sure if this is needed
    def set_blast_order_dict(self, taxa_blast_order):

        for i in range(len(taxa_blast_order)):
            self.blast_order_dict[taxa_blast_order[i]] = i

    def set_species_so_far(self, complete_reference_species):

        for species in complete_reference_species:
            self.species_so_far[species] = True

    def download_new_genomes(self):

        print("downloading new genomes...")
        taxa_list = list(self.taxids_to_taxa.keys())

        downloader = ServerGenomeDownloader(self.save_path, taxa_list,
                                            self.genome_storage_path)
        downloader.set_accessions_to_download()
        if len(downloader.get_accessions_to_download()) == 0:
            print("No new accessions to download!")
        else:
            downloader.download_genomes()

    def blast_genomes(self, expect_value: str):

        available_genome_data = read_species_data(self.genome_storage_path)
        print("blast_order_dict: " + str(self.blast_order_dict))

        for genome in available_genome_data:

            curr_ref_taxon = None
            curr_prio = float("inf")
            overarching_taxon = None
            for taxon in genome["lineage"]:
                if taxon in self.taxids_to_taxa:
                    if overarching_taxon is None:
                        overarching_taxon = self.taxids_to_taxa[taxon]
                    else:
                        print(genome["name"] + " is a member of multiple inputted taxa, will be placed in " + overarching_taxon)

                if taxon in self.blast_order_dict and self.blast_order_dict[taxon] < curr_prio:
                    curr_ref_taxon = taxon
                    curr_prio = self.blast_order_dict[taxon]

            if curr_ref_taxon is not None:

                blast_db = os.path.join(self.genome_storage_path, "blast_db", list_to_string(genome["name"].split(), "_"))

                reference_filepath = os.path.join(self.queries_path, curr_ref_taxon + ".fas")
                xml_out_path = os.path.join(self.save_path, "temp.xml")
                print("blasting " + genome["name"] + " against " + curr_ref_taxon + "...")
                os.system("nice -5 blastn -db " + blast_db +
                          " -outfmt 5 -evalue " + str(expect_value) +
                          " -word_size 11 -gapopen 10 -gapextend 6 -reward 5 "
                          "-penalty -4 "
                          "-num_threads 10 "
                          #"-task dc-megablast "
                          "-soft_masking false "
                          "-dust no "
                          "-query " + reference_filepath +
                          " > " + xml_out_path)

                # parse the blast results
                taxon_and_name = {"taxon": overarching_taxon, "name": genome["name"]}
                print("parsing...")
                self.parse_blast_xml(xml_out_path, taxon_and_name)
                print("done parsing")
                os.remove(xml_out_path)

    @abstractmethod
    def parse_blast_xml(self, file_to_analyze, curr_species):
        pass

class ServerExonGenomeBlaster(ServerGenomeBlaster):

    def parse_blast_xml(self, file_to_analyze, curr_species):

        parser = ExonBlastXMLParser(file_to_analyze, self.save_path,
                                    curr_species, on_server=True)
        parser.parse_blast_xml()

class ServerImprovedExonGenomeBlaster(ServerGenomeBlaster):

    def parse_blast_xml(self, file_to_analyze, curr_species):

        parser = ImprovedExonParser(file_to_analyze, self.save_path,
                                    curr_species, on_server=True)
        parser.parse_blast_xml()


class ServerFullGenomeBlaster(ServerGenomeBlaster):

    def __init__(self, save_path, queries_path, taxa_blast_order,
                 complete_reference_species, genome_storage_path,
                 taxa_to_codes, queries_to_genes_to_exons):
        super().__init__(save_path, queries_path, taxa_blast_order,
                         complete_reference_species, genome_storage_path,
                         taxa_to_codes)
        self.queries_to_genes_to_exons = queries_to_genes_to_exons

    def parse_blast_xml(self, file_to_analyze, curr_species):

        parser = FullBlastXMLParser(file_to_analyze, self.save_path,
                                    curr_species,
                                    self.queries_to_genes_to_exons)
        parser.parse_blast_xml()





    # I can probably integrate this by early next week
