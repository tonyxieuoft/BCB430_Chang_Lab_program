import os
from typing import Dict, List

from Basic_Tools.numeric_user_input import numeric_user_input
from Prepare_For_BLAST.get_longest_transcript import get_longest_transcript


def select_fill_in_manual(species_name: str, gene_name: str, gene_path: str,
                          subject_taxa: str) -> Dict or None:
    """
    Called when building a query sequence file for taxa, and the reference
    species being collected does not have a particular gene.

    :param gene_name: name of the missing gene
    :param gene_path: the filepath to the folder containing sequences for that
    gene
    :param species_name: name of the species missing the gene
    :param subject_taxa: the taxon that the query sequences currently being built will
    be BLASTed against
    :return: a dictionary storing the taxon and species of the user-selected
    species, or None if no species is selected
    """

    print("Gene: " + gene_name + " not found for " + species_name +
          ", which will be used to BLAST against " + subject_taxa)

    # given that the user has requested to fill in missing sequences manually,
    # asks the user whether the would like to fill in something for this
    # particular instance
    choose_alt_statement = "Would you like to use an alternative species for " + gene_name + "? Press 1 for yes, and 2 for no."
    choose_alt = numeric_user_input(1, 2, choose_alt_statement)

    if choose_alt == 1:

        taxa_present = os.listdir(gene_path)

        # displays subject_taxa and asks for user input
        print("Available subject_taxa:")
        for a in range(len(taxa_present)):
            print("(" + str(a) + ")" + " " + taxa_present[a])
        taxa_choice = numeric_user_input(0, len(taxa_present), "Enter a number to choose")

        # path to check out species
        taxa_path = os.path.join(gene_path, taxa_present[taxa_choice])

        # lists species that have a reference sequence for the gene
        species_present = os.listdir(taxa_path)

        # displays species and asks for user input
        print("Available reference species for this taxon:")
        for a in range(len(species_present)):
            print("(" + str(a) + ")" + " " + species_present[a])
        print("(" + str(len(species_present)) + ") None of the above, and move on")
        species_choice = numeric_user_input(0, len(species_present), "Enter a number to choose")

        # grabs the exons of the gene from the specified species
        if species_choice != len(species_present):
            print("obtaining from " + species_present[species_choice] + "...")
            return {"taxon": taxa_present[taxa_choice],
                    "species": species_present[species_choice]}

    else:
        print("moving on...")

    return None


def select_fill_in_auto(species_name, available_species: List[Dict],
                        lineage_dict: Dict) -> Dict:
    """
    Called when building a query sequence file for species_name, but is missing
    sequences for some genes. Automatically pulls from available species the
    species with the closest lineage (as stored in NCBI).

    :param species_name: name of species to build a query sequence file for
    :param available_species: a list of dictionaries, each storing a species
    that has the missing gene and the taxa the species is a part of, for species
    that have the missing gene.
    :param lineage_dict: a dictionary where the keys are species and the values
    are lists corresponding to the species' lineages.
    :return: a dictionary storing the species with the closest lineage to
    species_name
    """
    # Stores the lineage of species_name in hashmap format for quick retrieval.
    # Each key represents a taxon in the lineage, and its value is the
    # corresponding depth (0 for species level)
    host_lineage_hash = {}

    host_depth = 0
    for taxon in lineage_dict[species_name]:
        host_lineage_hash[taxon] = host_depth
        host_depth += 1

    # trace up the lineages of available species to determine how similar they
    # are to species_name. store the closest one to return later
    closest_depth = float('inf')
    closest_species = ""
    for organism in available_species:

        # trace up the species lineage, stopping when encountering the first
        # taxa shared with the host
        avail_depth = 0
        curr_taxon = lineage_dict[organism["species"]][avail_depth]
        while curr_taxon not in host_lineage_hash:
            avail_depth += 1
            curr_taxon = lineage_dict[organism["species"]][avail_depth]

        # if the depth is smaller than the closest depth so far, this species
        # is more simialr
        if host_lineage_hash[curr_taxon] < closest_depth:
            closest_depth = host_lineage_hash[curr_taxon]
            closest_species = organism

    # print("Species missing: " + species_name + " closest alternative: " + closest_species["species"])

    return closest_species


def single_query_assembly(autofill, species_name: str, taxa: str,
                          reference_seq_path: str, save_path: str,
                          lineage_dict) -> bool:
    """
    Assumes that exons for reference sequences have been pulled out and are in
    the folder format: general folder -> gene -> taxon -> species -> transcript.
    Concatenate all results for one species into a single file, and name the
    file after the taxon that the query will be blasted against. If no reference
    sequence of a species exists for a given gene, asks the user to input an
    alternative species to draw the sequence from (or finds one automatically).

    :param autofill: '1' for automatically pulling from the closest available
    species when a sequence is missing, '0' for manual user input.
    :param species_name: the species to concatenate results for
    :param reference_seq_path: the general folder containing reference sequences
    :param save_path: str
    :param taxa: the subject_taxa that the query will be blasted against
    :return: True iff the species has reference sequences for every
    specified gene. Also, a file will be created in the specified save_path
    """

    # query file to contain reference sequences to blast against a taxon
    file = open(os.path.join(save_path, taxa + ".fas"), "a")

    complete_species = True

    # follows through the folder structure
    for gene_folder in os.listdir(reference_seq_path):
        gene_path = os.path.join(reference_seq_path, gene_folder)

        # just in case no reference sequence exists for the species for a given
        # gene
        species_found = False
        available_species = []

        for taxa_folder in os.listdir(gene_path):
            taxa_path = os.path.join(gene_path, taxa_folder)

            for species_folder in os.listdir(taxa_path):
                species_path = os.path.join(taxa_path, species_folder)

                available_species.append({"taxon": taxa_folder,
                                          "species": species_folder})
                # if the species is the one we are looking for
                if species_folder.upper() == species_name.upper():

                    if len(os.listdir(species_path)) != 0:
                        # get the "best" transcript for the species, then append it
                        # to the query file we are building up
                        transcript_file = get_longest_transcript(species_path)
                        transcript_path = os.path.join(species_path, transcript_file)
                        file.write(open(transcript_path, "r").read() + "\n")
                        species_found = True

        if not species_found:
            complete_species = False
            # determine the alternative species to pull from
            if autofill == 1:
                # automatic
                to_select = select_fill_in_auto(species_name, available_species, lineage_dict)
            else:
                # manual
                to_select = select_fill_in_manual(species_name, gene_folder, gene_path, taxa)

            species_path = os.path.join(gene_path, to_select["taxon"], to_select["species"])
            transcript_file = get_longest_transcript(species_path)

            if transcript_file != "":
                transcript_path = os.path.join(species_path, transcript_file)
                file.write(open(transcript_path, "r").read() + "\n")
            else:
                print("Something went wrong")

    return complete_species

