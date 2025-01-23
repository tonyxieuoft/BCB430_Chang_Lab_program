import os
from typing import Dict, List

from Basic_Tools.lists_and_files import file_to_list


def auto_assign_taxa_to_ref(lineage_dct: Dict, ref_species: List[str],
                            overhead_taxon: str) -> List:
    """
    Given an overarching taxon (ex. cetacea, elasmobranchii), delegate
    sub-branches of the taxon for reference species to blast against. This is to
    ensure that query sequences are as similar to the subject as possible.

    :param lineage_dct: a dictionary where keys are species and values are lists
    corresponding to a species' lineage.
    :param ref_species: a list of reference species to delegate sub-assigned_taxa to
    :param overhead_taxon: the overarching taxon
    :return: a list of tuples each containing a reference species and the
    corresponding sub-taxon in the overhead_taxon it was assigned to. Order
    matters, in that species-taxon combos at the top of the list will be blasted
    first, and NOT be reblasted later.
    """
    lineage_assignments = {}
    result = []
    for ref in ref_species:
        ref_lineage = lineage_dct[ref]
        depth = 0
        # traces up the lineage and assigns them to itself until bumping into
        # a subject_taxa that has already been assigned
        while ref_lineage[depth] != overhead_taxon and \
                ref_lineage[depth] not in lineage_assignments:
            lineage_assignments[ref_lineage[depth]] = ref
            depth += 1

        if ref_lineage[depth] == overhead_taxon and \
                overhead_taxon not in lineage_assignments:
            # if the overhead taxon is reached for the first time
            lineage_assignments[overhead_taxon] = ref
            result.insert(0, [ref, overhead_taxon])
        else:
            # this happens after bumping into an already assigned subject_taxa
            # therefore, its subject_taxa is the one just below it
            if depth-1 >= 0:
                # if statement covers the edge case of a species duplicate in the taxon file
                # which would cause an index out of bounds error
                result.insert(0, [ref, ref_lineage[depth-1]])
            # insert at the front

    return result


def get_assignments(auto: int, lineage_dict: Dict[str, List[str]],
                    taxa_to_codes: Dict[str, str]) -> List:
    """
    Get the sub-branch assignments for each reference species.

    :param auto: '1' for automatic determination of assignments, '0' for
    manually entering assignments.
    :param lineage_dict: a dictionary where keys are species names and values
    are lists of assigned_taxa corresponding to the lineage of a species
    :param taxa_to_codes: a dictionary where keys are assigned_taxa and values are taxids
    :param ref_species: a dictionary where keys are assigned_taxa and values are
    reference species in the assigned_taxa
    :return: a list of reference -> assigned_taxa assignments
    """
    # automatic assignment of reference species to sub-branches of overhead assigned_taxa
    if auto == 1:

        # dictionary mapping taxids to the reference species contained within them
        code_to_ref_species = {}
        for ref_org in lineage_dict:
            for code in lineage_dict[ref_org]:
                if code not in code_to_ref_species:
                    code_to_ref_species[code] = [ref_org]
                else:
                    code_to_ref_species[code].append(ref_org)

        assignments = []
        # for each overhead taxon, get automatic assignments. concatenate them
        # at the end as a list of entries to BLAST.
        for taxon in taxa_to_codes:
            taxid = taxa_to_codes[taxon][0]

            # if no reference species are a part of the taxon
            if taxid not in code_to_ref_species:

                i = 1
                parent_id = taxa_to_codes[taxon][i]
                while parent_id not in code_to_ref_species:
                    i += 1
                    parent_id = taxa_to_codes[taxon][i]

                # find the ref species that's closest and assign it to the taxon
                assignments += [[code_to_ref_species[parent_id][0], taxid]]

            else:
                assignments += auto_assign_taxa_to_ref(lineage_dict,
                                                       code_to_ref_species[taxid],
                                                       taxid)

        return assignments

    # user themselves inputs a file, where it's reference sequence + taxon for
    # name
    else:
        # ex. r"C:\Users\tonyx\Downloads\concatenate_species.txt"
        correct_file = False
        assignments_arr = []
        while not correct_file:
            assignments = []
            assignments_file = input("Please enter a path to valid taxon assignment file: ")
            if not os.path.isfile(assignments_file):
                print("Invalid. File does not exist.")
            else:
                correct_file = True
                assignments_arr = file_to_list(assignments_file)
                for line in assignments_arr:
                    split_line = line.strip().split(",")
                    if len(split_line) != 2 or not split_line[1].isalnum():
                        print("Invalid. File format is incorrect")
                        correct_file = False
                        break
                    assignments.append(line.split(","))

        return assignments
