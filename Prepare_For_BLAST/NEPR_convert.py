import os
import shutil
from Basic_Tools.lists_and_files import make_unique_directory
from Prepare_For_BLAST.get_longest_transcript import get_longest_transcript

def convert_NEPR_directory(dir_path, save_path):
    """
    Convert a directory in NEPR format to BR (blast reference) format.

    :param dir_path: Path to the directory in NEPR format
    :param save_path: Save path
    """

    # general directory creation
    converted_path = make_unique_directory(save_path, "converted_NEPR")

    for gene in os.listdir(dir_path):

        old_gene_path = os.path.join(dir_path, gene)

        # gene directory creation
        os.mkdir(os.path.join(converted_path, gene))
        new_gene_dir = os.path.join(converted_path, gene)

        for taxon in os.listdir(old_gene_path):

            old_taxon_path = os.path.join(old_gene_path, taxon)

            for species in os.listdir(old_taxon_path):

                old_species_path = os.path.join(old_taxon_path, species)

                picked_file = get_longest_transcript(old_species_path)

                if picked_file != "":
                    # copy transcript file into gene directory and rename by its species
                    shutil.copyfile(os.path.join(old_species_path, picked_file),
                                    os.path.join(new_gene_dir, species + ".fas"))

if __name__ == "__main__":
    dir_path = "/Users/tonyx/Documents/chang_lab/NCBI_exon_pull_results24"
    save_path = "/Users/tonyx/Documents/chang_lab"
    convert_NEPR_directory(dir_path, save_path)



