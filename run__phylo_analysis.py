from Basic_Tools.numeric_user_input import numeric_user_input
from User_Interaction.file_acceptors import get_generic_directory
from test2.assess_similarity import Analyser

if __name__ == "__main__":

    print("Enter directory containing fasta files of unaligned sequences")
    dir_path = get_generic_directory()
    # /crun/tony.xie/Downloads/phylo_accel/GeMoMa_vision_results

    print("Enter reference sequence converted_NEPR dir")
    reference_path = get_generic_directory()
    # /crun/tony.xie/Downloads/official_results/converted_NEPR

    # just play around in
    wd = "/crun/tony.xie/Downloads/phylo_accel"
    analyser = Analyser("/crun/tony.xie/Downloads/official_results/alignments8", wd)
    # exon2 = Analyser(alignment_path3, wd)

    # full._align_against(exon1, csv_for_r)
    analyser.phylo_tree_generation(reference_path)
    # _global_aligner("bruhasdfaskjfewj", "asdfasifew")