from Basic_Tools.numeric_user_input import numeric_user_input
from User_Interaction.user_exon_pulling import enter_gene_filepath


class ExonPullerView:

    def __init__(self, ui):
        self.ui = ui

    def gene_filepath_view(self):

        discard_gene_query = 1
        if self.ui.refined_gene_query_filepath != "":
            print()
            print("Path of last gene query file used: " + self.ui.refined_gene_query_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_gene_query = numeric_user_input(1, 2, "")
            if discard_gene_query == 1:
                self.ui.gene_query_filepath = self.ui.refined_gene_query_filepath
        elif self.gene_query_filepath != "":
            print()
            print("Path of last gene query file used: " + self.gene_query_filepath)
            print("Enter 1 to use this file, or enter 2 to input the path of a different one.")
            discard_gene_query = numeric_user_input(1, 2, "")
        if self.gene_query_filepath == "" or discard_gene_query == 2:
            self.gene_query_filepath = enter_gene_filepath()

