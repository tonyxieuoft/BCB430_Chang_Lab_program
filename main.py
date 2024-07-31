import sys

from User_Interaction.UI import UI

if __name__ == "__main__":

    menu = UI()
    print("Welcome to the pipeline!")

    menu.email_view()
    menu.base_directory_view()

    while True:

        menu_choice = menu.main_view()

        if menu_choice == 1:
            menu.exon_puller_view()
        elif menu_choice == 2:
            menu.gene_description_refiner_view()
        elif menu_choice == 3:
            menu.genome_blaster_view()
        elif menu_choice == 4:
            menu.alignment_creator_view()
        elif menu_choice == 5:
            menu.discard_and_reiterate_view()
        else:
            sys.exit(0)


