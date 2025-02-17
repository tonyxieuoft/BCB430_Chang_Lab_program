import os

from Basic_Tools.lists_and_files import file_to_list

if __name__ == "__main__":

    alignment_dir = "/Users/tonyx/Documents/chang_lab/alignments30"
    save_dir = "/Users/tonyx/Documents/chang_lab/ordered_elasmobranchs3"

    order = file_to_list("/Users/tonyx/Documents/chang_lab/Elasmobranch_species_in_order.txt")

    os.mkdir(save_dir)

    for file in os.listdir(alignment_dir):
        path1 = os.path.join(alignment_dir, file)
        path2 = os.path.join(save_dir, file)

        f2 = open(path2, "w")

        for species in order:

            f1 = open(path1, "r")

            line = f1.readline()
            while line != "":
                while line != "" and line[0] != ">":
                    line = f1.readline()

                print(species)
                print(line.split()[1] + " " + line.split()[2])
                if (line.split()[1] + " " + line.split()[2]) == species:
                    f2.write(line)
                    line = f1.readline()
                    while line != "" and line[0] != ">":
                        f2.write(line)
                        line = f1.readline()
                else:
                    line = f1.readline()
                    while line != "" and line[0] != ">":
                        line = f1.readline()

            f1.close()




