import os

from Basic_Tools.lists_and_files import file_to_list

if __name__ == "__main__":

    path = "/Users/tonyx/Documents/chang_lab/major_comparisons/GNB2_giga_fixed_mafft2.fas"
    path2 = "/Users/tonyx/Documents/chang_lab/major_comparisons/GNB2_giga_fixed_mafft2_organized.fas"
    f2 = open(path2, "w")

    arr = file_to_list(path)

    curr_seq = ""
    name = ""
    dct = {}
    for line in arr:
        if line[0] == ">":

            if curr_seq != "":
                if name not in dct:
                    dct[name] = [curr_seq]
                else:
                    dct[name].append(curr_seq)

            name = line.split()[1] + line.split()[2]
            curr_seq = ""

        else:
            curr_seq += line

    dct[name].append(curr_seq)

    for name in dct:
        counter = 0
        for seq in dct[name]:
            f2.write(">" + name + str(counter) + "\n")
            counter += 1
            f2.write(seq + "\n")




