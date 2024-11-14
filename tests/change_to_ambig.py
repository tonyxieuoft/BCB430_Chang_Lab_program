import os

if __name__ == "__main__":

    dir_path = r"C:\Users\tonyx\Downloads\alignments11"
    save_path = r"C:\Users\tonyx\Downloads\cetacea_results_dolphin_ref"

    for alignment in os.listdir(dir_path):

        align_path = os.path.join(dir_path, alignment)
        f = open(align_path, "r").readlines()

        dest = open(os.path.join(save_path, alignment), "w")

        for line in f:
            if line[0] == ">":
                dest.write(line)
            else:
                to_write = ""
                gap_count = 0
                for i in range(len(line)):
                    if i == len(line) - 1 and gap_count > 0:
                        to_write += "N"*(gap_count + 1)
                    elif line[i] == "-":
                        gap_count += 1
                    else:
                        if gap_count > 0:
                            to_write += "N"*gap_count
                        else:
                            to_write += "-"*gap_count
                        gap_count = 0
                        to_write += line[i]
                dest.write(to_write)



