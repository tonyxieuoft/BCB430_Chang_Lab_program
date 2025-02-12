import os

if __name__ == "__main__":

    dir_path = r"/Users/tonyx/Documents/chang_lab/alignments28"
    save_path = r"/Users/tonyx/Documents/chang_lab/no_gaps"

    os.mkdir(save_path)

    for alignment in os.listdir(dir_path):

        align_path = os.path.join(dir_path, alignment)
        f = open(align_path, "r").readlines()

        dest = open(os.path.join(save_path, alignment), "w")

        for line in f:
            if line[0] == ">":
                dest.write(line)
            else:
                to_write = ""
                for i in range(len(line)):
                    if line[i] != "-":
                        to_write += line[i]
                dest.write(to_write)