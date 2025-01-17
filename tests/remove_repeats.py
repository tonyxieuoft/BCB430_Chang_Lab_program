import os

duplicates = ["Lagenorhynchus obliquidens",
              "Physeter catodon",
              "Lagenorhynchus acutus"]

if __name__ == "__main__":

    alignment_dir = "/Users/tonyx/Documents/chang_lab/alignments18"
    save_dir = "/Users/tonyx/Documents/chang_lab/1e-5_remove_dup"

    os.mkdir(save_dir)

    for file in os.listdir(alignment_dir):
        path1 = os.path.join(alignment_dir, file)
        path2 = os.path.join(save_dir, file)

        f1 = open(path1, "r")
        f2 = open(path2, "w")

        line = f1.readline()
        while line != "":
            while line != "" and line[0] != ">":
                line = f1.readline()
            if line.split()[1] + " " + line.split()[2] not in duplicates:
                f2.write(line)
                line = f1.readline()
                while line != "" and line[0] != ">":
                    f2.write(line)
                    line = f1.readline()
            else:
                line = f1.readline()
                while line != "" and line[0] != ">":
                    line = f1.readline()