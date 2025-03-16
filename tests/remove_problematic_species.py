import os

if __name__ == "__main__":

    red_problems = ["Malacoraja senta",
                    "Torpedo suessii",
                    "Squatina squatina"]

    alignment_dir = r"/Users/tonyx/Documents/chang_lab/alignments7"
    save_dir = "/Users/tonyx/Documents/chang_lab/1-1-4-1-e0.5-exon-raw-score-no_force"

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
            if line.split()[1] + " " + line.split()[2] not in red_problems:
                f2.write(">" + line.split()[1][0] + line.split()[2][:5] + "\n")
                line = f1.readline()
                while line != "" and line[0] != ">":
                    f2.write(line)
                    line = f1.readline()
            else:
                line = f1.readline()
                while line != "" and line[0] != ">":
                    line = f1.readline()





