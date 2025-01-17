import os
mysticetes = ["Erob",
              "Mnov",
              "Bacu",
              "Bmus",
              "Egla",
              "Bric",
              "Bphy",
              "Ejap",
              "Cmar",
              "Bbon"
              ]

if __name__ == "__main__":

    alignment_dir = "/Users/tonyx/Documents/chang_lab/1e-5_top_5_gone"
    save_dir = "/Users/tonyx/Documents/chang_lab/1e-5_final"

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
            if line.strip()[1:] in mysticetes:
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
        f1 = open(path1, "r")

        line = f1.readline()
        while line != "":
            while line != "" and line[0] != ">":
                line = f1.readline()
            if line.strip()[1:] not in mysticetes:
                f2.write(line)
                line = f1.readline()
                while line != "" and line[0] != ">":
                    f2.write(line)
                    line = f1.readline()
            else:
                line = f1.readline()
                while line != "" and line[0] != ">":
                    line = f1.readline()



