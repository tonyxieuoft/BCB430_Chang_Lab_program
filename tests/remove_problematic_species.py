import os

if __name__ == "__main__":
    red_problems = ["Pbla", "Pgan", "Zcav"]
    orange_problems = ["Pbla", "Pgan", "Zcav", "Pmin", "Bbon", "Ejap", "Cmar"]
    middle = ["Pbla", "Pgan", "Zcav", "Pmin", "Bbon"]

    alignment_dir = r"/Users/tonyx/Documents/chang_lab/1e-5_name_change"
    save_dir = r"/Users/tonyx/Documents/chang_lab/1e-5_top_5_gone"

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
            if line.strip()[1:] not in middle:
                f2.write(line)
                line = f1.readline()
                while line != "" and line[0] != ">":
                    f2.write(line)
                    line = f1.readline()
            else:
                line = f1.readline()
                while line != "" and line[0] != ">":
                    line = f1.readline()





