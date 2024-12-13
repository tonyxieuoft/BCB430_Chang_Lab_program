import os

if __name__ == "__main__":
    red_problems = ["Pbla", "Pgan", "Zcav"]
    orange_problems = ["Pbla", "Pgan", "Zcav", "Pmin", "Bbon", "Ejap", "Cmar"]
    middle = ["Pbla", "Pgan", "Zcav", "Pmin", "Bbon"]

    alignment_dir = r"C:\Users\tonyx\Downloads\improved_results_eval_1e-5"
    save_dir = r"C:\Users\tonyx\Downloads\plus_orange_removed_1e-5"

    for file in os.listdir(alignment_dir):

        path1 = os.path.join(alignment_dir, file)
        path2 = os.path.join(save_dir, file)
        f1 = open(path1, "r")
        f2 = open(path2, "w")

        line = f1.readline()
        while line != "":
            if line.strip()[1:] not in orange_problems:
                f2.write(line)
                line = f1.readline()
                f2.write(line)
            else:
                line = f1.readline()
            line = f1.readline()





