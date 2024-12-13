import os

if __name__ == "__main__":

    alignments_path = r"C:\Users\tonyx\Downloads\intermediate1"
    renamed_path = r"C:\Users\tonyx\Downloads\improved_results_eval_1e-8"

    os.mkdir(renamed_path)

    for file in os.listdir(alignments_path):
        f_read_name = os.path.join(alignments_path, file)
        f_write_name = os.path.join(renamed_path, file)

        f_read = open(f_read_name, "r")
        f_write = open(f_write_name, "w")

        line = f_read.readline()
        while line != "":
            if line[0] == ">":
                split_heading = line.split()
                print(split_heading[1] + " " + split_heading[2])
                last_abbrev_length = min(3, len(split_heading[2]))
                f_write.write(">" +
                              split_heading[1][0] +
                              split_heading[2][:last_abbrev_length] +
                              "\n")
            else:
                f_write.write(line)
            line = f_read.readline()

        print("file done")

        f_read.close()
        f_write.close()

