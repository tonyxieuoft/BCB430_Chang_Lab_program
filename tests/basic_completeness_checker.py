import os

from Basic_Tools.lists_and_files import file_to_list

if __name__ == "__main__":

    dir_path = r"/Users/tonyx/Documents/chang_lab/ordered_elasmobranchs3"
    save_path = r"/Users/tonyx/Documents/chang_lab/80coverage_elasmobranchs"

    os.mkdir(save_path)

    for alignment in os.listdir(dir_path):

        align_path = os.path.join(dir_path, alignment)
        f = file_to_list(align_path)

        print(alignment)

        dest = open(os.path.join(save_path, alignment), "w")

        i = 0
        while i < len(f):

            line = f[i]

            if line[0] == ">":

                header = line
                i += 1
                line = f[i]

                correct = 0
                for j in range(len(line)):
                    if line[j] != "-":
                        correct += 1

                if (correct / len(line)) > 0.8:
                    dest.write(header + "\n" + line + "\n")
                else:
                    name = header.split()[1] + " " + header.split()[2]
                    print(name)


            i += 1
