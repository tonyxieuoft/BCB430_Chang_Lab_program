import os

if __name__ == "__main__":

    in_path = r"C:\Users\tonyx\Downloads\out_hypanus"
    out_path = r"C:\Users\tonyx\Downloads\out_hypanus_trinity_blast"

    for file in os.listdir(in_path):

        in_file_path = os.path.join(in_path, file)

        f1 = open(in_file_path, "r")
        f2 = open(os.path.join(out_path, os.path.splitext(file)[0] + ".txt"), "w")

        f2.write(f1.read())
