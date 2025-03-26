import os

if __name__ == "__main__":

    blast_results_dir = input("Enter transcriptome blast results dir: ")
    pulled_sequence_dir = input("Enter an output directory (not made yet!) for pulled rnaseq full transcripts: ")
    blast_db = input("Enter path to trinity blast database: ")
    os.mkdir(pulled_sequence_dir)

    for blast_file in os.listdir(blast_results_dir):

        pulled_seq_path = os.path.join(pulled_sequence_dir, blast_file)
        blast_filepath = os.path.join(blast_results_dir, blast_file)
        blast_f = open(blast_filepath, "r")

        line = blast_f.readline()
        while (len(line.split()) == 0 or
               (line.split()[0] != "Sequences" or
                line.split()[1] != "producing")):
            line = blast_f.readline()

        line = blast_f.readline()
        line = blast_f.readline()

        if len(line.strip()) > 0:
            transcript_name = line.split()[0]
            e_value = float(line.split()[-1])

        while len(line.strip()) > 0 and e_value < 1e-75:

            print(e_value)

            os.system("blastdbcmd -db " + blast_db + " -entry " + transcript_name + " >> " +
                      pulled_seq_path)

            line = blast_f.readline()

            if len(line.strip()) > 0:
                transcript_name = line.split()[0]
                e_value = float(line.split()[-1])

        #/Users/tonyx/Documents/chang_lab/Good_Shark_References 2



