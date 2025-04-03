import os

if __name__ == "__main__":

    genome_fasta_dir = input("Enter path to genome fastas directory")

    for file in os.listdir(genome_fasta_dir):
        genome_fasta_path = os.path.join(genome_fasta_dir, file)
        os.system("nice -8 busco -i " + genome_fasta_path +
                   " -m genome -c 8 -l vertebrata_odb12 -o BUSCO_output/" + file)
