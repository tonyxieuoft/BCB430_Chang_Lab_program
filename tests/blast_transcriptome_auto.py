import os

if __name__ == "__main__":

    reference_cds_path = input("Enter reference full cds file")
    trinity_path = input("Enter path to trinity")

    output_dir = input("Enter path to a new directory (not yet created to store the output")

    os.mkdir(output_dir)

    for gene in os.listdir(reference_cds_path):
        gene_path = os.path.join(reference_cds_path, gene)

        for taxon in os.listdir(gene_path):
            taxon_path = os.path.join(gene_path, taxon)

            for species in os.listdir(taxon_path):
                species_path = os.path.join(taxon_path, species)

                for transcript in os.listdir(species_path):

                    transcript_path = os.path.join(species_path, transcript)

                    output_path = os.path.join(output_dir, gene + "_" + os.path.splitext(transcript)[0] + ".txt")

                    print(transcript_path)

                    os.system("nice -2 blastn -db '" + trinity_path +
                              "' -evalue 0.05 -word_size 11 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 "
                              "-num_threads 16 -query '" + transcript_path +
                              "' > '" + output_path + "'")


