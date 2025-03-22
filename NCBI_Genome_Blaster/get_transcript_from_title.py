import os

if __name__ == "__main__":

    holder = []
    reference_dir = "/Users/tonyx/Documents/chang_lab/good_elasmo_references"
    for gene in os.listdir(reference_dir):

        if gene != ".DS_Store":
            gene_path = os.path.join(reference_dir, gene)

            for taxon in os.listdir(gene_path):

                if taxon != ".DS_Store":
                    taxon_path = os.path.join(gene_path, taxon)

                    for species in os.listdir(taxon_path):

                        if species != ".DS_Store":
                            species_path = os.path.join(taxon_path, species)

                            for transcript_file in os.listdir(species_path):
                                if transcript_file != ".DS_Store":
                                    transcript_file = os.path.splitext(transcript_file)[0]
                                    transcript_file = transcript_file.split("_")[1] + "_" + transcript_file.split("_")[2]
                                    holder.append(transcript_file)
                                    break

    print(holder)