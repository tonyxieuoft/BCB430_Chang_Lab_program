import os
import shutil

if __name__ == "__main__":

    og_paths = ["/Users/tonyx/Documents/chang_lab/good_batoid_references",
                "/Users/tonyx/Documents/chang_lab/good_shark_references"]

    combined_dir_path = "/Users/tonyx/Documents/chang_lab/good_elasmo_references"
    os.mkdir(combined_dir_path)

    for og_path in og_paths:
        for gene_dir in os.listdir(og_path):

            combined_gene_path = os.path.join(combined_dir_path, gene_dir)
            if not os.path.isdir(combined_gene_path):
                os.mkdir(combined_gene_path)

            gene_path = os.path.join(og_path, gene_dir)
            for taxon_dir in os.listdir(gene_path):

                taxon_path = os.path.join(gene_path, taxon_dir)
                shutil.copytree(taxon_path, os.path.join(combined_gene_path, taxon_dir))




