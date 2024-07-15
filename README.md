# NCBI Exon Extraction and BLAST Pipeline

## Requirements

The following must be installed:
- Google Chrome
- ChromeDriver (https://sites.google.com/chromium.org/driver/downloads/)
- Selenium for Python (https://pypi.org/project/selenium/)
- BioPython (https://biopython.org/wiki/Download)

More specifically, BioPython is required to access the NCBI Entrez API, and Selenium is required for the web-driver based automation of NCBI BLAST.

## Usage

Ensure the requirements are met, clone the repository and run the main program. 

## Program Overview

The program offers automated functionality for two major use cases:
- Pulling exons for well-annotated genes and taxa from the NCBI Gene database
- Blasting gene exons of reference species against non-well annotated genomes to obtain full sets of sequences for given taxa.

When running the main program, the following will occur sequentially:
1) The program first asks the user for basic information (username, directory to download to)
2) (skippable) Given user-inputted lists of genes and taxa, the program pulls out all available exon information from NCBI and directs it to an output file. Contact with the NCBI server is achieved through the Entrez API. 
3) The program organizes seuqences pulled out in step 2 (or during a previous iteration of the program) into query files in preparation for BLAST.
4) Using the query files assembled in step 3, the program runs the NCBI BLAST program against genomes in the NCBI database for specified taxa. To do so, the program utilizes a Selenium-based webdriver to emulate a web user. 
5) Finally, the program concatenates the results and outputs alignments by gene.

Details for steps 2-5 are provided below. 

## Pulling exons from NCBI

The user provides two files on hand for this part of the program: one that lists taxa to pull sequences for, and another that lists queries that the program will use to search for specific genes. 

### Taxa file

The taxa file should be in .txt format, and must be organized the following way:
```
taxon1
taxon2
taxon3
...
```
Note that there are no headers. Each line contains one taxon. 

### Gene query file

The structure of input for the gene query file is more involved, and the file should be in either .txt or .csv format. As a .txt file, it must be organized the following way:
```
marker: gene1query1, marker: gene1query2,   ...
marker: gene2query1, marker: gene2query2,   ...
marker: gene3query1, marker: gene3:query2,  ...
...
```
where all queries for a gene are on the same line, separated by commas similar to a comma-separated-file format (.csv). Before each query is a marker denoting the type of query. Markers are separated from the queries they denote via colons (':'). 

The following markers are available:
- `g` : indicates that a query is an abbreviated gene name (eg. 'RHO', 'GRK7').  
- `d` : indicates that a query is a gene description. (eg. 'rhodopsin', 'G protein-coupled receptor kinase 7')

 For genes that are known under multiple possible abbreviated names or description, multiple instances of the same marker can be used in a given line.

 ### Output directory

 The results pulled out from the program are outputted in a layer of nested folders with the following structure: `General Folder -> Gene -> Taxon -> Species`. Each species folder contains fasta files that each correspond to a different transcript version. 

## Preparing query files for BLAST.

To prepare query files for BLAST, a folder of sequences mirroring the structure of directories outputted after pulling exons from NCBI must be provided. If the user blasts directly after pulling exons, the output folder of pulled exons will be used to compile the query files for BLAST. 

### Automatic assignment of reference species to sub-taxa

To ensure query sequences are as similar to the subject genome as possible yet altogether cover the entire taxa, the user can select the option for the program to automatically assign available reference query sequences to blast against subject genomes only within a sub-branch of the taxa they are most similar to. How the program does this given an overarching taxon to blast and reference species within that taxon is as follows:
1) Select an arbitary reference species S1 and assign it to the overarching taxon.
2) Select a different reference species S2 and assign it to the largest taxon *T* within its lineage such that *T* is not in the lineage of another already-selected reference species.
3) Repeat step 3 for species S3, S4 ... until all reference species have been assigned a taxon.

The order in which these sub-taxa will be blasted is the reverse order that they were assigned, and species that have been blasted are not blasted again. I'll use the following example to make the algorithm clearer: 
```
Utilizing the program to pull exons for the taxa 'Elasmobranchii', the user has reference sequences from three species:

(1) Carcharodon carcharias
(2) Amblyraja radiata
(3) Hemiscyllium ocellatum

The lineages for each species is as follows:

(1) Carcharodon -> Lamninae -> Alopiidae -> Lamniformes -> Galeoidea -> Galeomorphii -> Selachii -> Elasmobranchii
(2) Amblyraja -> Rajidae -> Rajiformes -> Batoidea -> Elasmobranchii
(3) Hemiscyllium -> Hemiscylliidae > Orectolobiformes -> Galeoidea -> Galeomorphii -> Selachii -> Elasmobranchii

Suppose for step 1 of the algorithm, the program arbirtrarily selects Amblyraja radiata and assigns it the
overarching taxon, 'Elasmobranchii.' For step 2 of the algorithm, it arbitrarily selects Hemiscyllium ocellatum
and assigns it 'Selachii', the largest taxon in its lineage not in the lineage of Amblyraja radiata. Finally,
Carcharodon carcharias is assigned 'Lamniformes', as Elasmobranchii, Selachii, Galeomorphii, and Galoeidea are
in Hemiscyllium ocellatum's lineage (and Elasmobranchii is also in Amblyraja radiata's lineage).

When it's time to BLAST, Carcharodon carcharias sequences are used first to query against Lamniforme genomes.
Then, as the genomes of species that have been already blasted are not blasted again Hemiscyllium ocellatum
sequences are used to query against non-Lamniforme Selachii genomes. Finally, Amblyraja radiata sequences query
non-Lamniforme, non-Selachii (which overlaps) genomes, and are thereby effectively blasted only against Batoidea. 
```
### Manual assignment of reference species to sub-taxa

Much more goes into phylogenetic analysis than purely clade and lineage information, and the algorithm only crudely estimates appropriate reference sequence for a given taxon. If the user is willing to spend more time, they can manually specify these assignments to increase accuracy. The format of the file used to give assignment commands is as follows:
```
reference_species1,sub_taxa1
reference_species2,sub_taxa2
reference_species3,sub_taxa3
...
```

### Filling in missing genes

For speed's sake, pulled sequences for all genes for a given reference species are combined into one query file before BLAST occurs. Sometimes, a reference species will be missing some user-specified genes. In cases where this occurs, the user can manually (or automatically) specify alternative species from which the missing gene sequences can be pulled from and used. No file is required to specify this; instead a response system is built directly into the program when a missing gene is detected.

## Automatic NCBI BLAST

No input files are required for this part. Everything will ahve been handled by the section of the pipeline right above that handles preparing for BLAST. The user can specify custom BLAST parameters, which include:
- expect threshold
- word size

Please note that a Selenium-based webdriver will be used to run this part of the program. A pop-up chrome tab will appear: DO NOT CLOSE IT, unless you wish to terminate the program. THe program emulates a web-user, and accesses NCBI BLAST directly from a web browser. When the program runs to completion, all downloaded files and opened tabs will be automatically closed. 

TODO: accessing BLAST via CLOUD services or locally, after eukaryotic genomes have been downloaded to the local server.

## After BLAST

No input files are required for this section. Here, results from BLAST are automatically concatenated into gene alignments. 

TODO: automate alignment algorithms like ClustalW





 
