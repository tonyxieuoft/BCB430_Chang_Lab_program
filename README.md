# NCBI Exon Extraction and BLAST Pipeline

## Requirements

The following must be installed:
- Google Chrome
- ChromeDriver (https://sites.google.com/chromium.org/driver/downloads/)
- Selenium for Python (https://pypi.org/project/selenium/)
- BioPython (https://biopython.org/wiki/Download)

Selenium and BioPython are both Python libraries, and can be installed from the command line using Python's package management tool `pip`. More information is available in the above links. 

Google Chrome, ChromeDriver and Selenium are required for the web-driver based automation of NCBI BLAST, whereas BioPython is required to pull sequences from the NCBI Entrez database. 

## Installation and Usage

1. Ensure the requirements are met.
2. Download the code. This can be achieved by cloning the repository via the following command:
```
git clone https://github.com/tonyxieuoft/NCBI_Gene_Extraction_Pipeline.git
```
3. `cd` into the program directory `NCBI_Gene_Extraction_Pipeline`.  
4. Run the main program via the command `python3 main.py`

## Program Overview

The program automates the extraction of gene sequences across different organisms for molecular evolutionary analysis. It offers functionality for two major use cases:
- Pulling exons for well-annotated genes and taxa from the NCBI Gene database
- Using NCBI BLAST to predict and extract gene sequences from non-well annotated genomes

After running the program, it first asks the user for their email, which is required for NCBI Entrez and acts as a point of contact if any issues arise. The user will also be prompted specify the full path to a directory to which the program can download files. For instance, if the user is running the program on Windows 10 and wishes the program to download files to their Downloads folder, they must specify the path `C:\Users\{user}\Downloads`. 

After the user enters in basic details, a main menu with five options appears:
1. Pull existing sequence data from the NCBI Gene database
2. Refine gene names and descriptions used to query the NCBI Gene database.
3. Run NCBI BLAST to pull exons from whole GENOMES
4. Concatenate gene sequences into alignment files
5. Quit the program

The rationale, required input and resultant output for options 1-4 are stated below. 

### Option 1. Pull existing sequence data from the NCBI Gene database

As its name suggests, the NCBI exon puller executes the major use case of pulling experimentally-derived or predicted sequences available on the NCBI Gene database.

#### Required Input

##### Gene Query File

The user will first be prompted by the program to provide a path to a file containing **gene names** and **descriptions** to query with. The file must be *tab-delimited*, in *.txt* format, and organized the following way:
```
marker:gene1query1   marker:gene1query2   ...
marker:gene2query1   marker:gene2query2   ...
marker:gene3query1   marker:gene3query2   ...
...
```
where all queries for a gene are on the same line, separated by tabs. The 

Before each query is a **marker** denoting the *type of query*. Markers are separated from the queries they denote via colons (:). The following markers are available:
- `g` : indicates that a query is an abbreviated gene name (eg. 'RHO', 'GRK7').  
- `d` : indicates that a query is a gene description. (eg. 'rhodopsin', 'G protein-coupled receptor kinase 7')
 For genes that are known under multiple possible abbreviated names or description, multiple instances of the same marker can be used in a given line.

**An Example.** If a user wished to extract gene sequences for 'rhodopsin' and 'G protein-coupled receptor kinase 7', they could use the following gene query file:
```
g:rho\td:rhodopsin
g:grk7\td:"G protein-coupled receptor kinase 7"
```

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





 
