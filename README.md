# NCBI Exon Extraction and BLAST Pipeline

##### Table of Contents

[Requirements](https://github.com/tonyxieuoft/NCBI_Gene_Extraction_Pipeline/blob/master/README.md#requirements)

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
2. Download the code. If Git is installed (https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), this can be achieved by cloning the repository on the command line with the following command:
```
git clone https://github.com/tonyxieuoft/NCBI_Gene_Extraction_Pipeline.git
```
3. `cd` into the program directory `NCBI_Gene_Extraction_Pipeline`.  
4. Run the main program via the command `python3 main.py`

## Program Overview

The program automates the extraction of gene sequences across different organisms for molecular evolutionary analysis. It offers functionality for two major use cases:
- Pulling exons for well-annotated genes and taxa from the NCBI Gene database
- Using NCBI BLAST to predict and extract gene sequences from non-well annotated genomes

After running the program, the user is first aked for their email, which is required for NCBI Entrez and acts as a point of contact if any issues arise. The user will also be prompted specify the full path to a directory to which the program can download files. For instance, if the user is running the program on the Windows command line and wishes the program to download files to their Downloads folder, they must specify the path `C:\Users\{user}\Downloads`. 

After the user enters in basic details, a main menu with five options appears:
1. Pull existing sequence data from the NCBI Gene database
2. Refine gene names and descriptions used to query the NCBI Gene database.
3. Run NCBI BLAST to pull exons from whole GENOMES
4. Concatenate gene sequences into alignment files
5. Quit the program

As a general rule of thumb, option #1 must be executed first, as all other options (except for 5) depend upon its results. The rationale, required input and resultant output for options 1-4 are stated below. 

### Option 1. Pull existing sequence data from the NCBI Gene database (NCBI Exon Puller)

As its name suggests, the NCBI exon puller automates the process of pulling experimentally-derived or predicted sequences available on the NCBI Gene database. The formats of required input and resultant output files/directories are stated below. 

#### Gene Query File (required input)

The user will first be prompted by the program to provide a path to a file containing **gene names** and **descriptions** to query the NCBI Gene database with. The file must be *tab-delimited*, in *.txt* format, and organized the following way:
```
marker:gene1query1   marker:gene1query2   ...
marker:gene2query1   marker:gene2query2   ...
marker:gene3query1   marker:gene3query2   ...
...
```
where all queries for a gene are on the *same* line, separated by tabs. 

Before each query is a **marker** denoting the *type of query*. Markers are separated from the queries they denote via colons (:). The following markers are available:
- `g` : indicates that a query is an abbreviated gene name (eg. 'RHO', 'GRK7').  
- `d` : indicates that a query is a gene description. (eg. 'rhodopsin', 'G protein-coupled receptor kinase 7')

For genes that are known under multiple possible abbreviated names or description, more than two queries and multiple instances of the same marker can be used in a given line.

**An Example.** If a user wished to extract gene sequences for 'rhodopsin', 'G protein-coupled receptor kinase 7', and 'ATP binding cassette subfamily A member 4' (abbrev. gene name ABCA4), they could use the following gene query file:
```
g:rho     d:rhodopsin
g:grk7    d:"G protein-coupled receptor kinase 7"
g:abca4   d:"ATP binding cassette subfamily A member 4"
```
There is no limit to the number of genes that can be inputted into the gene query file. Note that the file is case-insensitive to query *names* and *descriptions*, but markers must be lowercase. Therefore `g:RHO`, `g:rho`, and `g:Rho` yield the same results, but `G:rho` results in an error.  

The NCBI exon puller filters results based on **exact** matches, so provided gene names and descriptions must be identical to the ones stored in the NCBI database.

#### Taxa File (required input)

After providing the gene query file, the user will additionally be prompted to provide a path to a file containing **taxa of interest** to pull for. Again, the file must be in *.txt* format, with one taxon per line as follows:
```
taxon1
taxon2
taxon3
...
```
Note that the provided names must match a record in NCBI's taxonomy database. 

**An Example.** Suppose the user wanted to pull the genes in the previous example from bats, whales/dolphins/porpoises and humans, they could use the following taxa file:
```
chiroptera
cetacea
Homo sapiens
```
Again, like names and descriptions in the gene query file, the taxon names are *case-insensitive*. 

#### Exons or Full Sequences (required input)

The user is asked whether they would like exons or the full gene sequences to be pulled out. In both cases, only **coding regions** will be returned. Rather than a file, the user simply enters '1' for exons, and '2' for full sequences.

#### Selecting for Optimal Transcript Variants (required input)

*Multiple* transcript variants for the same gene and species may be pulled out in an iteration of the NCBI exon puller. The user will be asked whether they would like to keep all variants, or automatically select for the variant with "optimal length" and **discard the rest**. Read more about this method and its rationale in the *Methods and Algorithms* section. 

#### NCBI Exon Pull Results (output)

The results pulled out from the NCBI Exon Puller are contained within in a layer of nested folders with the following structure: 
```
General Folder -> Gene -> Taxon -> Species
```
Each species folder contains fasta files that each correspond to a different transcript variant.

The directory names at the `Gene` layer are based on the first names/descriptions of each line of the input gene query file. For instance, if the input gene query file was:
```
g:rho                 g:rh1    d:rhodopsin
d:cyclooxygenase-2    g:cox2
```
then the names of the Gene-level directories would be `rho` and `cyclooxygenase-2`. 

The directory names at the `Taxon` layer are the names from the input taxa file, and the directory names at the `Species` layer are the scientific names of species with available transcripts. 

Each *fasta file* is titled `{transcript_length}_{transcript_accession}.fas` and each *fasta heading* includes the gene, scientific name, transcript accession, genome accession, and interval (relative to +1) of a particular exon/full sequence. 

From here on forth, the structure of the directories produced by the NCBI Exon Puller will be referred to as the **NEPR** format.

### Option 2: Refine gene names and descriptions used to query the NCBI Gene database (Gene Description Refiner)

Oftentimes, gene names and descriptions for the same gene vary between different taxa. This inconsistent naming is especially prevelant for taxonomic groups that lack sequence annotation (eg. Elasmobranchii, the group containing sharks, rays and skates). 

For instance, CNGA3, a visual gene involved in the phototransduction cycle in cones, has the gene description `cyclic nucleotide-gated cation channel alpha 3, cone` for bony fish, but in elasmobranches it is labelled `cyclic nucleotide gated channel subunit alpha 3a`. In the same vein, the abbreviated gene name for CNGA3 in elasmobranchs is `CNGA3A`. Although the differences are slight (missing 'cation', addition of 'subunit', one letter addition), the NCBI exon puller filters for **exact** gene name and description matches (to avoid pulling extraneous gene sequences), and would fail to pull out the elasmobranch sequences had the gene description for fish been used. 

To alleviate the burden of manually checking the gene name/description naming conventions for every taxa of interest, the program can automatically identify alternative gene names/descriptions. First, the user must run an iteration of option #1 (the NCBI exon puller). Sequences acquired via the NCBI exon puller then act as queries for **sequence homology searches** of the NCBI **nucleotide (nt)** database. When highly similar transcripts are identified, their gene names and descriptions are extracted and appended to a new, "refined" gene name/description query file.  

The formats of required input and resultant output files/directories are stated below in greater detail:

#### Directory containing query sequences for homology search (required input)

The directory must be in ***NEPR (NCBI exon pull results) format**, and contain sequences for the genes the user wishes to refine names/descriptions for. Ideally, there should be a query sequence for each gene and taxon of interest. 

As mentioned, the easiest way to generate query sequences is to **first run an iteration of option #1** with the same genes and taxa of interest. If a previous iteration of the NCBI Exon Puller has been ran in a program session, the program will automatically suggest the filepath to its results directory as input for the gene description refiner. 

#### Original gene query file (required input)

This file must be in the same format as one used as input by the NCBI exon puller (see above). This file is used as a starting template for the refined gene query file. (below) 

#### Refined gene query file (output)

The refined file produced by the gene description refiner is in the *same format* as the original gene input file. The refined file begins as a *copy* of the original. After identifying alternative names/descriptions for genes in the original gene query file, the program appends these newly extracted names to the end of the gene's line in the refined gene query file. 

After the homology search is complete, the refined file should be used as the input gene query file in a **second** iteration of the NCBI exon puller to acquire sequences that were previously missed.  

### Option 3 part (a): Generating query files in preparation of blasting whole GENOMES (NCBI Genome Blaster)

The main purpose of the NCBI Genome Blaster is predict and extract gene sequences from unannotated whole genomes of species that do not directly have NCBI Gene database sequences.  

BLAST requires the provision of reference sequences to query subject genomes. For the current version of the program, these reference sequences must be in contained in a **NEPR** format directory and acquired from an iteration of the NCBI exon puller.  

Ideally, the phylogenetic distance between the query sequence and subject genome species' should be minimal to maximize the likelihood of accurate coding-region extraction. Therefore, it is unwise to use a single reference sequence to query subject genomes for an entire large taxon; this is especially true for genes under positive selection. Instead, each reference species should be assigned to BLAST subject genomes for a **small** taxon that they are closest to **within** the overarching taxon of interest. 

The program generates query files for BLAST based on the above considerations. The formats of required input and resultant output files/directories are stated below in greater detail:

#### Directory Containing Reference Sequences (required input)

As stated, this directory must be in **NEPR** format. If a previous iteration of the NCBI Exon Puller has ran during the program session, the program automatically suggests the path to its results directory as input for the NCBI Genome Blaster. The genes and taxa of interest within the input reference sequence directory will dictate the genes and taxa of interest whose sequences will be pulled out by the NCBI Genome Blaster. 

For instance, suppose a user specified an initial iteration of the NCBI exon puller (option #1) to pull out rhodopsin sequences for the taxon `Chiroptera`. Afterward, suppose they called option #3 and selected the results directory as the reference sequence directory. Then, the NCBI Genome Blaster would BLAST and extract rhodopsin sequences exclusively from Chiroptera subject genomes that did not have results from the NCBI Exon Puller. 

#### Automatic or Manual Assignment (required input)

The program asks for users to choose between automatic or manual assignment of reference species to sub-taxa within the overarching taxa of interest to BLAST. Automatic assignment assigns a reference species to every subject genome with the goal of maximizing query sequence to subject genome similarity. Read more about this method and its rationale in the *Methods and Algorithms* section. Alternatively, the user can choose to manually select the sub-taxa within the overarching taxa of interest that each reference species will BLAST against.  

#### Manual Assignment File (optional input)

The file used to give manual assignment commands must be in *.csv* or *.txt* format as follows:
```
reference_species1,sub_taxa1
reference_species2,sub_taxa2
reference_species3,sub_taxa3
...
```
The program interprets the file as instructions to blast query sequences for reference_species1 against subject genomes for species in sub_taxa1, then blast query sequences from reference_species2 against subject genomes for species in sub_taxa2, and so on and so forth. *Assignment matters*, in that species whose subject genomes that have been blasted will **not** be blasted again.

The **NCBI taxonomy taxid** of each sub-taxa, **not** their English names, must be provided. They can be obtained by searching the taxonomy browser at this link: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi. In contrast, the scientific names of each reference species must be provided (the names matching their names as printed on the NEPR input directory). 

**An Example.** Suppose the user has ran an iteration of the NCBI exon puller to obtain rhodopsin sequences for Chiroptera. As a result, they've obtained references sequences for three species: Molossus molossus, Desmodus rotundus, and Myotis daubentonii. 

They want to BLAST Molossus molossus query sequences against Molossus genomes, Desmodus rotundus query sequences against Phyllostomidae genomes, and Myotis daubentonii query sequences against all other genomes in Chiroptera. 

The NCBI taxonomy taxids for Molossus, Phyllostomidae, and Chiroptera are 27621, 9415, and 9397, respectively. Then, they can create the following manual assignment file:
```
Molossus molossus,27621
Desmodus rotundus,9415
Myotis daubentonii,9397
```
Based on the file, the program first blasts Molossus molossus query sequences against the taxon 27621 (Molossus). Then, the program blasts Desmodus rotundus query sequences against the taxon 9415 (Phyllostomidae). However, when it finally comes time to blast Myotis daubentonii query sequences against 9397 (Chiroptera), as the subject genomes in Molossus and Phyllostomidae have already been blasted (and Molossus, Phyllostomidae are sub-taxa of Chiroptera), only members of Chiroptera *not* in Molossus and Phyllostomidae will actually be blasted against by Myotis daubentonii. 

**An Alternative Example.** Suppose the user now wishes to blast Myotis daubentonii query againsts against all genomes in Chiroptera first, then blast Molossus molossus query sequences against Molossus genomes. They provide the following manual assignment file:
```
Myotis daubentonii,9397
Molossus molossus,27621
```
The program first blasts Myotis daubentonii query sequences against subject genomes for all of Chiroptera. But when the program executes the next line, as Molossus genomes were already blasted against by Myotis daubentonii, Molossus molossus query sequences are blasted against *no* subject genomes.  

#### Query File Directory (Output)

After automatic or manual assignment is complete, for each reference species, the program concatenates all of its sequences into one query file. Each query file's title is the sub-taxon its reference species is assigned to, and all query files are stored in a single query file directory. If a reference species is missing a particular gene sequence, the gap is filled in by the gene sequence of another reference species that shares the earliest most recent common ancestor.

This output largely serves no function to the user. Instead, it carries intermediate data for option 3(b) to run. 

### Option 3 part (b): Run NCBI BLAST to pull exons from whole GENOMES

As stated previously, the main purpose of the NCBI Genome Blaster is predict and extract gene sequences from un-annotated whole genomes of species that do not directly have NCBI Gene database sequences.

If option 3 part (b) is reached, most of the input work will have already been handled by option 3 part (a), with sub-taxa assignment and query file generation complete. 

The formats of the remaining required input and resultant output files/directories at this step are stated below in greater detail:

#### Remote or Local BLAST (required input)

The user is asked to choose between blasting remotely on the NCBI server or locally using only their computer's resources. 

During the automation of remote BLAST, the program emulates a web user and accesses the NCBI web application to obtain genome accessions and run BLAST processes. At any given moment, approximately ten BLAST processes run parallel to each other. 

In contrast, local BLAST requires the additional installment of the NCBI BLAST+ and NCBI Datasets command line applications. In the local BLAST workflow, subject genomes are sequentially downloaded, blasted against, then immediately deleted one-by-one. The process is slow and laborious, with the main bottleneck being the genome download step. Currently, because of its issues, local BLAST is not supported. Right now, I am working to create a "local server" version of the local BLAST that leverages the increased memory and processing power of local servers to improve upon these weaknesses. 

#### Expect Threshold (optional input)

The user is given the opportunity to enter a custom expect threshold, or stick with the default e-value of 0.05. As the program automatically only selects sequences with the lowest e-values, changing the expect threshold here does not make a significant difference, and may actually decrease BLAST efficiency. 

#### BLAST Results Directory (output)

The BLAST results directory has the same structure as **NEPR** directories. The gene and taxonomic directories from the input NEPR directory in option 3 part (a) (supplying query sequences) to are conserved. All naming conventions for fasta files and fasta headings remain the same. The only exception is the replacement of the "mRNA accession" section in the fasta heading with a "reference query mRNA accession" section, as the new BLAST-predicted sequences do not directly correspond to any RNA transcript sequences in NCBI Entrez.   

### Option 4: Concatenate gene sequences into alignment files

Once gene sequences are obtained from the NCBI Gene database and extracted from unannotated genomes via BLAST, this option offers the capacity to concatenate all sequences from different species for a particular gene into one gene alignment. If the gene sequences are split into exons, the program joins them together in the alignment. 

The formats of the remaining required input and resultant output files/directories at this step are stated below in greater detail:

#### NEPR / BLAST Result Directories (required input)

If the NCBI exon puller and/or Genome Blaster have previously ran during the program session, the program automatically suggests their results as input. Otherwise, the user must manually input paths to directories (in NEPR format) that contain gene sequences they wish to concatenate into alignments. 

#### Alignments Directory (output)

This directory contains alignment files in fasta format for each gene present in the input directory. 

## Methods and Algorithms

### Automatic assignment of reference species to sub-taxa in option 3(a)

To ensure query sequences are as similar to the subject genome as possible yet altogether cover the entire taxa, the user can select the option for the program to automatically assign available reference query sequences to blast against subject genomes only within a sub-branch of the taxa they are most similar to. 

How the program does this given an overarching taxon to blast and reference species within that taxon is as follows:
1) Select an arbitary reference species S1 and assign it to the overarching taxon.
2) Select a different reference species S2 and assign it to the largest taxon *T* within its lineage such that *T* is not in the lineage of another already-selected reference species.
3) Repeat step 3 for species S3, S4 ... until all reference species have been assigned a taxon.

The order in which these sub-taxa will be blasted is the reverse order that they were assigned, and species that have been blasted are not blasted again. 

**An Example.** Utilizing the program to pull exons for the taxa 'Elasmobranchii', suppose the user has reference sequences from three species:

1) Carcharodon carcharias
2) Amblyraja radiata
3) Hemiscyllium ocellatum

The lineages for each species is as follows:

1) Carcharodon -> Lamninae -> Alopiidae -> Lamniformes -> Galeoidea -> Galeomorphii -> Selachii -> Elasmobranchii
2) Amblyraja -> Rajidae -> Rajiformes -> Batoidea -> Elasmobranchii
3) Hemiscyllium -> Hemiscylliidae > Orectolobiformes -> Galeoidea -> Galeomorphii -> Selachii -> Elasmobranchii

Suppose for step 1 of the algorithm, the program arbirtrarily selects Amblyraja radiata and assigns it the overarching taxon, 'Elasmobranchii.' For step 2 of the algorithm, it arbitrarily selects Hemiscyllium ocellatum and assigns it 'Selachii', the largest taxon in its lineage not in the lineage of Amblyraja radiata. Finally, Carcharodon carcharias is assigned 'Lamniformes', as Elasmobranchii, Selachii, Galeomorphii, and Galoeidea are in Hemiscyllium ocellatum's lineage (and Elasmobranchii is also in Amblyraja radiata's lineage).

When it's time to BLAST, Carcharodon carcharias sequences are used first to query against Lamniforme genomes. Then, as the genomes of species that have been already blasted are not blasted again Hemiscyllium ocellatum sequences are used to query against non-Lamniforme Selachii genomes. Finally, Amblyraja radiata sequences query non-Lamniforme, non-Selachii (which overlaps) genomes, and are thereby effectively blasted only against Batoidea. 

**Please note** that much more goes into phylogenetic analysis than purely clade and lineage information, and the algorithm only roughly estimates appropriate reference sequence for a given taxon. If the user is willing to spend more time and has phylogenetic trees with molecuar distances on hand, they can manually specify these assignments to increase accuracy.

### More to come!



 
