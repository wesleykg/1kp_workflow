# Workflow for 1kp Assemblies

## How to Use

* `make download` will download and decompress assemblies and stats files
* `make split` will divide a specified genome in fasta format into
   individual genes
* `make blastdb` creates a blast database for each assembly
* `make search` performs a blast search on all databases using all 
   genes from a specied genome
* `make create` creates an alignment for each gene found in each 
   species
* `make align` automatically aligns the alignment created by 
  `make create`

Specify a genome by typing genome=CODE (See below for CODE details) 
after the above commands

## Roadmap

### Prepare Query Sequences
1. Download annotated genome and convert names to the following form:
   * CODE-geneName
2. Split each gene into its own file

### Prepare Local Blast Database
1. Download & decompress wanted 1kp SOAPdenovoTrans assembly
2. Filter out "bad" scaffolds using Transrate stats file
3. Create a blast database using the filtered 1kp assembly

### Search
1. Use `blastn` and `tblastx` to find a gene in the 1kp assembly from an  
   annotated genome
2. Create different versions of search that can be modified from the command  
   line (ex. blastp, megablast, hmmer3)

### Alignment
1. Collect search hits from each database into an alignment for each gene

___

#### Codes

* For 1kp assemblies, use assigned 4-letter code
* For annotated genomes, use first three letters of genus and first 
  two letters of species
