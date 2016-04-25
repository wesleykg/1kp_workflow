# Workflow for 1kp Assemblies

### Prepare Query Sequences
1. Download annotated genome and convert names to the following form:
 * CODE_geneName
2. Split each gene into its own file

### Prepare Local Blast Database
1. Download wanted 1kp SOAPdenovoTrans assembly
2. Convert the assembly into a blast database

### Search
1. Use blastn to find the gene from annotated genome in the 1kp assembly
2. Create different versions of search that can be modified from the command
line

#### Codes

* For 1kp assemblies, use assigned 4-letter code
* For annotated genomes, use first two letters of genus and first two letters of
species