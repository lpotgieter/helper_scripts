# SASUF workshop
*Author: Lizel Potgieter*
*[SLU Bioinformatics Infrastructure](slubi.se)*

These commands are all relative to our Even server. The tools can be installed on any Ubuntu machine and run without a problem. With this protocol, you can annotate effectors in any fungal genome and the immune genes in any plant genome, narrowing down candidate targets for breeding efforts!

Here it is also important to note that there are many other proteins that are involved in host-pathogen interactions that may elucidate a plant immune defense response. This is a *very* short workshop, so we can't look at all the options. Other options to consider are glycoside hydrolases, pectinases, and lipases. We will have a look at protein annotation as well, so you can look at functions on your own.

## Folder Overview
### cercsospora
[Source](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002742065.1/)

|  File |  Description |
|---|---|
|  GCF_002742065.1_CB0940_V2_genomic.fna |  Genomic fasta |  
| genomic.gff  | GFF indicating genes, intergeneic regions, introns, etc   | 
| protein.faa   |  Protein sequences |  
| rna.fna  |  RNA sequences |
| cerco_interpro  |  InterProScan results |    

### beta
[Source](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026745355.1/)

|  File |  Description |
|---|---|
|  GCF_026745355.1_EL10.2_genomic.fna |  Genomic fasta |  
| genomic.gff  | GFF indicating genes, intergeneic regions, introns, etc   | 
| protein.faa   |  Protein sequences |  
| rna.fna  |  RNA sequences | 
| beta_interpro  |  InterProScan results | 

### users
Here you will create your own working directories called by your name with the command:
~~~~
mkdir users/your_name
~~~~
For all of your sanities, please **DO NOT** called your directory *your_name* or I will delete it :)

To change to your directory
~~~~
cd /users/your_name
~~~~

## Workflow Overview
### *Cercospora beticola*
Let us see how many sequences we have in the cercospora/protein.faa file. We will collect all of this info later on, but it is cool to see how your dataset shrinks as you filter!
~~~~
grep ">" ../../cercospora/protein.faa | wc -l
~~~~

We first need to find the sequences with signal peptide regions.

~~~~
signalp6 --fastafile ../../cercospora/protein.faa  --organism eukarya --output_dir cerco_signalp --model_dir /bioinfo/signalp6_fast/signalp-6-package/models
~~~~

We'll just put the output in our working directory for now to make it a bit easier, and then count how many proteins have signal peptide regions
~~~~
cp cerco_signalp/processed_entries.fasta .
grep ">" processed_entries.fasta | wc -l
~~~~

To check whether there are transmembrane domains, we can run Phobius. Phobius identifies which region of the peptide is cytoplasmic, and which is non-cytoplasmic facing. The short report is easier to parse, but we will create both
~~~~
perl ../../phobius/phobius.pl processed_entries.fasta -short > phobius_output_short
~~~~

~~~~
perl ../../phobius/phobius.pl processed_entries.fasta > phobius_output_long
~~~~

We will then predict which proteins are potentially effectors
~~~~
python /bioinfo/EffectorP-3.0/EffectorP.py -f -i processed_entries.fasta -E cerco_effectors.fasta > cerco_effectorp
~~~~

We will then use DeepTMHMM, 
~~~~
biolib run DTU/DeepTMHMM --fasta  cerco_effectors.fasta
~~~~

We will also run Deeploc2 to further determine localisation
~~~~
deeploc2 --fasta cerco_effectors.fasta --output cerco_effectors_deeploc
~~~~

And then InterProScan (**Note: this has been done!** Currently I am still having permission errors when not running as sudo user- it is a work in progress)
~~~~
java -jar /bioinfo/interproscan-5.67-99.0/interproscan-5.jar -i protein.faa -o cerco_effectors_interpro
~~~~


We will download all of these files adn combine them in R to some more visual representations of our data.

We have so much redundancy in determining where things are located, whether they have transmembrane helices and so on since these are all predictions, and some programs will make a positive prediction while 3 others will not. This exercise just helps us narrow down which effectors are the ones we should be focusing on in the lab. Bioinformatics is cool, but real life validation is still **very** important!

## *Beta vulgaris subsp. vulgaris*
We are going to be looking for  NB-ARC domain (PF00931) domains within the dataset.

I completed the InterProScan analysis prior to the start of the workshop, as with the Cercospora dataset. Since we are looking for the particular domain PF00931 we pull all sequences annotated as such out with 
~~~~
grep "PF00931" ../../beta/beta_interpro/protein.faa.gff3 > pf00931
~~~~
In this file we have all lines that contain the pattern pf00931. To extract the sequences that have this identity, we first have to extract the first column of the extracted file
~~~~
awk -F'\t' '{print $1}' pf00931 > pf00931_names
~~~~
Then we use samtools to pull out the sequences from the original fasta by name
~~~~
xargs samtools faidx ../../beta/protein.faa < pf00931_names > pf00931.fasta
~~~~

## What you would do if you had multiple assemblies of each species and wanted to see what is unique to each genome
You would follow this exact pipeline, and then determine the singletons contained within each genome. One way of going about this would be to use OrthoMCL. A simpler way would be to download the subsets of each dataset, and input it into OrthoVenn3 and extract the sequences that are unique to each assembly. This is a great way to determine whether particular genes have been gained/lost that may contribute to pathogenicity/immunity

## Combining this with population genomics data
If you have WGS of a population, you can use the reference genome you have annotated to map reads from each individual to. You can compute various population genetics statistics like Tajima's D and Fst. You can then merge these statistics with the annotation to see which genes are under selection, for instance.
