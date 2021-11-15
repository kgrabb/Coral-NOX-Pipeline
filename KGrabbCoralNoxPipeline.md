# Coral NOX Bioinformatic Pipeline

**Author: Kalina Grabb**
*Created November, 2021*

## 

## Read me

This pipeline describes the steps to analyze coral genomes for the presence of NADPH Oxidase (NOX) enzymes. This pipeline requires Python3 v 3.6.5, biopython v1.77, and Bash.

This document is the accompanying methods for a paper currently in progress.

## 

## 

## Outline

1. [Download sequences](#1-Download-sequences) 
   1. [Download coral protein sequences](#1i-Download-coral-protein-sequences)
   2. [Download NOX reference sequences](#1ii-Download-NOX-reference-sequences)
   3. [Download transcriptomic sequences and translate to proteins](#1iii-Download-transcriptomic-sequences-and-translate-to-proteins)
2. [Query annotated coral sequences](#2-Query-annotated-coral-sequences)
   1. [Query annotated sequences for NOX proteins](#2i-Query-annotated-sequences-for-NOX-proteins)
3. [Identify NOX-like coral protein sequences in unannotated sequences](#3-Identify-NOX-like-Coral-protein-sequences-in-unannotated-sequences)
   1. [Create protein family (pfam) profiles using hidden markov models (hmm)  for each NOX-type ](#3i-Create-protein-family-(pfam)-profiles-using-hidden-markov-models-(hmm)-for-each-NOX-type )
   2. [Query unannotated coral protein sequences using hmmSearch](#3ii-Query-unannotated-coral-protein-sequences-using-hmmSearch)
   3. [Identify NOX-like coral sequences](#3iii-Identify-NOX-like-coral-sequences)
   4. [Verify identified NOX-like sequences against *nr* database](#3iv-Verify-identified-NOX-like-sequences-against-*nr*-database) 
5. [NOX Domain Analysis](#4-NOX-Domain-Analysis)
   1. [Download NOX-human domains from online pfam database](#4i-Download-NOX-human-domains-from-online-pfam-database)
   2. [Identify specific domains in NOX-like sequences](#4ii-Identify-specific-domains-in-NOX-like-sequences)
6. [Create a phylogenetic tree](#51-Create-a-phylogenetic-tree)
   1. [Align all coral NOX-like sequences](#5i-Align-all-coral-NOX-like-sequences)
   2. [Create a tree file](#5ii-Create-a-tree-file)
   3. [Overall bash script to execute for creating a tree file](#5iii-Overall-bash-script-to-execute-for-creating-a-tree-file)
   4. [Visualize tree using iTol Software](#5iv-Visualize-tree-using-iTol-Software)



## 

## 1. Download sequences 

### 1.i Download coral protein sequences

Coral protein sequences were obtained online from published papers and `reefgenomics.org`. Majority of coral sequences were unannotated, however, a few were annotated.

Majority of the annotated and some of the unannotated coral protein sequences came from [reefgenomics.org][reefgemonics] using the `wget` command in `Bash` command line:

``` bash
wget http://comparative.reefgenomics.org/faa/Coral/Porites_astreoides_peptides_100.final.clstr.faa
```

[reefgenomics]: http://comparative.reefgenomics.org/ "reefgenomics"

All unannotated coral sequences were stored in a file that we called `coralSeqDatabase.fas`.



### 1.ii Download NOX reference sequences

NOX reference sequences were downloaded from previously published NOX references in two papers:

- [Kawahara et al., 2007; *BMC Ecology and Evolution*][kawahara] 

  [kawahara]: https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-7-109	"Kawahara et al., 2007 "

- [Gandara et al., 2017; *BMC Ecology and Evolution][gandara]

  [gandara]: https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-017-0940-0	"Gandara et al., 2017"

All NOX reference sequences were stored in a file that we called `noxAllReference.fas` . We also split up each NOX type into respective files to be used to classify by NOX-type (e.g. all NOX1 reference sequence were stored in `nox1Reference.fas`).



### 1.iii Download transcriptomic sequences and translate to proteins

Coral transcription sequences were downloaded from online sources, including but not limited to [reefgenomics.org][reefgenomics]. They were then translated to peptides using [TransDecoder][transdecoder] v.2.0.1. The following `Bash` command was used to predict coding regions and extract long open reading frames:

``` bash
TransDecoder.LongOrfs -t coralTranscriptome.fa -m 30
```

[transdecoder]: https://github.com/TransDecoder/TransDecoder/wiki "transdecoder"

All unannotated translated coral sequences were added to our `coralSeqDatabase.fas` for use in our queries below



----



## 2. Query annotated coral sequences

### 2.i Query annotated sequences for NOX proteins

After downloading all annotated coral protein sequences, all sequences were queried for words that were associated with NOX and/or DUOX including *DUOX, NOX, NADPH*. This query was completed by querying within the folder that contained all of the protein sequences and the output was written to a text file named `NOXQueryOutputFile.txt`.

```bash
for j in DUOX NOX NADPH; 
	do echo "$j" >> ../NOXQueryOutputFile.txt; 
	for i in *; 
		do echo "$i" >> ../NOXQueryOutputFile.txt; 
		grep -ci $j $i >> ../NOXQueryOutputFile.txt; 
	done; 
done
```



-----



## 3. Identify NOX-like coral protein sequences in unannotated sequences

### 3.i Create protein family (pfam) profiles using hidden markov models (hmm) for each NOX-type 

Pfams are publicly available for many proteins from the [Pfam database][pfam]. However, for NOX, the only NOX pfam available was not specific to NOX type. Therefore, we chose to build our own pfams from the more specific NOX reference sequences that we had obtained.

To create a pfam, the first step is to align the reference sequences that are to make up the pfam. This was one by using multiple sequence alignments (MSAs) using MUSCLE v3.8.31. Here, we created a pfam specific to each NOX-type. As an example, we have written out the process for NOX1, using the NOX reference file `nox1Reference.fas` that we created in [Section 1.ii](#1ii-Download-NOX-reference-sequences). This process was then completed for each NOX-type:

``` bash
muscle -in nox1Reference.fas -out muscleAlignmentNox1.msa
```

Using the multiple sequence alignments, we then built a hidden markov model (HMM) profile for each NOX-type using the following script:

``` bash
hmmbuild noxSeqNox1.hmm muscleAlignmentNox1.msa
hmmpress noxSeqNox1.hmm
```

[pfam]:http://pfam.xfam.org/	"pfam"



### 3.ii Query unannotated coral protein sequences using hmmSearch

Using the created pfams for each NOX-type, the unannotated coral sequences were queried to identify NOX-like sequences for each NOX-type. 

``` bash
hmmsearch --tblout coralNOX1hmmSearchResults.txt noxSeqKNox1.hmm coralProteinSequences.fas
```

Coral NOX-like sequences identified through hmmSearch were filtered using an e-value < 10^-5^ as a cut-off. 

```bash
sed '/^#/d' coralNOX1hmmSearchResults.txt | awk '{if ($5<=10^-5) print}'  > ../coralNOX1hmmSearchResults_10-5.txt
```



### 3.iii Identify NOX-like coral sequences

To verify the NOX-like sequences, the headers from the hmmSearch were extracted and matched with their FASTA sequences:

```bash
cut -d' ' -f1 coralNOX1hmmSearchResults_10-5.txt | sort -u > coralNox1Headers.txt
```

The headers were then matched with their FASTA sequences from our `coralSeqDatabase.fas` using the `SeqIO` package within `biopython` through the following `python` script, which we call `biopythonSeqFasta`. The input filed necessary for this script are the `coralSeqDatabase.fas` created in [Section 1.i](#1i-Download-coral-protein-sequences) and [Section 1.iii](#1iii-Download-transcriptomic-sequences-and-translate-to-proteins), and the `coralNox1Headers.txt` file created in this section. The output file will store the coral sequences for the given NOX-type that corresponds with the headers in the `coralNox1Headers.txt` file. In this example, we call our output file`coralNOX1Sequences.fas` :

``` python
import Bio
from Bio import SeqIO
import pandas as pd
seqDBFile="coralSeqDatabase.fas"
coralHeadersFile="coralNox1Headers.txt"
outputFile="coralNOX1Sequences.fas"

seqs = pd.read_csv(coralHeadersFile)

count=0
total=0
outputHandle = (open(outputFile,"w"))
for record in SeqIO.parse(seqDBFile, "fasta"):
        total = total+1
        if seqs['coralSequence'].str.contains(record.id+"$").any():     
                count = count +1
                SeqIO.write(record,outputHandle,"fasta")
outputHandle.close()
print(str(count) + " records selected out of " + str(total))
```



### 3.iv Verify identified NOX-like sequences against *nr* database

The identified NOX-like coral sequences were verified against the NCBI *nr* protein database to confirm that the top hits were NOX enzymes. The *nr* database was loaded on our institutions high performance cluster and the following analysis was performed in `bash` command line: 

```bash
blastp -query coralNOX1Sequences.fas -db nr -out coralNox1Blastnr.txt -outfmt 6
```

Top hits were also investigated on the [NCBI website][ncbi] to ensure that the proteins were NOX-like. All coral sequences that had top hits with NOX-like proteins were identified as coral NOX-like proteins.

[ncbi]: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome "ncbi"



---



## 4. NOX Domain Analysis

### 4.i Download NOX-human domains from pfam database

Specific NOX domains for each NOX type within humans were identified via published papers such as [Kawahara et la., 2007][kawahara]. Pfam profiles were then obtained online from the [Pfam database][pfam] by downloading the `pfam.seed` file, for example for `PF08030`, which is the *NAD Binding 6* domain:

```bash
wget http://pfam.xfam.org/family/PF08030/alignment/seed/gzipped
```



### 4.ii Identify specific domains in NOX-like sequences

Specific domains in coral NOX-like sequences were identified using the `.seed` files and `hmmSearch`, similar to that described in [Section 3.i][#3i-Create-protein-family-(pfam)-profiles-using-hidden-markov-models-(hmm)-for-each-NOX-type ], [Section 3.ii][#3ii-Query-unannotated-coral-protein-sequences-using-hmmSearch], and [Section 3.iii][3iii-Identify-NOX-like-coral-sequences]. The only exception to the pipeline described in Sections 3.i - 3.iii is in the creation of the HMM profile. Since a `.seed` file was downloaded for the known NOX domains, the following script was used in `bash` to create a HMM profile:

```bash
hmmbuild PF08030.hmm PF08030.seed
hmmpress PF08030.hmm
```



---



## 5. Create a phylogenetic tree

### 5.i Align all coral NOX-like sequences

In order to create a phylogenetic tree, the output FASTA file from [Section 3.iii](#33-Identify-NOX-like-coral-sequences) `coralNOX1Sequences.fas` can be used. In the case where we want to create a tree for all of the coral NOX-like sequences, all of the sequences can be analyzed through the pipeline described in [Section 3](#3-Identify-NOX-like-coral-protein-sequences-in-unannotated-sequences), or the output FASTA files from [Section 3.iii](#33-Identify-NOX-like-coral-sequences) can be combined using the `cat` command:

```bash
cat coralNOX1Sequences.fas coralNOX2Sequences.fas coralNOX3Sequences.fas coralNOX4Sequences.fas coralNOX5Sequences.fas coralDUOXSequences.fas > coralAllSequences.fas
```

We used this technique in order to create a tree for all coral NOX-like sequences. We then aligned the sequences by making a MSA through the `bash` script that we call `MSAScript`:

```bash
module load bio
module load muscle
muscle -in coralAllSequences.fas -out muscleAlignmentCoralAll.msa
```



### 5.ii Create a tree file

Using the MSA file created above, the `AlignIO` package in `Biopython` was used to convert the MSA to a phylip-relaxed format using the following `python` script that we call `biopythonConvertScript`:

```python
import Bio
from Bio import AlignIO
inFile = "muscleAlignmentCoralAll.msa"
outFile = "muscleAlignmentCoralAll.phy"

AlignIO.convert(inFile, "fasta", outFile, "phylip-relaxed")
```

Using this `.phy` file, a tree file was then created also using the `AlignIO` package in `Biopython`.  The tree files were created using the default distance calculator (`identity`), the neighbor joining constructor (`nj`), and bootstrapped to `100`. The trees are output as a PhyloXML file (`.xml`) using the following `python` script that we call `biopythonCreateTreeScript`:

```python
import Bio
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.Consensus import *

inputFile = "muscleAlignmentCoralAll.phy"
outputFile = "treeFileCoralAll.xml"

alignment = AlignIO.read(inputFile, "phylip")

calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')
trees = bootstrap_trees(alignment, 100, constructor)

print(trees)

Phylo.write(trees, outputFile, 'phyloxml')
```



### 5.iii Overall bash script to execute for creating a tree file

As an option to create multiple tree files from different sequences, the following `bash` script can be run to execute the above scripts in one command. The scripts that are being executed within in this script are:

- `biopythonSeqFasta` from [Section 3.iii](#3iii-Identify-NOX-like-coral-sequences)
- `MSAScript` from [Section 5.i](#5i-Align-all-coral-NOX-like-sequences)
- `biopythonConvertScript` in [Section 5.ii](#5ii-Create-a-tree-file)
- `biopythonCreateTreeScript` in [Section 5.ii](#5ii-Create-a-tree-file)

The script we call `treePipelineScript` and can be executed in `bash` command line in the same folder as all of the other scripts and necessary files. All filenames need to be changed in the script files listed above:

```bash
#!/usr/bin/env python

chmod +x biopythonSeqFasta
python ./biopythonSeqFasta

chmod +x MSAScript
bash ./MSAScript

chmod +x biopythonConvertScript
python ./biopythonConvertScript

chmod +x biopythonCreateTreeScript
python ./biopythonCreateTreeScript

```



### 5.iv Visualize tree using iTol Software

To visualize the tree files, we used the [Interactive Tree of Life (iTol)][itol], which is a web-based platform. The PhyloXML file was uploaded into the web platform and amended accordingly to add colors for each NOX-type. 

[itol]: https://itol.embl.de/	"itol"