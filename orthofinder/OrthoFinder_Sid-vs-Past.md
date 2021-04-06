---
layout: post
title: Ortholog identification between transcriptome assemblies
date: '2021-03-18'
categories: Protocols
tags: [Transcriptome, Bioinformatics, Ortholog]
---

# Ortholog identification 
## Compare *Porites astreoides* and *Siderastrea sp.* transcriptome assemblies

**Orthologs:**
- genes in different species that evolved from a common ancestral gene
- are thought to retain the same function in the course of evolution

### Orthogroups & Orthologs

**Orthologs** are pairs of genes that descended from a single gene in the last common ancestor (LCA) of two species. 

An orthogroup is the extension of the concept of orthology to groups of species. 
An **orthogroup** is the group of genes descended from a single gene in the LCA of a group of species.

### Why Orthogroups

Orthogroups allow you to analyse all of your data

All of the genes in an orthogroup are descended from a single ancestral gene. Thus, **all the genes in an orthogroup started out with the same sequence and function.** As gene duplication and loss occur frequently in evolution, one-to-one orthologs are rare and limitation of analyses to on-to-one orthologs limits an analysis to a small fraction of the available data. By analysing orthogroups you can analyse all of your data.

##### Orthogroups allow you to define the unit of comparison

It is important to note that with orthogroups you choose where to define the limits of the unit of comparison. For example, if you just chose to analyse human and mouse in the above figure then you would have two orthogroups.
Orthogroups are the only way to identify orthologs.

Orthology is defined by phylogeny. It is not definable by amino acid content, codon bias, GC content or other measures of sequence similarity. Methods that use such scores to define orthologs in the absence of phylogeny can only provide guesses. The only way to be sure that the orthology assignment is correct is by conducting a phylogenetic reconstruction of all genes descended from a single gene the last common ancestor of the species under consideration. This set of genes is an orthogroup. Thus, the only way to define orthology is by analysing orthogroups.


### Method:  OrthoFinder
OrthoFinder is a software program for phylogenetic orthology inference.

[OrthoFinder GitHub](https://github.com/davidemms/OrthoFinder)
[OrthoFinder Tutorials](https://davidemms.github.io/)

[Running example OrthoFinder analysis](https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html)
[OrthoFinder best practices](https://davidemms.github.io/orthofinder_tutorials/orthofinder-best-practices.html)
[Exploring OrthoFinder results](https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html)
[What OrthoFinder provides](https://github.com/davidemms/OrthoFinder#what-orthofinder-provides)
[Orthogroups, orthologs, paralogs](https://github.com/davidemms/OrthoFinder/#orthogroups-orthologs--paralogs)


### References

[Emms, D.M., Kelly, S. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biol 16, 157 (2015)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2)

[Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)

If you use the OrthoFinder species tree then also cite:

[Emms D.M. & Kelly S. STRIDE: Species Tree Root Inference from Gene Duplication Events (2017), Mol Biol Evol 34(12): 3267-3278](https://doi.org/10.1093/molbev/msx259)

[Emms D.M. & Kelly S. STAG: Species Tree Inference from All Genes (2018), bioRxiv https://doi.org/10.1101/267914](https://www.biorxiv.org/content/10.1101/267914v1)


It takes as input the proteomes for the species you want to analyse and from these it automatically:

    - infers the orthogroups for your species
    - infers a complete set of rooted gene trees
    - infers a rooted species tree
    - infers all orthology relationships between the genes using the gene trees
    - infers gene duplication events and cross references them to the corresponding nodes on the gene and species trees
    - provides comparative genomics statistics for your species

OrthoFinder uses gene trees. This means you can check every ortholog relationship in the gene tree it came from. The use of gene trees gives very high ortholog inference accuracy. Despite OrthoFinder using a more rigorous, gene-tree based approach to ortholog inference it is incredibly fast.

The default, and fastest, version of OrthoFinder uses DIAMOND [24] for sequence similarity searches. These sequence similarity scores provide both the raw data for orthogroup inference [10] and for gene tree inference of these orthogroups using DendroBLAST [24]. The default implementation of OrthoFinder has been designed to enable a complete analysis with maximum speed and scalability using only gene sequences as input. However, OrthoFinder has also been designed to allow the use of alternative methods for tree inference and sequence search to accommodate user preferences. For example, BLAST [4] can be used for sequence similarity searches in place of DIAMOND. Similarly, gene trees do not need to be inferred using DendroBLAST. Instead, OrthoFinder can automatically infer multiple sequence alignments and phylogenetic trees using most user-preferred multiple sequence alignment and tree inference methods. Moreover, if the species tree is known prior to the analysis, this can also be provided as input, rather than inferred by OrthoFinder.


----------------------------------------------------------------------------------------------
## Transcriptome assemblies

*Important to use the longest isoform per transcript files (where multiple isoforms have been removed)*

### Porites astreoides
- Hybrid reference
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/hybridref/hybridreference.fasta

- Host (Kenkel) reference, longest isoform with >500 bp length threshold
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/


Just longest isoform assembly, or also filter for >500 bp length?

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/pastreoides_may2014/

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 29422_past_LongestContig.fasta
```

```
The total number of sequences is 29422
The average sequence length is 537
The total number of bases is 15806910
The minimum sequence length is 0
The maximum sequence length is 8171
The N50 is 640
Median Length = 632
contigs < 150bp = 899
contigs >= 500bp = 10714
contigs >= 1000bp = 2911
contigs >= 2000bp = 273
```

#### Add suffix
```
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Past 29422_past_LongestContig.fasta
```

cp 29422_past_LongestContig_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/

#### Filter by length threshold (>500 bp)

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_len_filter.py 500 500lnThresh 29422_past_LongestContig.fasta
```

Output: 
- Number of total seqs for 29422_past_LongestContig.fasta: 29422
- Number of seqs over 500 for 29422_past_LongestContig.fasta: 10715

**Best practice to rename with contig number**

grep -c ">" 29422_past_LongestContig500lnThresh.fasta
10714

mv 29422_past_LongestContig500lnThresh.fasta 10714_past_LongestContig_500ln.fasta


#### Add suffix
```
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Past 10714_past_LongestContig_500ln.fasta
```

```
cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/pastreoides_may2014/10714_past_LongestContig_500ln_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/
```

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 10714_past_LongestContig_500ln_suffixed.fasta
```

```
The total number of sequences is 10714
The average sequence length is 894
The total number of bases is 9585005
The minimum sequence length is 500
The maximum sequence length is 8171
The N50 is 908
Median Length = 510
contigs < 150bp = 0
contigs >= 500bp = 10714
contigs >= 1000bp = 2911
contigs >= 2000bp = 273
```


### Symbiodinium microadriaticum (A1)
- includes ITS2 type A1 symbiont assembly (González-Pech et al. 2019)
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/42652_SymTranscriptome_suffixed.fasta
*All single genes, example: Smic_CassKB8.gene1.mRNA1_Sym*

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 42652_SymTranscriptome_suffixed.fasta
```

```
The total number of sequences is 42652
The average sequence length is 1836
The total number of bases is 78347212
The minimum sequence length is 102
The maximum sequence length is 78486
The N50 is 2727
Median Length = 1836
contigs < 150bp = 4
contigs >= 500bp = 36346
contigs >= 1000bp = 25680
contigs >= 2000bp = 12612
```

#### Filter by length threshold (>500 bp)

```
/cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/fasta_len_filter.py 500 500lnThresh 42652_SymTranscriptome_suffixed.fasta
```

Output:
- Number of total seqs for 42652_SymTranscriptome_suffixed.fasta: 42652
- Number of seqs over 500 for 42652_SymTranscriptome_suffixed.fasta: 36346

**Best practice to rename with contig number**

grep -c ">" 42652_SymTranscriptome_suffixed500lnThresh.fasta
36345

```
mv 42652_SymTranscriptome_suffixed500lnThresh.fasta 36345_Symbiodinium_microadriaticum_A1_CassKB8_500ln.fasta
```

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 36345_Symbiodinium_microadriaticum_A1_CassKB8_500ln.fasta
```

```
The total number of sequences is 36345
The average sequence length is 2094
The total number of bases is 76136570
The minimum sequence length is 500
The maximum sequence length is 78486
The N50 is 2811
Median Length = 1572
contigs < 150bp = 0
contigs >= 500bp = 36345
contigs >= 1000bp = 25680
contigs >= 2000bp = 12612
```

```
cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/36345_Symbiodinium_microadriaticum_A1_CassKB8_500ln.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/
```

### Siderastrea (my de novo transcriptome assembly)
Hybrid reference
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/Sid_hybridref_final.fasta

- includes my de novo host assembly
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/19222_Sid_GoodCoral_500lnThresh_Final.fasta

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 19222_Sid_GoodCoral_500lnThresh_Final.fasta
```

```
The total number of sequences is 19222
The average sequence length is 1127
The total number of bases is 21667620
The minimum sequence length is 500
The maximum sequence length is 32630
The N50 is 1211
Median Length = 827
contigs < 150bp = 0
contigs >= 500bp = 19222
contigs >= 1000bp = 6693
contigs >= 2000bp = 1829
```

### Breviolum symbiont assembly
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 52508_Breviolum_Sid-radians_Avila-Medina_500lnThresh_LongestIso_renamed.fasta
```

```
The total number of sequences is 52508
The average sequence length is 1178
The total number of bases is 61857756
The minimum sequence length is 500
The maximum sequence length is 2459
The N50 is 1372
Median Length = 687
contigs < 150bp = 0
contigs >= 500bp = 52508
contigs >= 1000bp = 28949
contigs >= 2000bp = 4806
```

### Cladocopium symbiont assembly
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/hybridref_final/47262_Davies_Cladocopium_LongestContig-500_renamed.fasta

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 47262_Davies_Cladocopium_LongestContig-500_renamed.fasta
```

```
The total number of sequences is 47262
The average sequence length is 1364
The total number of bases is 64485324
The minimum sequence length is 500
The maximum sequence length is 18168
The N50 is 1588
Median Length = 1519
contigs < 150bp = 0
contigs >= 500bp = 47262
contigs >= 1000bp = 27124
contigs >= 2000bp = 7623
```


----------------------------------------------------------------------------------------------
## Ortholog identification

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/

```
enable_lmod
module load container_env orthofinder/2.5.2
```

When running any command from orthofinder, add "crun" in front of it, for example:
    crun orthofinder

```
crun orthofinder -h
```


**SIMPLE USAGE:**
Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  orthofinder [options] -f <dir>

Add new species in <dir1> to previous run in <dir2> and run new analysis
  orthofinder [options] -f <dir1> -b <dir2>

OPTIONS:
```
 -t <int>        Number of parallel sequence search threads [Default = 20]
 -a <int>        Number of parallel analysis threads
 -d              Input is DNA sequences
 -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'
                 [Default = dendroblast] #Multiple sequence alignment
 -S <txt>        Sequence search program [Default = diamond]
                 Options: blast, diamond, diamond_ultra_sens, blast_gz, mmseqs, blast_nucl
 -A <txt>        MSA program, requires '-M msa' [Default = mafft]
                 Options: mafft, muscle
 -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]
                 Options: fasttree, raxml, raxml-ng, iqtree
 -s <file>       User-specified rooted species tree
 -I <int>        MCL inflation parameter [Default = 1.5]
 -x <file>       Info for outputting results in OrthoXML format
 -p <dir>        Write the temporary pickle files to <dir>
 -1              Only perform one-way sequence search
 -X              Don't add species names to sequence IDs
 -y              Split paralogous clades below root of a HOG into separate HOGs
 -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)
 -n <txt>        Name to append to the results directory
 -o <txt>        Non-default results directory
 -h              Print this help text

WORKFLOW STOPPING OPTIONS:
 -op             Stop after preparing input files for BLAST
 -og             Stop after inferring orthogroups
 -os             Stop after writing sequence files for orthogroups
                 (requires '-M msa')
 -oa             Stop after inferring alignments for orthogroups
                 (requires '-M msa')
 -ot             Stop after inferring gene trees for orthogroups 

WORKFLOW RESTART COMMANDS:
 -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>
 -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>
 -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>
```


**BLAST all-v-all algorithm?**


From OrthoFinder best practice webpage:

*One problem that can arise with transcriptomes is when you start with ~100,000 transcripts per species. This can be computationally expensive and will likely result in a very large number of files being generated, so be careful in such cases.*


I decided to run OrthoFinder separately for hosts versus symbionts.


----------------------------------------------------------------------------------------------
## Host ortholog identification - part A

Comparing Sid and Past longest Contig assembly

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/

```
nano orthofinder_host.sh
```

```
#!/bin/bash -l

#SBATCH -o orthofinder_host.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=orthofinder_host

enable_lmod
module load container_env orthofinder/2.5.2

crun orthofinder -f /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig -d -t 24
```

```
sbatch orthofinder_host.sh
```

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/
```
cat orthofinder_host.txt
```
Ran blast_nucl all-versus-all

```
Writing orthogroups to file
---------------------------
OrthoFinder assigned 13060 genes (26.8% of total) to 4319 orthogroups. Fifty percent of all genes were in orthogroups with 1 or more genes (G50 was 1) and were contained in the largest 15581 orthogroups (O50 was 15581). There were 2026 orthogroups with all species present and 1472 of these consisted entirely of single-copy genes.
```

*OrthoFinder note:  In general it’s nice to see at least 80% of your genes assigned to orthogroups. Fewer than this means that you are probably missing orthology relationships that actually exist for some of the remaining genes, poor species sampling is the most likely cause for this.*
...I think this is for organisms that are better sequenced/ annotated than corals...


**Results:**
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/OrthoFinder/

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/OrthoFinder/Results_Apr02/Comparative_Genomics_Statistics/
There's a tab-delimited file called Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv. Like other “.tsv” files from OrthoFinder this is best viewed in a spreadsheet program like Excel.

#### Orthologues directory

One of the most common reasons for running OrthoFinder is to find the orthologue of a gene you’re interested in.

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/OrthoFinder/Results_Apr02/Orthologues/Orthologues_29422_past_LongestContig_suffixed/

File:  29422_past_LongestContig_suffixed__v__19222_Sid_GoodCoral_500lnThresh_Final.tsv


#### Orthogroups

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/OrthoFinder/Results_Apr02/Orthogroups/

Often we’re interested in group-wise species comparisons, that is comparisons across a clade of species rather than between a pair of species. The generalisation of orthology to multiple species is the orthogroup. Just like orthologues are the genes descended from a single gene in the last common ancestor of a pair of species **an orthogroup is the set of genes descended from a single gene in a group of species.** 

Each gene tree from OrthoFinder, for example the one above, is for one orthogroup. The orthogroup gene tree is the tree we need to look at if we want it to include all pairwise orthologues. And even though some of the genes within an orthogroup can be paralogs of one another, if we tried to take any genes out then we would also be removing orthologs too.

So if we want to do a comparison of the ‘equivalent’ genes in a set of species, we need to do the comparison across the genes in an othogroup. The orthogroups are in the file Orthogroups/Orthogroups.tsv. This table has one orthogroup per line and one spcies per column and is ordered from largest orthogroup to smallest.

#### Orthogroup Sequences

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/OrthoFinder/Results_Apr02/Orthogroup_Sequences/

For each orthogroup there is a FASTA file in Orthogroup_Sequences/ which contains the sequences for the genes in that orthogroup.




----------------------------------------------------------------------------------------------
## Host ortholog identification - part B

Comparing Sid and Past longest Contig assembly >500 bp length assemblies

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/

```
nano orthofinder_host_b.sh
```

```
#!/bin/bash -l

#SBATCH -o orthofinder_host_b.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=orthofinder_host_b

enable_lmod
module load container_env orthofinder/2.5.2

crun orthofinder -f /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/ -d -t 24
```

```
sbatch orthofinder_host_b.sh
```

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/
```
cat orthofinder_host_b.txt
```
Ran blast_nucl all-versus-all

```
Writing orthogroups to file
---------------------------
OrthoFinder assigned 9500 genes (31.7% of total) to 3005 orthogroups. Fifty percent of all genes were in orthogroups with 1 or more genes (G50 was 1) and were contained in the largest 8473 orthogroups (O50 was 8473). There were 1267 orthogroups with all species present and 926 of these consisted entirely of single-copy genes.
```

*OrthoFinder note:  In general it’s nice to see at least 80% of your genes assigned to orthogroups. Fewer than this means that you are probably missing orthology relationships that actually exist for some of the remaining genes, poor species sampling is the most likely cause for this.*
...I think this is for organisms that are better sequenced/ annotated than corals...


**Results:**
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/OrthoFinder/

There's a tab-delimited file called Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv. Like other “.tsv” files from OrthoFinder this is best viewed in a spreadsheet program like Excel.

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/OrthoFinder/Results_Apr02/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv ./

#### Orthologues directory

One of the most common reasons for running OrthoFinder is to find the orthologue of a gene you’re interested in.

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/OrthoFinder/Results_Apr02/Orthologues/Orthologues_10714_past_LongestContig_500ln_suffixed/

File:  10714_past_LongestContig_500ln_suffixed__v__19222_Sid_GoodCoral_500lnThresh_Final.tsv


#### Orthogroups

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/OrthoFinder/Results_Apr02/Orthogroups/

Often we’re interested in group-wise species comparisons, that is comparisons across a clade of species rather than between a pair of species. The generalisation of orthology to multiple species is the orthogroup. Just like orthologues are the genes descended from a single gene in the last common ancestor of a pair of species **an orthogroup is the set of genes descended from a single gene in a group of species.** 

Each gene tree from OrthoFinder, for example the one above, is for one orthogroup. The orthogroup gene tree is the tree we need to look at if we want it to include all pairwise orthologues. And even though some of the genes within an orthogroup can be paralogs of one another, if we tried to take any genes out then we would also be removing orthologs too.

So if we want to do a comparison of the ‘equivalent’ genes in a set of species, we need to do the comparison across the genes in an othogroup. The orthogroups are in the file Orthogroups/Orthogroups.tsv. This table has one orthogroup per line and one spcies per column and is ordered from largest orthogroup to smallest.

#### Orthogroup Sequences

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/OrthoFinder/Results_Apr02/Orthogroup_Sequences/

For each orthogroup there is a FASTA file in Orthogroup_Sequences/ which contains the sequences for the genes in that orthogroup.



----------------------------------------------------------------------------------------------
## Symbiont ortholog identification

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/

```
nano orthofinder_sym.sh
```

```
#!/bin/bash -l

#SBATCH -o orthofinder_sym.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=orthofinder_sym

enable_lmod
module load container_env orthofinder/2.5.2

crun orthofinder -f /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/ -d -t 24
```

```
sbatch orthofinder_sym.sh
```

**Results:**
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/OrthoFinder/


> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/

cat orthofinder_sym.txt
```
Writing orthogroups to file
---------------------------
OrthoFinder assigned 53080 genes (39.0% of total) to 21527 orthogroups. Fifty percent of all genes were in orthogroups with 1 or more genes (G50 was 1) and were contained in the largest 36505 orthogroups (O50 was 36505). There were 960 orthogroups with all species present and 458 of these consisted entirely of single-copy genes.
```
Running blast_nucl all-versus-all

#### Orthogroups

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/OrthoFinder/Results_Mar25/Orthogroups/

#### Orthologues directory

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/OrthoFinder/Results_Mar25/Orthologues/

#### Orthogroup sequences

For each orthogroup there is a FASTA file in Orthogroup_Sequences/ which contains the sequences for the genes in that orthogroup.
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Symbiont/OrthoFinder/Results_Mar25/Orthogroup_Sequences/

----------------------------------------------------------------------------------------------

