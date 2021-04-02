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


### Method:  OrthoFinder
OrthoFinder is a software program for phylogenetic orthology inference.

[OrthoFinder GitHub](https://github.com/davidemms/OrthoFinder)
[OrthoFinder Tutorials](https://davidemms.github.io/)

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

**Results:**
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig/OrthoFinder/


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

**Results:**
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/orthofinder_Sid-Past/orthofind_Host/Sid.v.Past_longestContig_500ln/OrthoFinder/


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


----------------------------------------------------------------------------------------------

