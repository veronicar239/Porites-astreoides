---
layout: post
title: Porites astreoides transcriptome data processing
date: '2020-05-27'
categories: Protocols
tags: [RNASeq, Bioinformatics]
---

# Coral transcriptome data processing

### 2-year transplant of *Porites astreoides* 
- Puerto Morelos, Mexico
- from low pH ojo (submarine discharge spring) and high pH control site to low pH ojo and control site
- Ana Martinez (PhD candidate) & Adina Paytan (PI) & Dan Barshis (PI)

### Transplant design (information from Ana)
- Corals were taken from 3 different origins: within the ojo (center), outside the ojo (control) and in the reef.
- Corals were then fixed in 3 transplant sites: within the ojo (Ojo LAJA Center and Ojo NORTE Center) and in a control site (Ojo Laja CONTROL).
- From each colony, 3 replicates (or individual cores) were taken

### Analyzed by Veronica Radice, Barshis Lab, Old Dominion University

## Data
- raw data (year 2 of transplant, fastq g-zipped) backed up on RC drive
- contains data for multiple species:  *Porites astreoides*, *Porites porites*, *Siderastrea*

> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2017_April/gslserver.qb3.berkeley.edu/170419_50SR_HS4K2A/Paytan

- other raw data from Paytan (year 1 data?) backed up on RC drive:

> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/2014_October/AdinaOA_2014_October_rawdata.tar.gz

### Files:

> /RC/group/rc_barshis_lab/taxonarchive/AdinaOA/Porites_astreoides/
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan

----------------------------------------------------------------------------------------------
## Rename sequencing files with shorter names with relevant sample ID information 
renamer_advbioinf.py

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan

renamingtable_yr2.txt

Example:  

OldName	| NewName |
--- | --- |
APDB20_A_NO_261_C_Pa_yr2_A1_S1_L001_R1_001.fastq | NO_261_C_Pa_yr2_R1.fastq
APDB20_B_NO_264_C_Pa_yr2_B1_S2_L001_R1_001.fastq | NO_264_C_Pa_yr2_R1.fastq

``` sh
nano renamer_yr2.sh
```

``` sh
#!/bin/bash

#$ -cwd
#$ -o outfile_renamer_yr2.txt -j y
#$ -S /bin/bash
#$ -q main
#$ -m beas
#$ -M a1fernan@odu.edu

/cm/shared/courses/dbarshis/barshislab/Ana/170419_50SR_HS4K2A/Paytan/renamer_advbioinf.py ./renamingtable_yr2.txt
```

----------------------------------------------------------------------------------------------
## Make adapter list 
Be aware of formatting (e.g., extra spaces, extra lines) (in Excel)
Ideally, list is made as .txt file

First column with individual .fastq file name (starting with '>'), second column with associated adapter sequence, no header line in file

If adapter sequence is short, then need to find the long version of each adapter sequence
- need long version of adapter sequence to do full adapter trimming
- use Illumina adapter sequence list (for specific sequencer) to find “index” adapter sequences
- for this project, TruSeq (Index 1-27) was used 

Secure copy adapter list file from local machine to folder in the cluster.

----------------------------------------------------------------------------------------------
## Adapter trimming-clipping-quality filtering for all files 
can do this in separate batches

> /cm/shared/courses/dbarshis/15AdvBioinf/scripts/
Trimclipfilterstatsbatch_advbioinf.py

```
nano TrimClipFilter.sh
```

```
#!/bin/bash -l

#SBATCH -o TrimClipFilter.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=TrimClipFilter

/cm/shared/courses/dbarshis/15AdvBioinf/scripts/Trimclipfilterstatsbatch_advbioinf.py adapterlist.txt *.fastq
```

```
sbatch TrimClipFilter.sh
```

Check output file for information
```
cat TrimClipFilter.txt
```
----------------------------------------------------------------------------------------------
## QA-QC
- run trim clipped filtered stats Python script on stats output file 
- check adapter clipping stats for each file
cat trimclipstats.txt

in P_ast folder (both NC and NO files - both Porites astreoides)
```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/filteringstats/trimclipstats.txt /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/filteringstats/trimclipstats.txt > trimclipstats_final.txt
```

```
nano filterstats.sh
```

```
#!/bin/bash -l

#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=filter-stats

/cm/shared/courses/dbarshis/15AdvBioinf/scripts/Schafran_trimstatstable_advbioinf.py trimclipstats_final.txt P_ast_trimclipstats_out.txt
```

```
sbatch filterstats.sh
```

copy output of trim clip stats onto local machine

--------------------------------------------------------------------------------------------
## Test Host reference transcriptomes for *Porites astreoides*

#### Download reefgenomics transcriptome
- download CDS (coding sequence file) for P. astreoides
[http://comparative.reefgenomics.org/datasets.html](http://comparative.reefgenomics.org/datasets.html)

#### Download Carly Kenkel's transcriptome
- from Assistant Professor Carly Kenkel, University of Southern California Dornsife [https://dornsife.usc.edu/labs/carlslab](https://dornsife.usc.edu/labs/carlslab)
- CDS file is the translated protein coding sequences based on blastx against the uniprot swissprot database 
- the script for the translation is Misha’s
- pastreoides_2014 folder > pastreoides_may2014 (Dropbox) 
- File:  past_CDS.fas
- Citation:  Gene expression under chronic heat stress in populations of the mustard hill coral (Porites astreoides) from different thermal environments (2013). *Molecular Ecology* [https://doi.org/10.1111/mec.12390](https://doi.org/10.1111/mec.12390)
- NB:  in my shared directory, file known as P_ast_Kenkel
- only about 70% complete, and most of the genes are only partials based on comparisons to the BUSCO gene set
- Carly Kenkel has actually stopped using this assembly herself as its not as complete as she would like 
[https://matzlab.weebly.com/data--code.html](https://matzlab.weebly.com/data--code.html) 


#### Download updated Kenkel transcriptome
- Carly Kenkel and her students currently mapping to a version of the Mansour et al 2016 multi-life stage assembly
- [https://academic.oup.com/gigascience/article/5/1/s13742-016-0138-1/2720993](https://academic.oup.com/gigascience/article/5/1/s13742-016-0138-1/2720993) 
- Carly recleaned, refiltered and reannotated using a slightly modified version of Sheila Kitchen’s pipeline 
- pipeline that Carly has used for other de novo assemblies 
- [https://academic.oup.com/gigascience/article/6/9/gix074/4091589](https://academic.oup.com/gigascience/article/6/9/gix074/4091589)  
- This is 86-87% complete (again, based on BUSCO comparison) and 59% annotated
- This host reference is derived from the multi-life stage transcriptome of Porites astreoides
- originally sequenced by Mansour et al (2016) doi: [10.1186/s13742-016-0138-1](10.1186/s13742-016-0138-1) (Project accession: GEHP0000000)
- The holobiont Trinity assembly was downloaded by C. Kenkel and refiltered to identify the host-specific contigs as described in Kenkel & Bay (2017) 
- [https://doi.org/10.1093/gigascience/gix074](https://doi.org/10.1093/gigascience/gix074)

##### P_ast_full_Kenkel
- Past_host_All_iso.fasta 
> /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel/

##### Notes: 
> Contigs were filtered for biological contamination using a series of blastx searches to the Acropora digitifera 
> and Symbiodinium kawagutii proteomes and NCBI's nr database. Host contigs were clustered into isogroups (=genes) and 
> annotated using blastx comparisons against the UniProt Swiss-Prot database. Correspondence tables (.tab)list association 
> between contigs and isogroups (seq2iso), gene name (iso2gene), as well as Gene Ontology (iso2go), and KEGG (iso2kegg) annotations.
> Annotations also include a file	of protein translations	and corresponding coding sequences extracted from the transcriptome 
> based on blastx hits to the UniProtKB Swiss-Prot database. Annotations were completed in 03/2019.
>
> See <TBD> for additional details regarding the assembly.
> These data are free to use without restriction,  but please cite above refs if using this resource. 
> 
> Contact ckenkel[at]usc[dot]edu with questions/concerns.

##### from Carly Kenkel (04/16/2020):
- BUSCO stats for this version of the Mansour assembly are below, for your reference:

```
gVolante Analysis: ver.1.2.1

Summary of the Submitted Job:
Job ID:    201904030119-YL6PQ9V898XJXWB8
Project name:    PastLarvae
Fasta file:    Past_host_All_iso.fasta.gz
Cut-off length for sequence statistics and composition:    1
Sequence type:    trans
Selected program:    BUSCO_v2/v3
Selected reference gene set:    Metazoa

Completeness Assessment Results:    
Total # of core genes queried:    978
# of core genes detected
Complete:    843 (86.20%)
Complete + Partial:    850 (86.91%)
# of missing core genes:    128 (13.09%)
Average # of orthologs per core genes:    1.60
% of detected core genes that have more than 1 ortholog:    43.77
Scores in BUSCO format:    C:86.2%[S:48.5%,D:37.7%],F:0.7%,M:13.1%,n:978

for details of these metrics, see
 BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.
 Simao FA, Waterhouse RM, Ioannidis P, Kriventseva EV, Zdobnov EM.
 Bioinformatics. 2015 Oct 1;31(19):3210-2.

Length Statistics and Composition:
# of sequences:    82805
Total length (nt):    152548553
Longest sequence (nt):    28297
Shortest sequence (nt):    400
Mean sequence length (nt):    1842
Median sequence length (nt):    1341
N50 sequence length (nt):    2650
L50 sequence count:    17987
# of sequences >   1K (nt):    50969 (61.6% of total number)
# of sequences >  10K (nt):    259 (0.3% of total number)
# of sequences > 100K (nt):    0 (0.0% of total number)
# of sequences >   1M (nt):    0 (0.0% of total number)
# of sequences >  10M (nt):    0 (0.0% of total number)
Base composition (%):    A:29.43, T:29.30, G:20.68, C:20.59, N:0.00, Other:0.00
GC-content (%):    41.27
# of sequences containing non-ACGTN (nt):    0

gVolante Citation:
gVolante for standardizing completeness assessment of genome and transcriptome assemblies
Nishimura O, Hara Y, Kuraku S.
Bioinformatics. 2017 Nov 15;33(22):3635-3637.
```


### make reference assembly folder
make separate folders for different transriptomes in refassembly folder
```
mkdir refassembly
mkdir P_ast_Kenkel
mkdir P_ast_reefgenomics
```

secure copy transcriptomes to directories

--------------------------------------------------------------------------------------------
## Check quality of assemblies
avg_cov_len_fasta_advbioinf.py

#### P_ast_Kenkel
- past_CDS.fas

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fas
```

```
The total number of sequences is 18763
The average sequence length is 377
The total number of bases is 7077774
The minimum sequence length is 51
The maximum sequence length is 3372
The N50 is 435
Median Length = 360
contigs < 150bp = 2226
contigs >= 500bp = 3820
contigs >= 1000bp = 658
contigs >= 2000bp = 36
```

#### P_ast_full_Kenkel
- Past_host_All_iso.fasta
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fasta
```

```
The total number of sequences is 82805
The average sequence length is 1842
The total number of bases is 152548553
The minimum sequence length is 400
The maximum sequence length is 28297
The N50 is 2650
Median Length = 1638
contigs < 150bp = 0
contigs >= 500bp = 74032
contigs >= 1000bp = 51002
contigs >= 2000bp = 27346
```

#### P_ast_reefgenomics
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fna
```

```
The total number of sequences is 15755
The average sequence length is 519
The total number of bases is 8186466
The minimum sequence length is 300
The maximum sequence length is 7197
The N50 is 519
Median Length = 582
contigs < 150bp = 0
contigs >= 500bp = 5583
contigs >= 1000bp = 896
contigs >= 2000bp = 48
```

--------------------------------------------------------------------------------------------
## Make transcriptomes mappable
- need bowtie build module
- creates 6 files for mapping

```
module load bowtie2/2.2.4
```

bowtie for reefgenomics mapping
```
bowtie2-build *.fna past_rg
```

bowtie for Kenkel mapping
```
bowtie2-build *.fas past_k
```

- test alignments against two different assemblies using few test samples
- test filtered versus unfiltered (clipped and trimmed) fastq files
- test Kenkel transcriptome versus reefgenomics

#### Test reefgenomics P. astreoides transcriptome
mapping to reefgenomics (rg) reference
```
nano mapreads_past_rg.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_past_rg.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_past_rg

module load bowtie2/2.2.4

bowtie2 --local -x past_rg -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/NC_291_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_291_C_Pa_nof_past_rg.sam -k 5\n
bowtie2 --local -x past_rg -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/NC_327_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_327_N_Pa_nof_past_rg.sam -k 5\n
bowtie2 --local -x past_rg -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/NO_261_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_261_C_Pa_nof_past_rg.sam -k 5\n

bowtie2 --local -x past_rg -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/NC_291_C_Pa_yr2_R1_clippedtrimmedfilterd.fastq -S NC_291_C_Pa_f_past_rg.sam -k 5\n
bowtie2 --local -x past_rg -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/NC_327_N_Pa_yr2_R1_clippedtrimmedfilterd.fastq -S NC_327_N_Pa_f_past_rg.sam -k 5\n
bowtie2 --local -x past_rg -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/NO_261_C_Pa_yr2_R1_clippedtrimmedfilterd.fastq -S NO_261_C_Pa_f_past_rg.sam -k 5\n
```

```
sbatch mapreads_past_rg.sh
```

```
cat bowtie2_past_rg.txt 
```

##### NOT FILTERED
```
11663192 reads; of these:
  11663192 (100.00%) were unpaired; of these:
    10695849 (91.71%) aligned 0 times
    820882 (7.04%) aligned exactly 1 time
    146461 (1.26%) aligned >1 times
8.29% overall alignment rate

17526362 reads; of these:
  17526362 (100.00%) were unpaired; of these:
    16132922 (92.05%) aligned 0 times
    1142989 (6.52%) aligned exactly 1 time
    250451 (1.43%) aligned >1 times
7.95% overall alignment rate

13101819 reads; of these:
  13101819 (100.00%) were unpaired; of these:
    11979623 (91.43%) aligned 0 times
    950394 (7.25%) aligned exactly 1 time
    171802 (1.31%) aligned >1 times
8.57% overall alignment rate
```

##### FILTERED
```
11451437 reads; of these:
  11451437 (100.00%) were unpaired; of these:
    10498837 (91.68%) aligned 0 times
    808049 (7.06%) aligned exactly 1 time
    144551 (1.26%) aligned >1 times
8.32% overall alignment rate

17237818 reads; of these:
  17237818 (100.00%) were unpaired; of these:
    15864746 (92.03%) aligned 0 times
    1126109 (6.53%) aligned exactly 1 time
    246963 (1.43%) aligned >1 times
7.97% overall alignment rate

12893340 reads; of these:
  12893340 (100.00%) were unpaired; of these:
    11786973 (91.42%) aligned 0 times
    936760 (7.27%) aligned exactly 1 time
    169607 (1.32%) aligned >1 times
8.58% overall alignment rate
```

#### Test Kenkel P. astreoides transcriptome
- initially tried the partial past_CDS file
- mapping to Kenkel (k) reference

```
nano mapreads_past_k.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_past_k.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_past_k

module load bowtie2/2.2.4

bowtie2 --local -x past_k -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/NC_291_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_291_C_Pa_nof_past_k.sam -k 5\n
bowtie2 --local -x past_k -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/NC_327_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_327_N_Pa_nof_past_k.sam -k 5\n
bowtie2 --local -x past_k -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/NO_261_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_261_C_Pa_nof_past_k.sam -k 5\n

bowtie2 --local -x past_k -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/NC_291_C_Pa_yr2_R1_clippedtrimmedfilterd.fastq -S NC_291_C_Pa_f_past_k.sam -k 5\n
bowtie2 --local -x past_k -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/NC_327_N_Pa_yr2_R1_clippedtrimmedfilterd.fastq -S NC_327_N_Pa_f_past_k.sam -k 5\n
bowtie2 --local -x past_k -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/NO_261_C_Pa_yr2_R1_clippedtrimmedfilterd.fastq -S NO_261_C_Pa_f_past_k.sam -k 5\n
```

```
sbatch mapreads_past_k.sh
```

```
cat bowtie2_past_k.txt 
```

##### NOT FILTERED
```
11663192 reads; of these:
  11663192 (100.00%) were unpaired; of these:
    9226577 (79.11%) aligned 0 times
    1245525 (10.68%) aligned exactly 1 time
    1191090 (10.21%) aligned >1 times
20.89% overall alignment rate

17526362 reads; of these:
  17526362 (100.00%) were unpaired; of these:
    13184983 (75.23%) aligned 0 times
    2104270 (12.01%) aligned exactly 1 time
    2237109 (12.76%) aligned >1 times
24.77% overall alignment rate

13101819 reads; of these:
  13101819 (100.00%) were unpaired; of these:
    10263420 (78.34%) aligned 0 times
    1515973 (11.57%) aligned exactly 1 time
    1322426 (10.09%) aligned >1 times
21.66% overall alignment rate
```

##### FILTERED
```
11451437 reads; of these:
  11451437 (100.00%) were unpaired; of these:
    9050624 (79.03%) aligned 0 times
    1227233 (10.72%) aligned exactly 1 time
    1173580 (10.25%) aligned >1 times
20.97% overall alignment rate

17237818 reads; of these:
  17237818 (100.00%) were unpaired; of these:
    12955421 (75.16%) aligned 0 times
    2074960 (12.04%) aligned exactly 1 time
    2207437 (12.81%) aligned >1 times
24.84% overall alignment rate

12893340 reads; of these:
  12893340 (100.00%) were unpaired; of these:
    10093502 (78.28%) aligned 0 times
    1495319 (11.60%) aligned exactly 1 time
    1304519 (10.12%) aligned >1 times
21.72% overall alignment rate
```

--------------------------------------------------------------------------------------------
## Test to compare different P. astreoides references
Alignment is very low for /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/Past_host_All_iso_suffixed.fasta

##### Test past_CDS_Kenkel_2014
- from pastreoides_may2014 folder
- File: past_CDS.fas (renamed KenkelTranscriptome.fasta)
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel

##### Test past_full_Kenkel_2014
- from pastreoides_may2014 folder
- File: past.fasta
/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014

##### Test Past_host_All_iso_2020
- most recent Kenkel reference transcriptome (personally shared)
- File: Past_host_All_iso.fasta
/cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/

#### Check quality of assemblies

Test past_CDS_Kenkel_2014
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py KenkelTranscriptome.fasta
```

```
The total number of sequences is 18763
The average sequence length is 377
The total number of bases is 7077774
The minimum sequence length is 51
The maximum sequence length is 3372
The N50 is 435
Median Length = 360
contigs < 150bp = 2226
contigs >= 500bp = 3820
contigs >= 1000bp = 658
contigs >= 2000bp = 36
```

Test past_full_Kenkel_2014
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py past.fasta
```

```
The total number of sequences is 30740
The average sequence length is 550
The total number of bases is 16907062
The minimum sequence length is 100
The maximum sequence length is 8171
The N50 is 661
Median Length = 653
contigs < 150bp = 899
contigs >= 500bp = 11640
contigs >= 1000bp = 3274
contigs >= 2000bp = 325
```

Test Past_host_All_iso_2020
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Past_host_All_iso.fasta
```

```
The total number of sequences is 82805
The average sequence length is 1842
The total number of bases is 152548553
The minimum sequence length is 400
The maximum sequence length is 28297
The N50 is 2650
Median Length = 1638
contigs < 150bp = 0
contigs >= 500bp = 74032
contigs >= 1000bp = 51002
contigs >= 2000bp = 27346
```

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/KenkelTranscriptome.fasta test_past_CDS_Kenkel_2014

bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/past.fasta test_past_full_Kenkel_2014

bowtie2-build /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/Past_host_All_iso.fasta test_Past_host_All_iso_2020
```

#### Compare 3 P. astreoides assemblies with subset of samples

##### Test past_CDS_Kenkel_2014
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/ 

```
nano mapreads_past.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_test_past_CDS_Kenkel_2014.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_test_past_CDS_Kenkel_2014

module load bowtie2/2.2.4

bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_294_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_294_C_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_282_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_282_C_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_319_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_319_C_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_318_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_318_N_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_272_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_272_N_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_289_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_289_La_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_320_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_320_La_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/test_past_CDS_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_292_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_292_La_Pa_nof_test_past_CDS_Kenkel_2014.sam -k 5\n
```

```
sbatch mapreads_past.sh
```

```
cat bowtie2_test_past_CDS_Kenkel_2014.txt 
```

##### Test past_full_Kenkel_2014
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014

```
nano mapreads_past.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_test_past_full_Kenkel_2014.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_test_past_full_Kenkel_2014

module load bowtie2/2.2.4

bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_294_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_294_C_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_282_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_282_C_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_319_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_319_C_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_318_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_318_N_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_272_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_272_N_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_289_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_289_La_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_320_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_320_La_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/test_past_full_Kenkel_2014 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_292_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_292_La_Pa_nof_test_past_full_Kenkel_2014.sam -k 5\n
```

```
sbatch mapreads_past.sh
```

```
cat bowtie2_test_past_full_Kenkel_2014.txt 
```

##### Test Past_host_All_iso_2020
> /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/

```
nano mapreads_past.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_test_Past_host_All_iso_2020.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_test_Past_host_All_iso_2020

module load bowtie2/2.2.4

bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_294_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_294_C_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_282_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_282_C_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_319_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_319_C_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_318_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_318_N_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_272_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_272_N_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_289_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_289_La_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_320_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_320_La_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/test_Past_host_All_iso_2020 -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_292_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_292_La_Pa_nof_test_Past_host_All_iso_2020.sam -k 5\n
```

```
sbatch mapreads_past.sh
```

```
cat bowtie2_test_Past_host_All_iso_2020.txt
```

##### Test ReefGenomics P. astreoides reference - secondary test
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/

```
module load bowtie2/2.2.4
bowtie2-build Porites_astreoides_cds_100.final.clstr.fna past_ReefGenomics
```

```
nano mapreads_past_rg_second-test.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_past_rg_second-test.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_past_rg_second-test

module load bowtie2/2.2.4

bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_294_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_294_C_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_282_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_282_C_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_319_C_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_319_C_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_318_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_318_N_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_272_N_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_272_N_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NO_289_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NO_289_La_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_320_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_320_La_Pa_nof_past_ReefGenomics.sam -k 5\n
bowtie2 --local -x /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_reefgenomics/past_ReefGenomics -U /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/NC_292_La_Pa_yr2_R1_clippedtrimmed_nofilter.fastq -S NC_292_La_Pa_nof_past_ReefGenomics.sam -k 5\n
```

```
sbatch mapreads_past_rg_second-test.sh
```

```
cat bowtie2_past_rg_second-test.txt 
```


Samples had greatest alignment with 2014 full P. astreoides assembly (file:  past.fasta; Kenkel et al. 2013), so using this for analysis.


--------------------------------------------------------------------------------------------
## Identify symbiont type of P. astreoides
- based on internal transcribed spacer 2 (ITS2) region, a multicopy genetic marker commonly used to analyse Symbiodinium diversity
- test ITS2 fasta files 
- Arif et al. 2014 [https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869)
- File:  ArifITS2_mec12869-sup-0001-FileS1.txt

make separate folders for different references in refassembly folder
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/Arif_ITS2

#### Check quality of assemblies
Arif_ITS2
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.txt
```

```
The total number of sequences is 433
The average sequence length is 279
The total number of bases is 120855
The minimum sequence length is 152
The maximum sequence length is 376
The N50 is 283
Median Length = 280
contigs < 150bp = 0
contigs >= 500bp = 0
contigs >= 1000bp = 0
contigs >= 2000bp = 0
```

--------------------------------------------------------------------------------------------
## Map all P. astreoides sequence files to Arif ITS2 (all clades) reference
- ran separate jobs from respective directories (subset-fastq_NC/QCFastqs/nofilter and subset-fastq_NO/QCFastqs/nofilter)
- set up parallel array alignment jobs using (Misha Matz's) bowtie parameters 
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter

##### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for Arif ITS2 clade mapping
```
bowtie2-build ArifITS2_mec12869-sup-0001-FileS1.txt Arif_ITS2
```

##### Mapping to Arif_ITS2 reference

```
nano mapreads_Arif_ITS2.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Arif_ITS2_NO.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Arif_ITS2_NO

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Arif_ITS2 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Arif_ITS2.sam -k 5\n; done
```

```
sbatch mapreads_Arif_ITS2.sh
```

now do the same for directory /subset-fastq_NO/QCFastqs/nofilter

--------------------------------------------------------------------------------------------
## Count expression - all reads mapped to Arif_ITS2 symbiont reference

```
nano countexpression_Arif_ITS2.sh 
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Arif_ITS2_NO.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_Arif_ITS2_NO

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/*.sam
```

```
sbatch countexpression_Arif_ITS2.sh
```

then do the same for subset-fastq_NC

--------------------------------------------------------------------------------------------
## Merge all *nof_Arif_ITS2_counts.txt files (NC and NO) into one big table 
- first need to add column to each .txt file with unique sample id
- so we can later identify which sample has which ITS2 sequences

Using regular expressions:

for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done

- For each file in the list, this will use sed to append to the end of each line a tab and the filename
- Using the -i flag with sed to perform a replacement in-place, overwriting the file
- Perform a substitution with s/PATTERN/REPLACEMENT/. 
- In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), 
- and $f is the filename, from the loop variable. 
- The s/// command is within double-quotes so that the shell can expand variables
- sed is most practical for pattern substitution and saving in-place. 

```
#!/bin/bash -l
for f in *_R1_nof_Arif_ITS2_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

--------------------------------------------------------------------------------------------
## Concatenate all NC and NO Arif_ITS2_counts.txt files into one merged file
NC folder
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/

```
nano append_filename.sh
```

```
#!/bin/bash -l

#SBATCH -o append_filename.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=append_filename

for f in *_R1_nof_Arif_ITS2_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

```
sbatch append_filename.sh
```

Then concatenate all NC files
```
cat *_R1_nof_Arif_ITS2_counts.txt > merged_NC_ITS2_counts.txt
```

NO folder
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/

```
nano append_filename.sh
```

```
sbatch append_filename.sh
```

```
cat *_R1_nof_Arif_ITS2_counts.txt > merged_NO_ITS2_counts.txt
```

```
cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/merged_NC_ITS2_counts.txt /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/
```

merge together
```
cat merged_NC_ITS2_counts.txt merged_NO_ITS2_counts.txt > /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/merged_Past_ITS2_counts.txt
```

#### Outcome of symbiont ITS2 clademapping
For n=40 Porites astreoides samples with matched reads:
- 87% reads mapped to I sequences
  - BLAST confirmed 100% match to former Symbiodinium sp. partial rRNA genes and ITS2
- After removing mapping to I sequences
  - 83% mapped to A
  - 9% mapped to C


--------------------------------------------------------------------------------------------
## Map 4 dominant ITS2 clade (Arif) contigs (contig with highest number of reads for one of the 4 dominant clades)
- took ITS2 contigs across clade A, F, C, I
- ran separate jobs from respective directories (subset-fastq_NC/QCFastqs/nofilter and subset-fastq_NO/QCFastqs/nofilter)
- set up parallel array alignment jobs using (Misha Matz's) bowtie parameters 
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter

##### Make file mappable
need bowtie build module - creates 6 files for mapping

```
module load bowtie2/2.2.4
```

bowtie mapping
```
bowtie2-build Past_mapped_Arif-ITS2_4-dominant-clades.txt top-4-clade_ITS2
```

#### Mapping top-4-clade_ITS2 contigs

```
nano mapreads_top-4-clade_ITS2.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_top-4-clade_ITS2_NC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_top-4-clade_ITS2_NC

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x top-4-clade_ITS2 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_top-4-clade_ITS2.sam -k 5\n; done
```

```
sbatch mapreads_top-4-clade_ITS2.sh
```

now do the same for directory /subset-fastq_NO/QCFastqs/nofilter

#### Count expression - all reads mapped to top-4-clade_ITS2

```
nano countexpression_top-4-clade_ITS2.sh 
```

```
#!/bin/bash -l

#SBATCH -o countexpression_top-4-clade_ITS2_NC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_top-4-clade_ITS2_NC

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/*_nof_top-4-clade_ITS2.sam
```

```
sbatch countexpression_top-4-clade_ITS2.sh
```

then do the same for subset-fastq_NC

--------------------------------------------------------------------------------------------
## Parse expression to big table 
- for 4 dominant ITS2 clade contigs
- copy all *_nof_top-4-clade_ITS2_counts.txt files (all files NO and NC) to one folder
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/merged_NO-NC_counts/

##### Make genelist.txt
format should be as follows:

GeneName |
--- |
contig1name |
contig2name |
contig3name |
contig4name |

secure copy from local machine

```
nano ParseExpression.sh
```

```
#!/bin/bash -l

#SBATCH -o ParseExpression.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=ParseExpression

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/ParseExpression2BigTable_advbioinf.py genelist.txt merged_top-4-ITS2-clades_counts.txt nomatch **_nof_top-4-clade_ITS2_counts.txt
```

```
sbatch ParseExpression.sh
```

#### Output 
merged_top-4-ITS2-clades_counts.txt

GeneName | Sum_Unique total reads of all samples (n=40) |
--- | --- |
GS_A4.1 | 30204
GS_F4.2a | 805
GS_I1 | 408932
LJ_C161 | 4557

Need to check clade I sequences:
- sequence unique to clade I or other ITS2 clades?
- check alignment, BLAST

**All samples matched to at least one clade contig.**

--------------------------------------------------------------------------------------------
## Samtools 
- view alignment
- copy all *_nof_top-4-clade_ITS2.sam files (all files NO and NC) to one folder
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/merged_NO-NC_sam/

[http://biobits.org/samtools_primer.html#Tutorial](http://biobits.org/samtools_primer.html#Tutorial)
[https://github.com/davetang/learning_bam_file](https://github.com/davetang/learning_bam_file)
[https://davetang.org/wiki/tiki-index.php?page=SAMTools](https://davetang.org/wiki/tiki-index.php?page=SAMTools)

- SAM files consist of two types of lines: headers and alignments. 
- Headers begin with @, and provide meta-data regarding the entire alignment file. 
- Alignments begin with any character except @, and describe a single alignment of a sequence read against the reference genome.
- Two factors are most important to note. 
- First, each read in a FASTQ file may align to multiple regions within a reference genome, and an individual read can therefore result in multiple alignments. 
- In the SAM format, each of these alignments is reported on a separate line. 
- Second, each alignment has 11 mandatory fields, followed by a variable number of optional fields.

```
module avail samtools/
module load samtools/
samtools   # see which version is running
```

#### Convert the SAM files to BAM
- This is an important prerequisite, as all the downstream steps, including the identification of genomic variants and visualization of reads, require BAM as input.
- To convert from SAM to BAM, use the "samtools view" command
```
samtools view -b -S -o alignments/sim_reads_aligned.bam alignments/sim_reads_aligned.sam
```

> -b: indicates that the output is BAM.
>
> -S: indicates that the input is SAM.
>
> -o: specifies the name of the output file

BAM files are stored in a compressed, binary format, and cannot be viewed directly 
- .bams are binary versions of .sams so are just more memory efficient 
- sorted .bams are sorted by contig, not read, so that all reads that match a particular contig are grouped together 
- which is what we want for an alignment/genotyping/sequence analysis 
- (in this case, we don't care for a gene expression counts analysis)

You can use the same view command to display all alignments. For example, running:
```
samtools view alignments/sim_reads_aligned.bam | more
```
will display all your reads in the unix more paginated style.

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/merged_NO-NC_sam/

```
nano sam_view-sort-index.sh
```

```
#!/bin/bash -l

#SBATCH -o sam_view-sort-index.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=sam_view-sort-index

module load samtools/1.1

for i in *_nof_top-4-clade_ITS2.sam; do `samtools view -bS $i > ${i%.sam}_unsorted.bam`; done

for i in *unsorted.bam; do samtools sort -O bam -T ${i%_unsorted.bam} $i > ${i%_sorted.bam}.bam
samtools index ${i%_sorted.bam}.bam
done
```

```
sbatch sam_view-sort-index.sh
```

#### Dan's script - for individual sample
```
#!/bin/bash -l

#SBATCH -o sam_view-sort-index_djb.txt
#SBATCH -n 1
#BATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=sam_view-sort-index

module load samtools/1.1

samtools view -bS NO_289_La_Pa_yr2_R1_nof_top-4-clade_ITS2.sam
samtools sort -O bam -T djbNO_289_La_Pa_yr2_R1_nof_top-4-clade_ITS2 djbNO_289_La_Pa_yr2_R1_nof_top-4-clade_ITS2_unsorted.bam > djbNO_289_La_Pa_yr2_R1_nof_top-4-clade_ITS2_sorted.bam
samtools index djbNO_289_La_Pa_yr2_R1_nof_top-4-clade_ITS2_sorted.bam
```

#### Sort and index bam files
- There are two options for sorting BAM files: by read name (-n), and by genomic location (default). 
- As our goal is to call genomic variants, and this requires that we “pile-up” all matching reads within a specific genomic location, we sort by location

```
samtools sort
samtools index 
```

#### Visualizing Reads
For the final step, we will use the SAMtools tview command to view our simulated reads and visually compare them to the reference genome. 
```
samtools tview *unsorted.bam.bam Past_mapped_Arif-ITS2_4-dominant-clades.txt
```
In this view, the first line shows the genome coordinates, the second line shows the reference sequence, and the third line shows the consensus sequence determined from the aligned reads. 

Throughout tview, a '.' indicates a match to the reference genome

Alternate (combining multiple steps):
```
samtools view -b -h aligned_reads.sam > aligned_reads.bam
```

#### Use bioSyntax to prettify your output.
**in VS code**
```
samtools view *unsorted.bam.bam | sam-less
```

#### Outcome ITS2 information for deciding to map to which symbiont "clade"

##### CLUSTAL
- checked for reverse complement of sequence of clade C contig (C161)
- then see proper alignment with other clades

##### BLAST
- max target seq: 1000
- is portion of clade I sequence (end?) unique to I or matches other ITS2 clade seq?
- clade I is matching at very end and very beginning
- perfect match with clade A

## Conclusion: 
**Our samples have majority ITS2 clade A sequences**
- now we will test multiple clade A transcriptomes


--------------------------------------------------------------------------------------------
## Testing (former ITS2 clade A) transcriptomes
##### Mapping against S_tridacnidorum_CCMP2592_CDS
- test Raúl González-Pech et al. 2019 *Symbiodinium tridacnidorum*
- S.tridacnidorum_CCMP2592.CDS.fna 
- ITS2 type A3
[https://www.biorxiv.org/content/10.1101/783902v1](https://www.biorxiv.org/content/10.1101/783902v1)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_tridacnidorum_CCMP2592_CDS/

#### Check quality of assemblies
S_tridacnidorum_CCMP2592_CDS
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fna
```

```
The total number of sequences is 45474
The average sequence length is 2033
The total number of bases is 92471373
The minimum sequence length is 150
The maximum sequence length is 57981
The N50 is 2916
Median Length = 4590
contigs < 150bp = 0
contigs >= 500bp = 40321
contigs >= 1000bp = 30253
contigs >= 2000bp = 15298
```

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /home/vradice/P_ast/refassembly/S_tridacnidorum_CCMP2592_CDS/S.tridacnidorum_CCMP2592.CDS.fna S_tridacnidorum_CCMP2592_CDS
```

#### Mapping to S_tridacnidorum_CCMP2592_CDS reference

```
nano mapreads_S_tridacnidorum_CCMP2592_CDS.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_S_tridacnidorum_CCMP2592_CDS.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_S_tridacnidorum_CCMP2592_CDS

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x S_tridacnidorum_CCMP2592_CDS -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_S_tridacnidorum_CCMP2592.sam -k 5\n; done
```

```
sbatch mapreads_S_tridacnidorum_CCMP2592_CDS.sh
```

```
head bowtie2_S_tridacnidorum_CCMP2592_CDS.txt
cat bowtie2_S_tridacnidorum_CCMP2592_CDS.txt
```

## Count expression 

```
nano countexpression_S_tridacnidorum_CCMP2592.sh 
```

```
#!/bin/bash -l

#SBATCH -o countexpression_S_tridacnidorum_CCMP2592.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_S_tridacnidorum_CCMP2592

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /home/vradice/P_ast/clippedtrimmed_nofilter/*_nof_S_tridacnidorum_CCMP2592.sam
```

```
sbatch countexpression_S_tridacnidorum_CCMP2592.sh
```

--------------------------------------------------------------------------------------------
## Testing (former ITS2 clade A) transcriptomes
Mapping against *Symbiodinium microadriaticum*
- González-Pech et al. 2019
- S.microadriaticum_CassKB8.CDS.fna
[https://www.biorxiv.org/content/10.1101/800482v1.full](https://www.biorxiv.org/content/10.1101/800482v1.full)
> /home/vradice/P_ast/refassembly/S_microadriaticum_CassKB8/

#### Check quality of assemblies
### S_microadriaticum_CassKB8

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fna
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

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /home/vradice/P_ast/refassembly/S_microadriaticum_CassKB8/S.microadriaticum_CassKB8.CDS.fna S_microadriaticum_CassKB8
```

#### Mapping to S_microadriaticum_CassKB8 reference
```
nano mapreads_S_microadriaticum_CassKB8.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_S_microadriaticum_CassKB8.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_S_microadriaticum_CassKB8

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x S_microadriaticum_CassKB8 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_S_microadriaticum_CassKB8.sam -k 5\n; done
```

sbatch mapreads_S_microadriaticum_CassKB8.sh

#### Count expression 
```
nano countexpression_S_microadriaticum_CassKB8.sh 
```

```
#!/bin/bash -l

#SBATCH -o countexpression_S_microadriaticum_CassKB8.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_S_microadriaticum_CassKB8

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /home/vradice/P_ast/clippedtrimmed_nofilter/*_nof_S_microadriaticum_CassKB8.sam
```

```
sbatch countexpression_S_microadriaticum_CassKB8.sh
```

--------------------------------------------------------------------------------------------
## Testing (former ITS2 clade A) transcriptomes
- Mapping against *Symbiodinium microadriaticum* CCMP2467
- Genome originally from Aranda et al. 2016 (isolated from culture) [https://www.nature.com/articles/srep39734](https://www.nature.com/articles/srep39734)
- Chen et al. 2019 "Revised genome sequences and annotations of six Symbiodiniaceae taxa"
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
  - reannotated the genome - *Symbiodinium microadriaticum*
  - Symbiodinium_microadriaticum.CDS.fna
    [https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CCMP2467/

#### Check quality of assemblies
S_microadriaticum_CCMP2467
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fna
```

```
The total number of sequences is 29728
The average sequence length is 2217
The total number of bases is 65922655
The minimum sequence length is 130
The maximum sequence length is 51201
The N50 is 3294
Median Length = 1002
contigs < 150bp = 1
contigs >= 500bp = 26755
contigs >= 1000bp = 20390
contigs >= 2000bp = 11116
```

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CCMP2467/Symbiodinium_microadriaticum.CDS.fna S_microadriaticum_CCMP2467
```

#### Mapping to S_microadriaticum_CCMP2467 reference
```
nano mapreads_S_microadriaticum_CCMP2467.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_S_microadriaticum_CCMP2467.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_S_microadriaticum_CCMP2467

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x S_microadriaticum_CCMP2467 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_S_microadriaticum_CCMP2467.sam -k 5\n; done
```

```
sbatch mapreads_S_microadriaticum_CCMP2467.sh
```

#### Count expression 
```
nano countexpression_S_microadriaticum_CCMP2467.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_S_microadriaticum_CCMP2467.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_S_microadriaticum_CCMP2467

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter/*_nof_S_microadriaticum_CCMP2467.sam
```

```
sbatch countexpression_S_microadriaticum_CCMP2467.sh
```

--------------------------------------------------------------------------------------------
## Testing (former ITS2 clade A) transcriptomes
- mapping against *Symbiodinium tridacnidorum* 
- original Shoguchi et al. 2018 genome [https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4857-9](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4857-9)
- Chen et al. 2019
  - re-annotated this genome
  - Symbiodinium_tridacnidorum.CDS.fna
    [https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /home/vradice/P_ast/refassembly/S_tridacnidorum_Chen-Shoguchi

#### Check quality of assemblies
S_tridacnidorum_Chen-Shoguchi
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fna
```

```
The total number of sequences is 25808
The average sequence length is 1424
The total number of bases is 36771813
The minimum sequence length is 150
The maximum sequence length is 15846
The N50 is 1890
Median Length = 651
contigs < 150bp = 0
contigs >= 500bp = 20976
contigs >= 1000bp = 13739
contigs >= 2000bp = 4990
```

#### Make file mappable
need bowtie build module - creates 6 files for mapping
module load bowtie2/2.2.4

bowtie for mapping
```
bowtie2-build /home/vradice/P_ast/refassembly/S_tridacnidorum_Chen-Shoguchi/Symbiodinium_tridacnidorum.CDS.fna S_tridacnidorum_Chen-Shoguchi
```

#### Mapping to S_tridacnidorum_Chen-Shoguchi reference
```
nano mapreads_S_tridacnidorum_Chen-Shoguchi.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_S_tridacnidorum_Chen-Shoguchi.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_S_tridacnidorum_Chen-Shoguchi

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x S_tridacnidorum_Chen-Shoguchi -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_S_tridacnidorum_Chen-Shoguchi.sam -k 5\n; done
```

```
sbatch mapreads_S_tridacnidorum_Chen-Shoguchi.sh
```

#### Count expression 
```
nano countexpression_S_tridacnidorum_Chen-Shoguchi.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_S_tridacnidorum_Chen-Shoguchi.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_S_tridacnidorum_Chen-Shoguchi

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /home/vradice/P_ast/clippedtrimmed_nofilter/*_nof_S_tridacnidorum_Chen-Shoguchi.sam
```

```
sbatch countexpression_S_tridacnidorum_Chen-Shoguchi.sh
```

--------------------------------------------------------------------------------------------
## Testing (former ITS2 clade A) transcriptomes
- mapping against *Symbiodinium microadriaticum* 04-503SCI.03
- S.microadriaticum_04-503SCI.03.CDS.fna
- González-Pech et al. 2019
  [https://www.biorxiv.org/content/10.1101/800482v1.full](https://www.biorxiv.org/content/10.1101/800482v1.full)
> /home/vradice/P_ast/refassembly/S_microadriaticum_04-503SCI.03/

#### Check quality of assemblies
S_microadriaticum_04-503SCI.03
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py *.fna
```

```
The total number of sequences is 38462
The average sequence length is 1814
The total number of bases is 69777578
The minimum sequence length is 109
The maximum sequence length is 75977
The N50 is 2748
Median Length = 4623
contigs < 150bp = 9
contigs >= 500bp = 32203
contigs >= 1000bp = 22567
contigs >= 2000bp = 11237
```

#### Make file mappable
need bowtie build module - creates 6 files for mapping
module load bowtie2/2.2.4

bowtie for mapping
```
bowtie2-build /home/vradice/P_ast/refassembly/S_microadriaticum_04-503SCI.03/S.microadriaticum_04-503SCI.03.CDS.fna S_microadriaticum_04-503SCI.03
```

#### Mapping to S_microadriaticum_04-503SCI.03 reference
```
nano mapreads_S_microadriaticum_04-503SCI.03.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_S_microadriaticum_04-503SCI.03.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_S_microadriaticum_04-503SCI.03

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x S_microadriaticum_04-503SCI.03 -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_S_microadriaticum_04-503SCI.03.sam -k 5\n; done
```

```
sbatch mapreads_S_microadriaticum_04-503SCI.03.sh
```

#### Count expression 
```
nano countexpression_S_microadriaticum_04-503SCI.03.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_S_microadriaticum_04-503SCI.03.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_S_microadriaticum_04-503SCI.03

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /home/vradice/P_ast/clippedtrimmed_nofilter/*_nof_S_microadriaticum_04-503SCI.03.sam
```

```
sbatch countexpression_S_microadriaticum_04-503SCI.03.sh
```

--------------------------------------------------------------------------------------------
## Compare 4 clade A transcriptomes

#### Check number of reads aligning exactly 1 time:

Name | Isolated from | AVG_Proportion_Single_aligned | AVG_NumSingleAligned |
--- | --- | --- | --- |
S_microadriaticum_CassKB8 | Cassiopea (Hawai'i) | 0.0112 | 158675
S_microadriaticum_04-503SCI.03 | Orbicella faveolata (Florida) | 0.0107 | 152557
S_microadriaticum_CCMP2467 | Culture | 0.0078 | 108401
S_tridacnidorum_CCMP2592 | Heliofungia actiniformis (Great Barrier Reef) | 0.0048 | 66033
S_tridacnidorum_Chen-Shoguchi | Tridacna crocea (Okinawa, 1980; then cultured) | 0.0026 | 36603


**Decided to use S_microadriaticum_CassKB8 transcriptome for symbiont reference.**

--------------------------------------------------------------------------------------------
## Add suffix to fasta reference sequence names
addsuffixtofastaseqnames.py
> /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/

#### Add suffix Past to host reference
```
nano addsuffixtofastaseqnames.sh
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_Kenkel.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=addsuffixtofastaseqnames

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Past /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/Past_host_All_iso.fasta
```

```
sbatch addsuffixtofastaseqnames.sh
```

output:
> /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/Past_host_All_iso_suffixed.fasta

#### Add suffix Sym to Symbiodinium reference
> cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/

```
mv S.microadriaticum_CassKB8.CDS.fna SymTranscriptome.fasta
```

```
nano addsuffixtofastaseqnames.sh 
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_CassKB8.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=addsuffixtofastaseqnames

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/addsuffixtofastaseqnames.py Sym SymTranscriptome.fasta
```

```
sbatch addsuffixtofastaseqnames.sh
```

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py SymTranscriptome_suffixed.fasta
```
The total number of sequences is 42652. In this case, it's the number of genes.

Rename fasta file based on # of contigs
```
mv SymTranscriptome_suffixed.fasta 42652_SymTranscriptome_suffixed.fasta
```

##### Different file size of KenkelTranscriptome_suffixed versus original KenkelTranscriptome
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel

but the number of sequences did not change
```
ls -alh
```

Also difference in file size of updated Kenkel full transcriptome
```
-rwxrwxrwx 1 hpc-0225 users 153M Apr 16 11:51 Past_host_All_iso.fasta
-rwxrwxrwx 1 vradice  users 150M Apr 17 09:26 Past_host_All_iso_suffixed.fasta
```

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py KenkelTranscriptome.fasta
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py KenkelTranscriptome_suffixed.fasta
```

##### Check the first 100 names to compare the original vs suffixed versions
```
grep '>' | head -100 KenkelTranscriptome_suffixed.fasta
grep '>' | head -100 KenkelTranscriptome.fasta
```
original Kenkel transcriptome had long header text - all good

--------------------------------------------------------------------------------------------
## Hybrid reference transcriptome
Concatenate symbiont transcriptome with host transcriptome to make mapping file
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/hybridref/

```
cat /cm/shared/courses/dbarshis/barshislab/referenceomes/Porites_astreoides/Past_host_All_iso_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/SymTranscriptome_suffixed.fasta > hybridreference.fasta
```

```
ls -alh
```

```
-rwxrwxrwx  1 vradice users 227M Apr 17 09:30 hybridreference.fasta
```

--------------------------------------------------------------------------------------------
## Map all the P_ast samples to concatenated (host plus symbiont) hybridreference transcriptome
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

#### Make file mappable
need bowtie build module - creates 6 files for mapping
module load bowtie2/2.2.4

bowtie for mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/hybridref/hybridreference.fasta hybridreference
```

Porites astreoides: NC_ and NO_ samples
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter
> *_Pa_yr2_R1_clippedtrimmed_nofilter.fastq

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NO/QCFastqs/nofilter
> *_Pa_yr2_R1_clippedtrimmed_nofilter.fastq

For both NO and NC files:
```
nano copy-fastq.sh
```

```
#!/bin/bash -l

#SBATCH -o copy-fastq_NC.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=copy-fastq

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/subset-fastq_NC/QCFastqs/nofilter/*_Pa_yr2_R1_clippedtrimmed_nofilter.fastq /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/
```

```
sbatch copy-fastq.sh
```

#### Mapping to hybridreference
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref

```
nano mapreads_hybridreference.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_hybridreference.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_hybridreference

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x hybridreference -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_hybridreference.sam -k 5\n; done
```

```
sbatch mapreads_hybridreference.sh
```

--------------------------------------------------------------------------------------------
## Final Hybrid reference transcriptome
- use past.fasta (full) host reference

#### Add suffix to host reference
```
addsuffixtofastaseqnames.py Past past.fasta
```

```
nano addsuffixtofastaseqnames.sh
```

```
#!/bin/bash -l

#SBATCH -o addsuffixtofastaseqnames_past.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=addsuffixtofastaseqnames

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_Kenkel/addsuffixtofastaseqnames.py Past /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/past.fasta
```

```
sbatch addsuffixtofastaseqnames.sh
```

The total number of sequences is 30740
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py past_suffixed.fasta
```

Rename fasta file based on # of contigs:
```
mv past_suffixed.fasta 30740_past_suffixed.fasta
```

--------------------------------------------------------------------------------------------
## Concatenate host reference (Past_full) and symbiont reference to make mapping file
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/hybridref

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/past_suffixed.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/SymTranscriptome_suffixed.fasta > hybridreference.fasta
```

```
ls -alh
```

```
-rwxrwxrwx  1 vradice users 95M Apr 22 15:10 hybridreference.fasta
```

--------------------------------------------------------------------------------------------
## Map all the P_ast samples to concatenated (host plus symbiont) hybridreference transcriptome
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref

#### Make file mappable
need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

bowtie for mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/hybridref/hybridreference.fasta hybridreference
```

#### Mapping to hybridreference
```
nano mapreads_hybridreference.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_hybridreference.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_hybridreference

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x hybridreference -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_hybridref.sam -k 5\n; done
```

```
sbatch mapreads_hybridreference.sh
```

```
head -n 30 bowtie2_hybridreference.txt
```

--------------------------------------------------------------------------------------------
## Need seq2iso tables for host and symbiont
- add Carly's past_seq2iso.tab as a -g argument to the countxpression_SB_advbioinf.py script
- Sandrine (SB) re-wrote the script to be able to merge counts into isogroups at the same time as counting if you supply a seq2iso table

#### Use associated seq2iso table (host)
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/pastreoides_may2014

show count of how many lines in table
```
wc -l past_seq2iso.tab 
grep "" -c past_seq2iso.tab 
```
30740 lines in past_seq2iso.tab

#### Add "_Past" suffix to contig names (column 1) and gene names (column 2) and change seq2iso table to tab delimited

Use regular expressions:
> Find:
>
> (\w+) (\w+)
> 
> Replace:
>
> $1_Past\t$2_Past
>
> (Visual Studio Code uses .net syntax)
> 
> Unix syntax equivalent:
>
> \1_Past\t\2_Past

New file:  past_seq2iso_suffixed.tab
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/

#### Create symbiont symbiont_seq2iso table
- symbiont reference transcriptome:  S.microadriaticum_CassKB8.CDS.fna
- genes only
- make table with same genename_sym in column 1, and column 2 (no column headers)
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/

```
grep ">" 42652_SymTranscriptome_suffixed.fasta > sym_seq2iso.tab
```

show count of how many lines in table
```
wc -l sym_seq2iso.tab
```
42652

The sym_seq2iso.tab has one column with suffixed gene names
- for example:  
>Smic_CassKB8.gene1.mRNA1_Sym
>Smic_CassKB8.gene2.mRNA1_Sym
>Smic_CassKB8.gene3.mRNA1_Sym

Write Python script to find and replace
- re.sub(regex, replacement, subject)

> Find: 
> (>)(\w+.\w+.\w+)
> 
> Replace: 
> $2\t$2
> 
> equivalent:  \2\t\2

```
print(re.sub(r'(>)(\w+.\w+.\w+)', r'\2\t\2', line))
```
for loop

New file:  sym_seq2iso_suffixed.txt
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/

```
mv sym_seq2iso_suffixed.txt sym_seq2iso_suffixed.tab
```

--------------------------------------------------------------------------------------------
## Concatenate host and symbiont seq2 iso table
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref

```
cat /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/past_seq2iso_suffixed.tab /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/sym_seq2iso_suffixed.tab > hybrid_seq2iso.tab
```

--------------------------------------------------------------------------------------------
## Count expression 
Generate counts file for all samples
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

```
nano countexpression_hybridref.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_hybridref.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=countexpression_hybridref

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/*_nof_hybridref.sam -g /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/hybrid_seq2iso.tab
```

```
sbatch countexpression_hybridref.sh
```

```
grep "Sym" -c NO_289_La_Pa_yr2_R1_nof_hybridref_counts.txt
```
42652 lines

search and show which lines have "Sym" to see example of matched contigs
```
grep -n "Sym" NO_289_La_Pa_yr2_R1_nof_hybridref_counts.txt
```

```
grep "Past" -c NO_289_La_Pa_yr2_R1_nof_hybridref_counts.txt
```
29422 lines

```
grep "Sym" -c hybrid_seq2iso.tab
```
42652 lines

```
wc -l hybrid_seq2iso.tab
```
73391 lines

/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014
```
wc -l past_seq2iso_suffixed.tab
```
30740 lines
multiple contigs per gene


/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/S_microadriaticum_CassKB8/
```
wc -l sym_seq2iso_suffixed.tab
```
42651 sym_seq2iso_suffixed.tab
maybe this command thinks that the first line is a header 

```
grep "Sym" -c sym_seq2iso_suffixed.tab 
```
42652


/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref
```
head -n 20 NC_291_C_Pa_yr2_R1_nof_hybridref_counts.txt
tail -n 20 NC_291_C_Pa_yr2_R1_nof_hybridref_counts.txt
grep "Sym" -c NC_291_C_Pa_yr2_R1_nof_hybridref_counts.txt
grep -n "Sym" NC_291_C_Pa_yr2_R1_nof_hybridref_counts.txt
```

```
cat match_counts.txt
```

--------------------------------------------------------------------------------------------
## Parse expression to table
hybrid_seq2iso.tab

#### Make genelist.txt
Format should be:

GeneName |
--- |
gene1name |
gene2name |
gene3name |
gene4name |

secure copy from local machine

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

```
nano ParseExpression.sh
```

```
#!/bin/bash -l

#SBATCH -o ParseExpression.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=ParseExpression

/cm/shared/courses/dbarshis/18AdvBioinf/scripts/ParseExpression2BigTable_advbioinf.py genelist.txt merged_hybridref_counts.txt nomatch *_nof_hybridref_counts.txt
```

```
sbatch ParseExpression.sh
```

#### Output 
merged_hybridref_counts.txt

--------------------------------------------------------------------------------------------



***checking symbiont typing results with revised Arif ITS2 database***

## Identify dominant symbiont type
- based on internal transcribed spacer 2 (ITS2) region, a multicopy genetic marker commonly used to analyse Symbiodinium diversity
- map samples against ITS2 BLAST database
- Arif et al. 2014 [https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869)
    - Corrigendum, with updated ITS2 database file, published 17 July 2019 [https://onlinelibrary.wiley.com/doi/10.1111/mec.14956](https://onlinelibrary.wiley.com/doi/10.1111/mec.14956)
- File:  mec14956-sup-0001-files1_corrigendum.fasta

scp mec14956-sup-0001-files1_corrigendum.fasta vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

#### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py mec14956-sup-0001-files1_corrigendum.fasta
```

```
The total number of sequences is 400
The average sequence length is 376
The total number of bases is 150400
The minimum sequence length is 376
The maximum sequence length is 376
The N50 is 376
Median Length = 376
contigs < 150bp = 0
contigs >= 500bp = 0
contigs >= 1000bp = 0
contigs >= 2000bp = 0
```

--------------------------------------------------------------------------------------------
## Map all samples to Arif ITS2 (all clades) reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

### Make file mappable
Need bowtie build module - creates 6 files for mapping
```
module load bowtie2/2.2.4
```

Bowtie files for Arif ITS2 clade mapping
```
bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/mec14956-sup-0001-files1_corrigendum.fasta Arif_ITS2_corrigendum
```

##### Mapping to Arif_ITS2_corrigendum reference

```
nano mapreads_Arif_ITS2_corrigendum.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Arif_ITS2_corrigendum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=bowtie2_Arif_ITS2_corrigendum

module load bowtie2/2.2.4

for i in *_clippedtrimmed_nofilter.fastq ; do bowtie2 --rg-id ${i%_clippedtrimmed_nofilter.fastq} \
--rg SM:${i%_clippedtrimmed_nofilter.fastq} \
--local -x Arif_ITS2_corrigendum -U $i \
> ${i%_clippedtrimmed_nofilter.fastq}_nof_Arif_ITS2_corrigendum.sam -k 5\n; done
```

```
sbatch mapreads_Arif_ITS2_corrigendum.sh
```

--------------------------------------------------------------------------------------------


--------------------------------------------------------------------------------------------
## Count expression - all reads mapped to Arif_ITS2 symbiont reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

```
nano countexpression_Arif_ITS2_corrigendum.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Arif_ITS2_corrigendum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_Arif_ITS2_corrigendum_Past

/cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/*_Arif_ITS2_corrigendum.sam
```

```
sbatch countexpression_Arif_ITS2_corrigendum.sh
```

#### copy output into local folder on computer
scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/match_counts.txt ./


--------------------------------------------------------------------------------------------
## Merge all *nof_Arif_ITS2_counts.txt files into one big table 
- first need to add column to each .txt file with unique sample id
- so we can later identify which sample has which ITS2 sequences

Using regular expressions:

for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done

- For each file in the list, this will use sed to append to the end of each line a tab and the filename
- Using the -i flag with sed to perform a replacement in-place, overwriting the file
- Perform a substitution with s/PATTERN/REPLACEMENT/. 
- In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), 
- and $f is the filename, from the loop variable. 
- The s/// command is within double-quotes so that the shell can expand variables
- sed is most practical for pattern substitution and saving in-place. 


```
nano append_filename.sh
```

```
#!/bin/bash -l

#SBATCH -o append_filename.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=append_filename

for f in *_Arif_ITS2_corrigendum_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

```
sbatch append_filename.sh
```

#### Concatenate files
```
cat *_nof_Arif_ITS2_corrigendum_counts.txt > merged_Past_ITS2_counts_Arif-corrigendum.txt
```

scp vradice@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/merged_Past_ITS2_counts_Arif-corrigendum.txt ./

## Outcome symbiont clade mapping
- 86% reads mapped to I sequences
  - BLAST confirmed 100% match to former Symbiodinium sp. partial rRNA genes and ITS2
- After removing mapping to I sequences
  - 96% mapped to A


--------------------------------------------------------------------------------------------

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/mapped_hybridref/

Total number of reads per file
```
grep -c @K00364 *.fastq
```

Number of singly aligned reads
```
grep "aligned exactly 1 time" bowtie2_hybridreference.txt
```

Host genes from reference
```
grep -c _Past merged_hybridref_counts.txt
```
29422

Symbiont genes from reference
```
grep -c _Sym merged_hybridref_counts.txt
```
42652


--------------------------------------------------------------------------------------------

