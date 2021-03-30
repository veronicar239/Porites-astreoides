# Longest isoform (isotig) per gene (isogroup)
## Carly Kenkel's *Porites astreoides* reference

### Standard filtering criteria of transcriptome assemblies

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 30740_past_suffixed.fasta
```

*Reference assembly has <500 bp contigs, and multiple isoforms per isogroup*
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


----------------------------------------------------------------------------------------------
*Original reference transcriptome assembly*
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/pastreoides_may2014/

### first check file and edit to tab delimited if needed

*past_seq2iso.tab copy is space delimited*

make into /t delimited .txt file

```
cat past_seq2iso.txt | tr ' ' '\t' > past_seq2iso_tDelim.txt
```

### filter longest contigs

```
nano longestContig.sh
```

```
#!/bin/bash -l

#SBATCH -o longestContig.txt
#SBATCH -n 6
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=longestContig

python2 /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Porites_astreoides/Paytan/P_ast/refassembly/P_ast_full_Kenkel_2014/pastreoides_may2014/fasta_longest_contig_per_gene_seq2genetable.py --infasta 30740_past.fasta --seq2gene past_seq2iso_tDelim.txt --outfilename past_LongestContig.fasta
```

```
sbatch longestContig.sh
```

```
grep ">" -c past_LongestContig.fasta
```
29422

**Best practice to rename with contig number**
```
cp past_LongestContig.fasta 29422_past_LongestContig.fasta
```

#### Check assembly
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




----------------------------------------------------------------------------------------------
Alternate

Another way to get isoform sequence length information


#### First need to get length information about isotigs, the select longest isotig per isogroup.

```
grep ">" past_CDS.fas > past_CDS_headers.txt
```

Need: 
- longest isotig of each isogroup
- "GE xxxxxx" isogroup (there are no duplicates and each was linked to its own gene)

```
grep ">" past.fasta > past_headers.txt
```

grep ">" -c past.fasta
30740

Match contigs from past.fasta with longest isoform contig list


#### extract subset of sequences based on list of names in .txt file

Usage:
```
seqtk subseq   genes.fasta  subsetIDs.txt   > genes_subset.fastq
```

module load seqtk/1.0
[seqtk](https://github.com/lh3/seqtk)
```
enable_lmod
module load container_env seqtk
seqtk subseq past.fasta Past_longestIso_headers.txt > past_longestIso.fasta
```

