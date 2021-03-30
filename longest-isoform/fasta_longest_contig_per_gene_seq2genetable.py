#! /usr/bin/env python

import sys, os, argparse
# Extracts the longest contig per transcript/gene from a fasta file

def read_fasta(file): #returns a dictionary with the seq name as the key and the sequence as the entry/item
	fin = open(file, 'r')
	count=0
	
	contigs={}
	seq=''
	for line in fin:
		line=line.strip()
		if line and line[0] == '>':                #indicates the name of the sequence
			count+=1
			if count>1:
				contigs[name]+=seq
			name=line[1:]
			contigs[name]=''
			seq=''
		else:
			seq +=line
	fin.close()
	return contigs

def main():
	#to parse command line
	usage = "Extracts the longest contig per transcript/gene from a fasta file, written by DJ Barshis"
	p = argparse.ArgumentParser(usage)


	#Input/output files
	p.add_argument('-i', '--infasta', help='Assembly.fasta to parse for longest contig per gene/transcript')
	p.add_argument('-s', '--seq2gene', help='Sequence to gene table to use for parsing, tab delimited, sequence then gene, no header')
	p.add_argument('-o', '--outfilename', default='LongestContig.fasta', help='New, parsed fastafile with longest Contig per gene/transcript [default LongestContig.fasta]')

	#General
	args = p.parse_args()

	#Parsing seqtable first column sequence or contigname, second column genename
	seqtable = open(args.seq2gene, 'r')
	Seq2GeneDict={}
	for Line in seqtable:
		Line=Line.rstrip()
		cols=Line.split('\t')
		Seq2GeneDict[cols[0]]=cols[1]
	print('Read in seq2gene table')
	seqtable.close()

	#Parsing fasta file into a dictionary with the sequence name as the key and the sequence as the item/entry
	FastaContents = read_fasta(args.infasta)
	print('Read in fasta file to parse')

	#Generating longest Sequence/contig dictionary with the gene as the key and a list with the [contig name, sequence length]
	LongestContigs = {}

	for Seqname in FastaContents.keys():
		Shortname=Seqname.split(' ')[0]
		Length = len(FastaContents[Seqname])
		try:
			if LongestContigs[Seq2GeneDict[Shortname]][1]<Length:
					LongestContigs[Seq2GeneDict[Shortname]]=[Seqname, Length]
		except KeyError:
			LongestContigs[Seq2GeneDict[Shortname]]=[Seqname, Length]
	print('created longest Gene dictionary')
	print('now writing out longest isoform .fasta file')

	outfile = open(args.outfilename, 'w')
	for Gene in LongestContigs.keys():
		outfile.write('>%s\n%s\n'%(LongestContigs[Gene][0],FastaContents[LongestContigs[Gene][0]]))
	outfile.close()

if __name__ == '__main__':
	main()