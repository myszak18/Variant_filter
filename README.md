# Variant_filter

This scrpit is based on variant filter scrpit fpfilter.pl by Dan Koboldt and David Larson (https://github.com/ckandoth/variant-filter).
It adds new functionality: filtering insertions and deletion.

# USAGE: #

## 1. Vcf preparation ##

a) Script accepts variant calls in vcf format. It assumes that vcf contains one allele per line. If it is not the case,
one can format vcf with vcfbreakmulti utility from vcflib (https://github.com/ekg/vcflib).

b) Double allele issue: if there is both single nucleotide variant and indel present in one position, after using vcfbreakmulti the line with single nucleotide variant may look like this:

chr1 8789695 . AC AT. Such cases must be fixed (the script to do so is in preparation) like this:  
chr1 8789696 . C  T. 

c) Insertion and deletion format.
Script assumes insetions are present in this format:  
chr1 8789695 . A AT  
Deletions can be present in either:  
- 0 format: chr2 5858778 . AC A  
- 1 format: chr2 5858779 . C - .  
Readcount file preparation is more complicated if deletions are present in 0 format (see readcount file preparation below). 

## 2. Readcount file preparation. ##

a) Variant position preperation (var file)

If you have 1 format based indel, variant position extraction is just the same as in original fpfilter ie:   

`perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' snvs.vcf > snvs.var`  

However, if you have 0 format indels, you have to their position and change it to 1 format (the script to do so is in preparation) so that chr2 5858778 . AC A in your file becomes chr2 5858779 5858779 in var file.  

b) Readcount file preparation

Install readcount script: https://github.com/genome/bam-readcount 

Run it like follows:
`bam-readcount -q1 -b15 -w1 -l snvs.var -f reference.fasta example.bam > snvs.readcount`

## 3. Script example usage: ##

a) basic usage:

`python variant_filter.py example.vcf snvs.readcount --output_file new.fpfilter`  
This assumes indels are in 0 format coordinates. Output file name is optional.

b) changing to 1 format coordinates:

`python variant_filter.py example.vcf snvs.readcount --output_file new.fpfilter --indel_type 1_format`

c) changing filter parameters:

The default parameteres are like in original fpfilter.pl script. You can change all of them for example:
`python variant_filter.py example.vcf snvs.readcount --output_file new.fpfilter --min_depth 20 ` changes minimal variant depth to 20.
