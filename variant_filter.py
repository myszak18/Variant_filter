#!/usr/bin/python

"""
Created on Fri Sep 04 10:54:59 2015

@author: Myszak
"""

import argparse
import os.path
import sys

def main():

    parser = argparse.ArgumentParser(description='''Filtering variants, including insertions and deletions.
    This script assumes that there is only one variant call per line. If you have multiple alleles per line,
    separate them with vcfbreakmulti utility from vcflib (https://github.com/ekg/vcflib).''')

    min_depth = 8

    # Minimums for variant allele fraction, and the number of variant-supporting reads
    min_var_frac = 0.05
    min_var_count = 3
    # Minimum avg relative distance of variant from start/end of read
    min_read_pos = 0.10
    # Minimum representation of variant allele on each strand
    min_strandedness = 0.01
    # Maximum difference of mismatch quality sum between var/ref reads (paralog filter)
    max_mmqs_diff = 50
    # Maximum mismatch quality sum of reference-supporting reads
    max_var_mmqs = 100
    # Maximum difference of mapping quality between variant and reference reads
    max_mapqual_diff = 30
    # Maximum difference of average supporting read length between variant and reference reads (paralog filter)
    max_readlen_diff = 25
    # Minimum average distance to effective 3prime end of read (real end or Q2) for variant-supporting reads
    min_var_dist_3 = 0.2

    parser.add_argument('var_file', type=str)
    parser.add_argument('readcount_file', type=str)
    parser.add_argument('--output_file', type=str,default="snvs.fpfilter")
    parser.add_argument("--indel_type",choices=["0_format","1_format"],default="0_format",
                        help="""0_format example chr2 5858778 . AC A \n
                        1_format example chr2 5858779 . C -""")


    parser.add_argument('--min_depth', default=min_depth, type=int)
    parser.add_argument('--min_var_frac', default=min_var_frac, type=float)
    parser.add_argument('--min_var_count', default=min_var_count, type=int)
    parser.add_argument('--min_read_pos', default=min_read_pos, type=float)
    parser.add_argument('--min_strandedness', default=min_strandedness, type=float)
    parser.add_argument('--max_mmqs_diff', default=max_mmqs_diff, type=int)
    parser.add_argument('--max_var_mmqs', default=max_var_mmqs, type=int)
    parser.add_argument('--max_mapqual_diff', default=max_mapqual_diff, type=int)
    parser.add_argument('--max_readlen_diff', default=max_readlen_diff, type=int)
    parser.add_argument('--min_var_dist_3', default=min_var_dist_3, type=float)

    args = parser.parse_args()

    if os.path.isfile(args.var_file)==False and os.path.isfile(args.readcount_file)==False:
        sys.exit("var-file and readcount-file not valid!")
    elif os.path.isfile(args.var_file)==False:
        sys.exit("var-file not valid")
    elif os.path.isfile(args.readcount_file)==False:
        sys.exit("readcount-file not valid")

    readcounts_by_position={}

    readcount=open(args.readcount_file)

    readcount_content=readcount.readlines()

    readcount.close()

    ## Load the read counts into a hash for quick lookup ##

    for line in readcount_content:
        [chrom,position]=line.strip("\n").split("\t")[0:2]
        key="%s:%s" %(chrom,position)
        readcounts_by_position[key]=line.strip("\n")

    # Open the output file for writing, and write a header line with column names
    output=open(args.output_file,"w")

    output.write("#CHROM\tPOS\tREF\tVAR\tDEPTH\tRAF\tVAF\tFILTER\tFILTER_DETAILS\n")

    # Initialize all variant fail/pass counters to zero

    counters=["num_variants","num_fail_depth","num_fail_pos","num_fail_strand","num_fail_varcount",
        "num_fail_varfrac","num_fail_mmqs","num_fail_var_mmqs","num_fail_mapqual","num_fail_readlen",
        "num_fail_dist3","num_no_readcounts","num_pass_filter"]

    stats= {key: 0 for key in counters}

    #Parse vcf file

    vcf=open(args.var_file)

    vcf_content=vcf.readlines()

    vcf.close()

    for line in vcf_content:
        # Skip comment lines
        if line.startswith("#"):
            pass
        else:
            line_list=[]
            fields=line.split("\t")
            [chrom,position,id_var,ref,alt]=fields[0:5]
            #For SNVs we can only get reads for position listed first in vcf file
            if len(ref)==len(alt)!=1 and alt!="-":
                ref=ref[0]
                alt=alt[0]
            elif len(ref)<len(alt):
                alt="+"+alt[len(ref):]
            elif len(ref)>len(alt) or alt=="-":
                if args.indel_type=="0_format":
                    alt="-"+ref[len(alt):]
                    #Our readcount file should take into account position change in indels
                    position=str(int(position)+1)
                elif args.indel_type=="1_format":
                    alt="-"+ref
            stats["num_variants"]=stats["num_variants"]+1
            key="%s:%s" %(chrom,position)
            readcounts=readcounts_by_position.get(key,0)
            if readcounts!=0:
                if "-" in alt:
                    #Reference base change is necessary in case on indels since the position is different that the one in vcf file
                    ref=readcounts.split("\t")[2]
                ref_result=read_counts_by_allele(readcounts,ref)
                alt_result=read_counts_by_allele(readcounts,alt)
                # Proceed only if readcounts are available for reference and variant alleles
                if(ref_result and alt_result):
                    [total_depth,ref_count, ref_map_qual, ref_base_qual, ref_semq, ref_plus, ref_minus, ref_pos,
                                 ref_subs, ref_mmqs, ref_q2_reads, ref_q2_dist, ref_avg_rl, ref_dist_3]=[float(x) for x in ref_result.split("\t")]
                    [var_count, var_map_qual, var_base_qual, var_semq, var_plus, var_minus, var_pos,
                        var_subs, var_mmqs, var_q2_reads, var_q2_dist, var_avg_rl, var_dist_3]=[float(x) for x in alt_result.split("\t")[1:]]
                    if total_depth>0:
                        mismatch_qualsum_diff = var_mmqs - ref_mmqs
                        # Determine map qual diff between ref/var reads
                        mapqual_diff = ref_map_qual - var_map_qual
                        # Determine difference in average supporting read length
                        readlen_diff = ref_avg_rl - var_avg_rl
                        # Set max strandedness cutoff
                        max_strandedness = 1 - args.min_strandedness
                        vaf="{0:.4f}".format(var_count/(total_depth))
                        raf="{0:.4f}".format(ref_count/(total_depth))
                        line_list=[chrom,position,ref,alt,str(total_depth),raf,vaf]
                        # Set conservative default values for reference and variant strandedness
                        [ref_strandedness, var_strandedness] =[0.5, 0.5]
                        # Determine reference strandedness
                        if ref_plus + ref_minus  > 0:
                            ref_strandedness = ref_plus / ( ref_plus + ref_minus )
                        # Determine variant strandedness
                        if var_plus + var_minus  > 0:
                            var_strandedness = var_plus / ( var_plus + var_minus )
                        if (var_count>0):
                                ## FAILURE: Average distance of variant from clipped read ends ##
                            if(var_pos < args.min_read_pos):
                                comment="%.2f < %.2f\n" %(var_pos,args.min_read_pos)
                                line_list.extend(["ReadPos",comment])
                                stats["num_fail_pos"]=stats["num_fail_pos"]+1
                                ## FAILURE: Variant is strand-specific but reference is NOT strand-specific ##
                            elif((var_strandedness < args.min_strandedness or var_strandedness > max_strandedness) and (ref_strandedness >= args.min_strandedness and ref_strandedness <= max_strandedness)):
                                comment="Ref=%.2f Var=%.2f MinMax=[%.2f,%.2f]\n" %(ref_strandedness,var_strandedness,args.min_strandedness,max_strandedness)
                                line_list.extend(["Strandedness",comment])
                                stats["num_fail_strand"]=stats["num_fail_strand"]+1
                                ## FAILURE: Variant allele count does not meet minimum ##
                            elif(var_count < args.min_var_count):
                                comment="%.2f < %.2f\n" %(var_count,args.min_var_count)
                                line_list.extend(["VarCount",comment])
                                stats["num_fail_varcount"]=stats["num_fail_varcount"]+1
                                ## FAILURE: Read depth is too low to proceed onto next few filters ##
                            elif(total_depth < args.min_depth):
                                comment="%.2f < %.2f\n" %(total_depth,args.min_depth)
                                line_list.extend(["Low_Depth",comment])
                                stats["num_fail_depth"]=stats["num_fail_depth"]+1
                                ## FAILURE: Variant allele fraction does not meet minimum ##
                            elif(float(vaf) < args.min_var_frac):
                                comment="%.2f < %.2f\n" %(float(vaf),args.min_var_frac)
                                line_list.extend(["VarFrac",comment])
                                stats["num_fail_varfrac"]=stats["num_fail_varfrac"]+1
                                ## FAILURE: Paralog filter for sites where variant allele mismatch-quality-sum is significantly higher than the reference allele MMQS
                            elif( mismatch_qualsum_diff > args.max_mmqs_diff ):
                                comment="%.2f - %.2f = %.2f > %.2f\n" %(var_mmqs, ref_mmqs,mismatch_qualsum_diff,args.max_mmqs_diff)
                                line_list.extend(["MismatchQualsum",comment])
                                stats["num_fail_mmqs"]=stats["num_fail_mmqs"]+1
                                ## FAILURE: Mapping quality difference exceeds allowable maximum ##
                            elif(mapqual_diff > args.max_mapqual_diff):
                               comment="%.2f - %.2f = %.2f > %.2f\n" %(ref_map_qual, var_map_qual,mapqual_diff,args.max_mapqual_diff)
                               line_list.extend(["MapQual",comment])
                               stats["num_fail_mapqual"]=stats["num_fail_mapqual"]+1
                                ## FAILURE: Read length difference exceeds allowable maximum ##
                            elif(readlen_diff > args.max_readlen_diff ):
                                comment="%.2f - %.2f = %.2f > %.2f\n" %(ref_avg_rl, var_avg_rl,readlen_diff,args.max_readlen_diff)
                                line_list.extend(["ReadLen",comment])
                                stats["num_fail_readlen"]=stats["num_fail_readlen"]+1
                                ## FAILURE: Avg distance from 3' ends of reads is lower than allowed minimum ##
                            elif(var_dist_3 < args.min_var_dist_3):
                                comment="%.2f < %.2f\n" %(var_dist_3,args.min_var_dist_3)
                                line_list.extend(["VarDist3",comment])
                                stats["num_fail_dist3"]=stats["num_fail_dist3"]+1
                                ## FAILURE: Mismatch quality sum indicative of errors from misalignment ##
                            elif( var_mmqs > args.max_var_mmqs ):
                                comment="%.2f > %.2f\n" %(var_mmqs,args.max_var_mmqs)
                                line_list.extend(["VarMMQS",comment])
                                stats["num_fail_var_mmqs"]=stats["num_fail_var_mmqs"]+1
                                ## SUCCESS: Passes all filters above ##
                            else:
                                line_list.extend(["PASS"," \n"])
                                stats["num_pass_filter"]=stats["num_pass_filter"]+1
                        else:
                            line_list.extend(["NoVariantReads"," \n"])
                            stats["num_no_readcounts"]=stats["num_no_readcounts"]+1
                    else:
                        line_list=[chrom,position,ref,alt,"0","0","0","NoReadCounts"," \n"]
                        stats["num_no_readcounts"]=stats["num_no_readcounts"]+1
                else:
                    line_list=[chrom,position,ref,alt,"0","0","0","NoReadCounts"," \n"]
                    stats["num_no_readcounts"]=stats["num_no_readcounts"]+1
            else:
                line_list=[chrom,position,ref,alt,"0","0","0","NoReadCounts"," \n"]
                stats["num_no_readcounts"]=stats["num_no_readcounts"]+1
            if len(line_list)>0:
                output.write(("\t").join(line_list))


    output.close()
    print "%i variants" %(stats["num_variants"])
    print "%i had a position near the ends of most supporting reads (position < %.2f)" %(stats["num_fail_pos"], args.min_read_pos)
    print "%i had strandedness < %.2f (most supporting reads are in the same direction)" %(stats["num_fail_strand"], args.min_strandedness)
    print "%i had var_count <%i (not enough supporting reads)" %(stats["num_fail_varcount"],args.min_var_count)
    print "%i had depth < %i " %(stats["num_fail_depth"],args.min_depth)
    print "%i had var_frac < %.2f (low-fraction variants are likely artifacts or from crosstalk between samples in the same lane)" %(stats["num_fail_varfrac"],min_var_frac)
    print "%i had mismatch qualsum difference > %i (likely a result of paralogous misalignments)" %(stats["num_fail_mmqs"],args.max_mmqs_diff)
    print "%i had variant MMQS > %i (likely a result of paralogous misalignments)" %(stats["num_fail_var_mmqs"],args.max_mmqs_diff)
    print "%i had mapping quality difference > %i " %(stats["num_fail_mapqual"],args.max_mapqual_diff)
    print "%i had read length difference > %i " %(stats["num_fail_readlen"],args.max_readlen_diff)
    print "%i had var_distance_to_3' < %.2f (illumina errors are more frequent at the 3' ends of reads)" %(stats["num_fail_dist3"],args.min_var_dist_3)
    print "%i passed all filters" %(stats["num_pass_filter"])
    print "%i had no readcounts for the variant allele\n" %(stats["num_no_readcounts"])





## read_counts_by_allele - Retrieve relevant read counts for a certain allele #
def read_counts_by_allele(line,allele):
    lineContents=line.strip("\n").split("\t")
    total_depth=lineContents[3]
    for item in lineContents[5:]:
        alleleContents=item.split(":")
        if alleleContents[0]==allele:
            numAlleleContents=alleleContents
            if len(numAlleleContents)<8:
                return("")
            else:
                numAlleleString=("\t").join(numAlleleContents[1:])
                return_string="%s\t%s" %(total_depth,numAlleleString)
                return(return_string)



if __name__ == '__main__':
    main()





