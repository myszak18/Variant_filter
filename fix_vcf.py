# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 14:36:30 2015

@author: Myszak
"""

import argparse
import os.path
import sys

def main():
    
    parser = argparse.ArgumentParser(description='''The purpose of this script is to fix multiple alleles''')
    
    parser.add_argument('var_file', type=str)
    parser.add_argument('output_file', type=str)
    parser.add_argument("--fix_indels",choices=[True,False],default=False)
    
    args = parser.parse_args()
    
    if os.path.isfile(args.var_file)==False:
        sys.exit("var-file is not valid!")
        
    output=open(args.output_file,"w")
    
    var_file=open(args.var_file)
    
    content_vcf=var_file.readlines()
    
    var_file.close()
    
    for line in content_vcf:
        if line.startswith("#"):
            output.write(line)
        else:
            data=line.split("\t")
            ## VCF file must contain at least 5 fields: chr position id ref alt
            if len(data)<5:
                sys.exit("VCF file not valid")
            else:
                ref=data[3]
                alt=data[4]
                position=int(data[1])
                ## Fix indels if this option is set to true
                if args.fix_indels:
                    if len(ref)>len(alt):
                        ref=alt.lstrip(alt)
                        position=str(position+len(alt))
                        alt="-"
                        data[1]=position
                        data[3:5]=[ref,alt]
                else:
                    if len(ref) != len(alt) or len(ref)==len(alt)==1:
                        output.write(line)
                    