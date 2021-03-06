#!/usr/bin/python

import argparse
import os.path
import sys

def main():

    parser = argparse.ArgumentParser(description='''The purpose of this script is to fix multiple alleles''')

    parser.add_argument('var_file', type=str)
    parser.add_argument('output_file', type=str)
    parser.add_argument("--fix_indels",choices=["True","False"],default="False")

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
                if len(ref)>len(alt) and args.fix_indels=="True":
                    ref=ref[len(alt):]
                    position=str(position+len(alt))
                    alt="-"
                    data[1]=position
                    data[3:5]=[ref,alt]
                    output.write(("\t").join(data))
                elif len(ref) != len(alt) or len(ref)==len(alt)==1:
                    output.write(line)
                else:
                    alt_list=list(alt)
                    ref_list=list(ref)
                    #getting indices of the same item
                    comparison=sameindex(alt,ref)
                    for number in comparison:
                        del alt_list[number]
                        del ref_list[number]
                    for i, (item1,item2) in enumerate(zip(alt_list,ref_list)):
                        #We have to separate each alternate allele
                        new_index=list(find_all(alt,item1))
                        other_index=list(find_all(ref,item2))
                        index_final=list(set(other_index).intersection(new_index))
                        if len(index_final)==1:
                            new_position=str(position+index_final[0])
                            data[1]=new_position
                            data[3:5]=[item2,item1]
                            output.write(("\t").join(data))
    output.close()

def sameindex(string1, string2):
    number_list=[]
    for i, (char1, char2) in enumerate(zip(string1, string2)):
        if char1 == char2:
            number_list.append(i)
    number_list.sort(reverse=True)
    return number_list
    
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

if __name__ == '__main__':
    main()



