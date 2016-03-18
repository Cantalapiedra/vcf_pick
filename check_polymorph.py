#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys

fields = {}
other_gtps = {}
parent_1_i=-1
parent_2_i=-1
total_variants = 0
total_missing = 0
total_polymorphic = 0
total_homozygots = 0
total_indels = 0
total_snps = 0

# Checks if a genotype is hetero
def is_hetero(genotype):
    retvalue = False
    alleles = genotype.split("/")
    
    if alleles[0] == alleles[1]: retvalue = False
    else: retvalue = True
    
    return retvalue

# Checks if the variant of a genotype is reference-like or alternative
def parse_var(parent_gtp, ref, alt):
    retvalue = parent_gtp
    alleles = parent_gtp.split("/")
    retlist = []
    
    for allele in alleles:
        if allele == "0":
            retlist.append(ref)
        elif allele == "1":
            retlist.append(alt)
        else:
            retlist.append(allele)
    
    retvalue = "/".join(retlist)
    
    return retvalue

# Parse other genotypes to output stats
# TODO: input arguments could be a list of parent genotypes
# and a list of other lines to parse to obtain the stats
# The genotypes not in either list would be ignored (e.g.: other parentals)
def parse_others(data, other_gtps, parent_1_gtp, parent_2_gtp):
    num_miss = 0
    num_het = 0
    num_hom = 0
    parent_1_like = 0
    parent_2_like = 0
    alike = 0
    total_dp = 0
    avg_dp = 0
    
    for k in other_gtps:
        #print data[k]
        k_data = data[k].split(":")
        k_gtp = k_data[0]
        
        if "." in k_gtp: num_miss += 1
        elif is_hetero(k_gtp): num_het+=1
        else: num_hom+=1
        
        if k_gtp == parent_1_gtp: parent_1_like += 1
        elif k_gtp == parent_2_gtp: parent_2_like += 1
        else: alike += 1
        
        total_dp += int(k_data[2])
    
    avg_dp = total_dp / len(other_gtps)
    
    return (num_miss, num_hom, num_het, parent_1_like, parent_2_like, alike, total_dp, avg_dp)

inputfilename=sys.argv[1]
parent_1=sys.argv[2]
parent_2=sys.argv[3]
allow_het=sys.argv[4] # yes no
allow_miss="no"

sys.stderr.write("Parents to compare: "+parent_1+" and "+parent_2+"\n")
sys.stderr.write("parse_vcf: Processing "+inputfilename+"\n")

for i, line in enumerate(open(inputfilename, 'r')):
    # Header
    if (i==0):
        data = line[1:].strip().split()
        for j, field in enumerate(data):
            #sys.stderr.write(str(j)+" is "+field+"\n")
            fields[j] = field
            if (field==parent_1): parent_1_i=j
            elif (field==parent_2): parent_2_i=j
            elif (j>8): other_gtps[j] = field
        
        #sys.stderr.write("Parent "+parent_1+" is field "+str(parent_1_i)+"\n")
        #sys.stderr.write("Parent "+parent_2+" is field "+str(parent_2_i)+"\n")
        #sys.stderr.write("Other genotypes "+str(other_gtps.keys())+"\n")
        
        sys.stdout.write(parent_1+"\t"+parent_1+"_cov\t"+parent_1+"_het\t"+\
                         parent_2+"\t"+parent_2+"_cov\t"+parent_2+"_het\t"+\
                         "contig\tpos\tREF\tALT\tqual\ttype\n")
    # Variant
    else:
        total_variants += 1
        data = line.strip().split()
        #print data
        parent_1_data = data[parent_1_i].split(":")
        parent_2_data = data[parent_2_i].split(":")
        parent_1_gtp = parent_1_data[0]
        parent_2_gtp = parent_2_data[0]
        
        if "." in parent_1_gtp or "." in parent_2_gtp:
            total_missing += 1
            if allow_miss == "no": continue
        
        if parent_1_gtp == parent_2_gtp: continue
        
        total_polymorphic += 1
        
        # Check heteros
        parent_1_het = is_hetero(parent_1_gtp)
        parent_2_het = is_hetero(parent_2_gtp)
        if parent_1_het or parent_2_het:
            if allow_het == "no": continue
        else:
            total_homozygots += 1
        
        parent_1_dp = parent_1_data[2]
        parent_2_dp = parent_2_data[2]
        
        ## Remaining information
        chrom = data[0]
        pos = data[1]
        ref = data[3]
        alt = data[4]
        qual = data[5]
        if "INDEL" in data[7]:
            vartype = "indel"
            total_indels += 1
        else:
            vartype = "snp"
            total_snps += 1
            
        parent_1_var = parse_var(parent_1_gtp, ref, alt)
        parent_2_var = parse_var(parent_2_gtp, ref, alt)
        
        (num_miss, num_hom, num_het, parent_1_like, parent_2_like, alike, total_dp, avg_dp) \
        = parse_others(data, other_gtps, parent_1_gtp, parent_2_gtp)
        
        sys.stdout.write(parent_1_var+"\t"+parent_1_dp+"\t"+str(parent_1_het)+\
                         "\t"+parent_2_var+"\t"+parent_2_dp+"\t"+str(parent_2_het)+\
                         "\t"+chrom+"\t"+pos+"\t"+ref+"\t"+alt+"\t"+qual+"\t"+vartype+\
                         "\t"+str(num_miss)+"\t"+str(num_hom)+"\t"+str(num_het)+\
                        "\t"+str(parent_1_like)+"\t"+str(parent_2_like)+"\t"+str(alike)+\
                        "\t"+str(total_dp)+"\t"+str(avg_dp)+"\n")
    
sys.stderr.write("\nTotal variants: "+str(total_variants)+"\n")
sys.stderr.write("Total missing: "+str(total_missing)+"\n")
sys.stderr.write("Total polymorphic: "+str(total_polymorphic)+"\n")
sys.stderr.write("Total homozygots: "+str(total_homozygots)+"\n")
sys.stderr.write("Total SNPs: "+str(total_snps)+"\n")
sys.stderr.write("Total INDELs: "+str(total_indels)+"\n")
sys.stderr.write("\nFinished.\n")
