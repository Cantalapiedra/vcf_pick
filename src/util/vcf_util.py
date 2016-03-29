#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

# VCF version tested
##fileformat=VCFv4.1
VCF_CONTIG_COL = 0
VCF_POS_COL = 1
VCF_REF_COL = 3
VCF_ALT_COL = 4
VCF_INFO_COL = 7
VCF_SAMPLES_INI_COL = 9
# snpEff version tested
##SnpEffVersion="4.0b (build 2014-08-28), by Pablo Cingolani"
SNPEFF_FIELD = "EFF="
SNPEFF_EFF_SEP = ","
SNPEFF_AA_POS = 3
SNPEFF_GENE_POS = 5
SNPEFF_ISOF_POS = 8
# symbols
HETERO = "h"
MISS = "."
REF = "0"

def parse_header(line_data, genotypes_dict, \
                 samples_filename = "", samples_list = [], samples_translation = "", samples_trans_dict = {}):
    
    for j, genotype in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
        new_genotype = {"name":genotype, "index":j}
        if samples_translation == "":
            new_genotype["good_name"] = genotype
        else:
            new_genotype["good_name"] = samples_trans_dict[genotype]
        
        if samples_filename == "":
            if samples_translation == "":
                samples_list.append(genotype)
            else:samples_list.append(samples_trans_dict[genotype])
                
        else:
            if new_genotype["good_name"] not in set(samples_list): continue
        
        genotypes_dict[j] = new_genotype
    
    return

##

def parse_vcf_header_file(vcf_header_file, genotypes_dict, \
                          samples_filename, samples_list, samples_translation, samples_trans_dict):
    header_found = False
    # If a headers file is specified, record names of genotypes
    if vcf_header_file == "":
        pass
    else:
        try:
            headers_file = open(vcf_header_file, 'r')
            for header in headers_file:
                line_data = header.strip().split()
                parse_header(line_data, genotypes_dict, \
                             samples_filename, samples_list, samples_translation, samples_trans_dict)
                header_found = True
                break
            
        except Exception:
            raise
        finally:
            headers_file.close()
    
    return header_found


##

def parse_alleles(line_data, genotypes_dict):
    alleles = {}
    
    for j, line_line in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
        
        if j not in genotypes_dict: continue
        
        genotype = line_line.split(":")[0]
        
        if len(genotype)==3:
            if genotype[0] != genotype[2]:
                allele = HETERO
            elif genotype[0] == genotype[2]:
                allele = genotype[0]
            else:
                allele = genotype
            
        else:
            raise Exception("Genotype with more fields than expected: "+genotype+".")
        
        if allele in alleles:
            alleles[allele].append(j)
        else:
            alleles[allele] = [j]
    
    return alleles

##

def calculate_maf(alleles, num_genotypes):
    maf = 1.0
    
    for allele in alleles:
        if allele == HETERO or allele == MISS:
            continue
        else:
            curr = len(alleles[allele])
            if allele == REF: curr += 1
            curr = curr * 1.0 / num_genotypes
            
            if curr < maf: maf = curr
    
    return maf

##

 # Heteros --> missing; missing --> rm; monomorph --> rm
def preprocess_variant(alleles, max_heteros, max_missing, show_monomorph, maf = 0.3):
    retValue = True
    
    #print str(alleles)+"\t"+str(retValue)
    
    alleles_keys = set(alleles)
    
    if len(alleles_keys) <= 1 and not show_monomorph:
            retValue = False
    else:
        num_genotypes = sum([len(allele_genotypes) for allele_genotypes in alleles.values()])+1
        
        ## Heteros --> missing
        if HETERO in alleles_keys:
            percent_heteros = len(alleles[HETERO]) * 1.0 / num_genotypes
        else:
            percent_heteros = 0
        
        if percent_heteros > max_heteros:
            if MISS not in alleles:
                alleles[MISS] = []
            alleles[MISS] = alleles[MISS] + alleles[HETERO]
            del alleles[HETERO]
            alleles_keys = set(alleles)
            num_genotypes = sum([len(allele_genotypes) for allele_genotypes in alleles.values()])+1
        
        ## Missing -> rm
        if MISS in alleles_keys:
            percent_miss = len(alleles[MISS]) * 1.0 / num_genotypes
        else:
            percent_miss = 0
        
        if percent_miss > max_missing:
            retValue = False
        elif len(alleles_keys) <= 1 and not show_monomorph:
            retValue = False
        else:
            ## MAF
            alleles_maf = calculate_maf(alleles, num_genotypes)
            if alleles_maf < maf:
                retValue = False
    
    #print "\t"+str(alleles)+"\t"+str(retValue)
    
    return retValue

##