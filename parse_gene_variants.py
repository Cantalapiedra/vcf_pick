#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys, traceback
from optparse import OptionParser

GENES = 0
ISOFS = 1

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

def _print_parameters(options):
    options_dict = options.__dict__
    sys.stderr.write("Options:\n")
    for option in options_dict:
        if options_dict[option]:
            sys.stderr.write("\t"+option+"="+options_dict[option]+"\n")
    
    return

## Argument parsing
__usage = "usage: parse_variants.py [VCF_FILE] [OPTIONS]"

optParser = OptionParser(__usage)

optParser.add_option('-H', '--header', action='store', dest='vcf_header', type='string', \
                     help='VCF header if VCF_FILE does not include header.')

optParser.add_option('-g', '--genes_list', action='store', dest='genes_list', type='string', \
                     help='')

optParser.add_option('-i', '--isof_list', action='store', dest='isof_list', type='string', \
                     help='')

optParser.add_option('-s', '--samples_translation', action='store', dest='samples_translation', type='string', \
                     help='')

optParser.add_option('-f', '--output_format', action='store', dest='output_format', type='string', \
                     help='summary or tabular')

#optParser.add_option('-v', '--verbose', action='store_true', dest='verbose', help='More information printed.')

(options, arguments) = optParser.parse_args()

if not arguments or len(arguments)==0:
    optParser.exit(0, "You may wish to run '-help' option.\n")

vcf_filename = arguments[0] # This is mandatory

if options.vcf_header: vcf_header = options.vcf_header
else: vcf_header = ""

if options.genes_list:
    if options.isof_list:
        raise Exception("Only a list of genes or a list of isoforms a list of contigs should be provided, not both.")
    query_file = options.genes_list
    query_type = GENES
elif options.isof_list:
    query_file = options.isof_list
    query_type = ISOFS
else: query_file = ""

if query_file == "": raise Exception("Either a list of genes or isoforms is needed.")

if options.samples_translation: samples_translation = options.samples_translation
else: samples_translation = ""

if options.output_format: output_format = options.output_format
else: output_format = ""

###############
###############

def print_variants(variants_dict, genotypes_dict, alleles_dict, output_fmt = "tabular"):
    
    # Prepare output lines
    outputs_list = []
    for variant_id in variants_dict:
        variant = variants_dict[variant_id]
        
        output_list = []
        for eff in variant["effs"]:
            if query_type == GENES:
                output_list = [variant_id, eff["eff_isof"], eff["eff_type"], \
                                variant["ref"], variant["alt"], eff["eff_aa"], \
                                variant["contig"], variant["pos"]]
            elif query_type == ISOFS:
                output_list = [variant_id, eff["eff_isof"], eff["eff_type"], \
                                variant["ref"], variant["alt"], eff["eff_aa"], \
                                variant["contig"], variant["pos"]]
            else:
                sys.stderr.write("Unrecognized query type "+query_type+"\n")
            
            outputs_list.append(output_list)
        
    # Sort the list of variants to output
    # by isoform, contig and position
    outputs_list = sorted(outputs_list, key=lambda x: (x[1], x[6], int(x[7])))
    
    alleles_list = sorted([allele for allele in alleles_dict], \
        cmp=lambda x,y: -1 if x.isdigit() and y.isdigit() and x<y \
                    else 1 if x.isdigit() and y.isdigit() and x>=y \
                    else -1 if x.isdigit() and (not y.isdigit()) \
                    else 1 if (not x.isdigit) and y.isdigit() else 0)   
    
    ### Print output
    genotypes_index_list = sorted(genotypes_dict, key=lambda x: genotypes_dict[x]["good_name"])
    
    # Header
    sys.stdout.write("#\tisof\teffect\tref\talt\tchange\tcontig\tpos")
    if output_fmt == "tabular":
        for genotype in genotypes_index_list:
            sys.stdout.write("\t"+genotypes_dict[genotype]["good_name"]+"\t")
            
    sys.stdout.write("\n")
    
    # Data
    for output_line in outputs_list:
        var_id = output_line[0]
        isof_id = output_line[1]
        fields = [str(field) for field in output_line[1:]]
        sys.stdout.write(">\t"+"\t".join(fields))
        
        if output_fmt == "summary":
            sys.stdout.write("\n")
            
        elif output_fmt == "detail":
            sys.stdout.write("\n")
            variant_alleles = variants_dict[var_id]['alleles']
            
            for allele in (set(alleles_list) & set(variant_alleles)):
                sys.stdout.write(">>\t")
                sys.stdout.write(allele+"\t")
                genotypes_allele = variant_alleles[allele].split(":")
                genotypes_allele = [genotypes_dict[int(genotype)]["good_name"] for genotype in genotypes_allele]
                
                genotypes_allele = sorted(genotypes_allele)
                sys.stdout.write("\t".join(genotypes_allele))
                sys.stdout.write("\n")
                
        elif output_fmt == "tabular":
            for genotype in genotypes_index_list:
                sys.stdout.write("\t"+genotypes_dict[genotype][var_id]+"\t")
            sys.stdout.write("\n")
        
        else:
            raise Exception("Unrecognized output format "+output_fmt+".")
    
    return

###############
###############

_print_parameters(options)

#### Parse queries file
####
query_list = []
if query_file == "":
    raise Exception("No queries especified.")
else:
    try:
        sys.stderr.write("Parsing queries file...\n")
        queries_file = open(query_file, 'r')
        for line in queries_file:
            line_data = line.strip().split()
            query_list.append(line_data[0])
        
    except Exception:
        sys.stderr.write("An error ocurred while parsing queries file.\n")
        traceback.print_exc()
    finally:
        queries_file.close()

#### Parse samples translation
####
samples_dict = {}
if samples_translation == "":
    pass
else:
    try:
        sys.stderr.write("Parsing samples translation file...\n")
        samples_file = open(samples_translation, 'r')
        for sample in samples_file:
            sample_data = sample.strip().split()
            samples_dict[sample_data[0]] = sample_data[1]
        
    except Exception:
        sys.stderr.write("An error ocurred while parsing samples translation file.\n")
        traceback.print_exc()
    finally:
        samples_file.close()

#### Parse headers file
####
genotypes_dict = {}
# If a headers file is specified, record names of genotypes
if vcf_header == "":
    pass
else:
    try:
        sys.stderr.write("Parsing VCF header file...\n")
        headers_file = open(vcf_header, 'r')
        for header in headers_file:
            line_data = header.strip().split()
            for j, genotype in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
                new_genotype = {"name":genotype, "index":j}
                if samples_translation == "":
                    new_genotype["good_name"] = genotype
                else:
                    new_genotype["good_name"] = samples_dict[genotype]
                genotypes_dict[j] = new_genotype
            
            break
        
    except Exception:
        sys.stderr.write("An error ocurred while parsing headers file.\n")
        traceback.print_exc()
    finally:
        headers_file.close()

#### Parse VCF file
####
variants_dict = {}
total_variants = 0
total_effects = 0
alleles_dict = {"h":""}
try:
    sys.stderr.write("Parsing VCF file...\n")
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        
        if line.startswith("##"): continue
        
        line_data = line.strip().split()
        
        # If a header is found, and if no header file was specified
        # record names of genotypes
        if line.startswith("#") and vcf_header == "":
            sys.stderr.write("Parsing VCF header...\n")
            for j, genotype in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
                new_genotype = {"name":genotype, "index":j}
                if samples_translation == "":
                    new_genotype["good_name"] = genotype
                else:
                    new_genotype["good_name"] = samples_dict[genotype]
                genotypes_dict[j] = new_genotype
                
            continue
        
        contig = line_data[VCF_CONTIG_COL]
        pos = line_data[VCF_POS_COL]
        ref = line_data[VCF_REF_COL]
        alt = line_data[VCF_ALT_COL]
        
        info = line_data[VCF_INFO_COL]
        # The split(";")[0] is to avoid fields beyond EFF= (like LOF=)
        effs_list = info.split(SNPEFF_FIELD)[1].split(";")[0].split(SNPEFF_EFF_SEP)
        
        used_variant = 0
        
        for eff in effs_list:
            eff_data = eff.split("(")
            eff_fields = eff_data[1].split("|")
            
            eff_type = eff_data[0]
            eff_gene = eff_fields[SNPEFF_GENE_POS]
            eff_isof = eff_fields[SNPEFF_ISOF_POS]
            
            if query_type == GENES and (eff_gene not in query_list): continue
            if query_type == ISOFS and (eff_isof not in query_list): continue
            
            total_effects += 1
            if used_variant == 0:
                used_variant = 1
                total_variants+=1
                
                variant_dict = {'var_id':total_variants, 'contig':contig, 'pos':pos, 'ref':ref, 'alt':alt, \
                                'alleles':{}, 'effs':[]}
                variants_dict[total_variants] = variant_dict
            
            eff_aa = eff_fields[SNPEFF_AA_POS]
            eff_dict = {'eff_gene':eff_gene, 'eff_type':eff_type, 'eff_isof':eff_isof, 'eff_aa':eff_aa}
            variant_dict["effs"].append(eff_dict)
        
        if used_variant == 0: continue
        
        for snpeff_data in info.split(SNPEFF_FIELD)[1].split(";")[1:]:
            sys.stderr.write("*snpEff other annotations: "+snpeff_data+"\n")
        
        alleles = variant_dict['alleles']
        for j, line_line in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
            genotype = line_line.split(":")[0]
            
            if len(genotype)==3:
                if genotype[0] == genotype[2]:
                    allele = genotype[0]
                elif (genotype[0] == "0" and genotype[2] == "1") or \
                    (genotype[0] == "1" and genotype[2] == "0"):
                    allele = "h"
                else:
                    allele = genotype
                
                if allele not in alleles_dict: alleles_dict[allele] = ""
            else:
                raise Exception("Genotype with more fields than expected: "+genotype+".")
            
            if allele in alleles:
                alleles[allele] += ":"+str(j)
            else:
                alleles[allele] = str(j)
            
            genotypes_dict[j][variant_dict["var_id"]] = allele
        
except Exception as e:
    sys.stderr.write("An error ocurred while parsing VCF file.\n")
    traceback.print_exc()
finally:
    vcf_file.close()

#### Output
####
print_variants(variants_dict, genotypes_dict, alleles_dict, output_format)

sys.stderr.write("Total variants read: "+str(total_variants)+"\n")
sys.stderr.write("Total effects read: "+str(total_effects)+"\n")
sys.stderr.write("Finished.\n")

## END