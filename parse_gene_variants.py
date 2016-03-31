#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys, traceback
from optparse import OptionParser

from src.util.vcf_util import *
from src.util.util import *

GENES = 0
ISOFS = 1

def _print_parameters(options):
    options_dict = options.__dict__
    sys.stderr.write("Options:\n")
    for option in options_dict:
        if options_dict[option]:
            sys.stderr.write("\t"+option+"="+str(options_dict[option])+"\n")
    
    return

## Argument parsing
__usage = "usage: parse_genes_variants.py [VCF_FILE] [OPTIONS]"

optParser = OptionParser(__usage)

optParser.add_option('-H', '--header', action='store', dest='vcf_header', type='string', \
                     help='VCF header if VCF_FILE does not include header.')

optParser.add_option('-g', '--genes_list', action='store', dest='genes_list', type='string', \
                     help='')

optParser.add_option('-i', '--isof_list', action='store', dest='isof_list', type='string', \
                     help='')

optParser.add_option('-v', '--variants_list', action='store', dest='variants_list', type='string', \
                     help='')

optParser.add_option('-t', '--samples_translation', action='store', dest='samples_translation', type='string', \
                     help='')

optParser.add_option('-f', '--output_format', action='store', dest='output_format', type='string', \
                     help='summary or tabular')

optParser.add_option('-s', '--samples', action='store', dest='samples_filename', type='string', \
                     help='')

optParser.add_option('--max_missing', action='store', dest='max_missing', type='float', help='')

optParser.add_option('--max_heteros', action='store', dest='max_heteros', type='float', help='')

optParser.add_option('--maf', action='store', dest='maf', type='float', help='')

optParser.add_option('-m', '--show_monomorph', action='store_true', dest='show_monomorph', help='')

optParser.add_option('-b', '--biallelic', action='store_true', dest='biallelic', help='')

optParser.add_option('-n', '--numeric', action='store_true', dest='numeric', help='')

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

if options.variants_list:
    variants_file = options.variants_list
else: variants_file = ""

if options.samples_translation: samples_translation = options.samples_translation
else: samples_translation = ""

if options.output_format: output_format = options.output_format
else: output_format = "tabular"

if options.samples_filename: samples_filename = options.samples_filename
else: samples_filename = ""

if options.max_missing or options.max_missing == 0.0: max_missing = options.max_missing
else: max_missing = 1.0

if options.max_heteros or options.max_heteros == 0.0: max_heteros = options.max_heteros
else: max_heteros = 1.0

if options.maf or options.maf == 0.0: maf = options.maf
else: maf = 0.0

if options.show_monomorph: show_monomorph = options.show_monomorph
else: show_monomorph = False

if options.biallelic: biallelic = "bi"
else: biallelic = "mono"

if options.numeric: numeric = True
else: numeric = False

###############
###############

def print_variants_genes(variants_dict, genotypes_dict, samples_list, output_fmt = "tabular", \
                         biallelic = "mono", numeric = False):
    
    # Prepare output lines
    outputs_list = []
    for variant_id in variants_dict:
        variant = variants_dict[variant_id]
        
        output_list = []
        for eff in variant["effs"]:
            variant_effs = variant["effs"][eff]
            for i, variant_eff_isof in enumerate(variant_effs["eff_isof"]):
                variant_eff_aa = variant_effs["eff_aa"][i]
                if query_type == GENES:
                    output_list = [variant_id, variant_eff_isof, eff, ",".join(variant["eff_x"]), \
                                    variant["ref"], variant["alt"], variant_eff_aa, \
                                    variant["contig"], variant["pos"]]
                elif query_type == ISOFS:
                    output_list = [variant_id, variant_eff_isof, eff, ",".join(variant["eff_x"]), \
                                    variant["ref"], variant["alt"], variant_eff_aa, \
                                    variant["contig"], variant["pos"]]
                else:
                    sys.stderr.write("Unrecognized query type "+query_type+"\n")
                
                outputs_list.append(output_list)    
        
    # Sort the list of variants to output
    # by isoform, contig and position
    outputs_list = sorted(outputs_list, key=lambda x: (x[1], x[7], int(x[8])))  
    
    ### Print output
    genotypes_index_list = sorted(genotypes_dict, key=lambda x: genotypes_dict[x]["good_name"])
    
    # Header
    sys.stdout.write("#\tisof\teffect\teff_other\tref\talt\tchange\tcontig\tpos")
    if output_fmt == "tabular":
        for sample in samples_list:
            sys.stdout.write("\t"+sample)
            
    sys.stdout.write("\n")
    
    # Genotypes
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
            
            for allele in variant_alleles:
                sys.stdout.write(">>\t")
                sys.stdout.write(allele+"\t")
                genotypes_allele = variant_alleles[allele]
                genotypes_allele = [genotypes_dict[int(genotype)]["good_name"] for genotype in genotypes_allele]
                
                genotypes_allele = sorted(genotypes_allele)
                sys.stdout.write("\t".join(genotypes_allele))
                sys.stdout.write("\n")
                
        elif output_fmt == "tabular":
            for sample in samples_list:
                for genotype in genotypes_dict:
                    if genotypes_dict[genotype]["good_name"] == sample:
                        if biallelic == "bi" or not numeric:
                            sys.stdout.write("\t"+genotypes_dict[genotype][var_id])
                        else: # biallelic = "mono" and numeric:
                            genotype_var_allele = genotypes_dict[genotype][var_id]
                            if genotype_var_allele == MISS:
                                genotype_var_allele = -1
                            elif genotype_var_allele == HETERO:
                                genotype_var_allele = 0.5
                            else:
                                genotype_var_allele = int(genotype_var_allele) * 1.0
                            sys.stdout.write("\t"+str(genotype_var_allele))
            sys.stdout.write("\n")
        
        else:
            raise Exception("Unrecognized output format "+output_fmt+".")
    
    return

###############
###############

_print_parameters(options)

#### Parse queries file
####
query_list = parse_queries_file(query_file)

#### Variants to show
####
if variants_file != "":
    variants_list = parse_queries_file(variants_file, keys=(1,2))

#### Parse samples translation
####
samples_trans_dict = parse_samples_translation(samples_translation)

#### Parse samples list
####
samples_list = parse_samples_list(samples_filename)

#### Parse headers file
####
genotypes_dict = {}
header_found = parse_vcf_header_file(vcf_header, genotypes_dict, \
                                     samples_filename, samples_list, samples_translation, samples_trans_dict)

#### Parse VCF file
####
variants_dict = {}
total_variants = 0
snpeff_field_found = False
try:
    sys.stderr.write("Parsing VCF file...\n")
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        
        if line.startswith("##"): continue
        
        line_data = line.strip().split()
        
        # If a header is found, and if no header file was specified
        # record names of genotypes
        if line.startswith("#") and vcf_header == "":
            parse_header(line_data, genotypes_dict, \
                         samples_filename, samples_list, samples_translation, samples_trans_dict)
            header_found = True
            continue
        
        if not header_found:
            raise Exception("No header found nor provided for VCF data.")
        
        contig = line_data[VCF_CONTIG_COL]
        pos = line_data[VCF_POS_COL]
        
        if variants_file != "" and not [contig, pos] in variants_list: continue
        
        variant_dict = {'var_id':-1, 'contig':contig, 'pos':pos, \
                        'ref':line_data[VCF_REF_COL], 'alt':line_data[VCF_ALT_COL], \
                        'alleles':{}, 'effs':{}, 'eff_x':[]}
        
        info = line_data[VCF_INFO_COL]
        # The split(";")[0] is to avoid fields beyond EFF= (like LOF=)
        if (not snpeff_field_found) and SNPEFF_FIELD in info:
            snpeff_field_found = True
        
        if not SNPEFF_FIELD in info:
            sys.stderr.write("WARNING: VCF record withouth snpEff field ("+SNPEFF_FIELD+") found: "+line)
            continue
        
        effs_list = info.split(SNPEFF_FIELD)[1].split(";")[0].split(SNPEFF_EFF_SEP)
        
        if not parse_effects_genes(effs_list, query_type, query_list, total_variants, variant_dict):
            continue
        
        total_variants+=1
        var_id = total_variants
        variant_dict["var_id"] = var_id
        
        for snpeff_data in info.split(SNPEFF_FIELD)[1].split(";")[1:]:
            variant_dict["eff_x"].append(snpeff_data[:3])
        if len(variant_dict["eff_x"]) == 0: variant_dict["eff_x"] = ["-"]
        
        alleles = parse_alleles(line_data, genotypes_dict, biallelic)
        
        variant_dict['alleles'] = alleles
        
        if preprocess_variant(alleles, max_heteros, max_missing, show_monomorph, maf, biallelic):
            variants_dict[var_id] = variant_dict
            for allele in alleles:
                for j in alleles[allele]:
                    genotypes_dict[j][var_id] = allele

    if not snpeff_field_found:
        raise Exception("snpEff field ("+SNPEFF_FIELD+") \
                        should be present to parse genes info.\n")
    
    #### Output
    ####
    print_variants_genes(variants_dict, genotypes_dict, samples_list, output_format, biallelic, numeric)
    
    sys.stderr.write("Total variants read: "+str(total_variants)+"\n")
    sys.stderr.write("Finished.\n")

except Exception as e:
    sys.stderr.write("An error ocurred while parsing VCF file.\n")
    traceback.print_exc()
finally:
    vcf_file.close()

## END