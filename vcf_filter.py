#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys, traceback
from optparse import OptionParser

from src.output import *

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
__usage = "usage: vcf_contigs_variants.py [VCF_FILE] [OPTIONS]"

optParser = OptionParser(__usage)

#optParser.add_option('-H', '--header', action='store', dest='vcf_header', type='string', \
#                     help='VCF header if VCF_FILE does not include header.')

optParser.add_option('-c', '--contigs_list', action='store', dest='contigs_list', type='string', \
                     help='')

optParser.add_option('-v', '--variants_list', action='store', dest='variants_list', type='string', \
                     help='')

optParser.add_option('-g', '--genes_list', action='store', dest='genes_list', type='string', \
                     help='')

optParser.add_option('-i', '--isof_list', action='store', dest='isof_list', type='string', \
                     help='')

optParser.add_option('-s', '--samples', action='store', dest='samples_filename', type='string', \
                     help='')

optParser.add_option('-t', '--samples_translation', action='store', dest='samples_translation', type='string', \
                     help='')

(options, arguments) = optParser.parse_args()

if not arguments or len(arguments)==0:
    optParser.exit(0, "You may wish to run '-help' option.\n")

vcf_filename = arguments[0] # This is mandatory

#if options.vcf_header: vcf_header = options.vcf_header
#else: vcf_header = ""

if options.contigs_list:
    query_file = options.contigs_list
else: query_file = ""

if options.variants_list:
    variants_file = options.variants_list
else: variants_file = ""

if options.genes_list:
    if options.isof_list:
        raise Exception("Only a list of genes or a list of isoforms a list of contigs should be provided, not both.")
    genes_file = options.genes_list
    query_type = GENES
elif options.isof_list:
    genes_file = options.isof_list
    query_type = ISOFS
else: genes_file = ""

if options.samples_filename: samples_filename = options.samples_filename
else: samples_filename = ""

if query_file == "" and variants_file == "" and samples_filename == "" and genes_file == "":
    raise Exception("Either a list of contigs or of variants is required.")

if options.samples_translation: samples_translation = options.samples_translation
else: samples_translation = ""

###############
###############

_print_parameters(options)

try:
    genotypes_dict = {}
    variants_dict = {}
    
    sys.stderr.write("Parsing filter files...\n")
    
    #### Parse queries file
    ####
    query_list = parse_queries_file(query_file)
    
    #### Variants to show
    ####
    variants_list = parse_queries_file(variants_file, keys=(1,2))
    
    #### Genes or isofoms
    ####
    genes_list = parse_queries_file(genes_file)
    
    #### Parse samples list
    ####
    samples_list = parse_samples_list(samples_filename)
    
    #### Parse samples translation
    ####
    samples_trans_dict = parse_samples_translation(samples_translation)
    
    #### Parse headers file
    ####
    #header_found = parse_vcf_header_file(vcf_header, genotypes_dict, \
    #                samples_filename, samples_list, samples_translation, samples_trans_dict)
    
    #### Parse VCF file
    ####
    total_variants = 0
    total_output = 0
    vcf_file = open(vcf_filename, 'r')
    sys.stderr.write("Parsing VCF file...\n")
    for line in vcf_file:
        
        if line.startswith("##"):
            sys.stdout.write(line.strip()+"\n")
            continue
        
        line_data = line.strip().split()
        
        # Header
        if line.startswith("#"):# and vcf_header == "":
            header_found = True
            if len(samples_list)>0:
                samples_fields = filter_header(line_data, samples_list, samples_translation, samples_trans_dict)
                print_filtered_samples(line_data, samples_fields)
            else:
                sys.stdout.write(line.strip()+"\n")
            continue
        
        if not header_found:
            raise Exception("No header found in VCF file.") # A header is needed
        
        total_variants += 1
        
        contig = line_data[VCF_CONTIG_COL]
        
        if query_file != "" and not contig in query_list: continue
        
        pos = line_data[VCF_POS_COL]
        
        if variants_file != "" and not [contig, pos] in variants_list: continue
        
        if genes_file != "":
            variant_genes_isof = get_variant_genes_isof(line_data, query_type)
            one_gene_isof_found = False
            for gene_isof in variant_genes_isof:
                if gene_isof in genes_list:
                    one_gene_isof_found = True
                    break
            if not one_gene_isof_found: continue
        
        if len(samples_list)>0:
            print_filtered_samples(line_data, samples_fields)
        else:
            sys.stdout.write(line.strip()+"\n")
        
        total_output += 1
    
    sys.stderr.write("Total variants read: "+str(total_variants)+"\n")
    sys.stderr.write("Total variants retained: "+str(total_output)+"\n")
    sys.stderr.write("Finished.\n")

except Exception as e:
    print e
    traceback.print_exc()
finally:
    vcf_file.close()

## END