#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

## TODO:
# Make use of cyvcf https://github.com/arq5x/cyvcf based on https://github.com/jamescasbon/PyVCF

import sys, traceback
from optparse import OptionParser

from src.output import *

def _print_parameters(options):
    options_dict = options.__dict__
    sys.stderr.write("Options:\n")
    for option in options_dict:
        if options_dict[option]:
            sys.stderr.write("\t"+option+"="+str(options_dict[option])+"\n")
    
    return

## Argument parsing
__usage = "usage: parse_contigs_variants.py [VCF_FILE] [OPTIONS]"

optParser = OptionParser(__usage)

optParser.add_option('-H', '--header', action='store', dest='vcf_header', type='string', \
                     help='VCF header if VCF_FILE does not include header.')

optParser.add_option('-c', '--contigs_list', action='store', dest='contigs_list', type='string', \
                     help='')

optParser.add_option('-v', '--variants_list', action='store', dest='variants_list', type='string', \
                     help='')

optParser.add_option('-s', '--samples', action='store', dest='samples_filename', type='string', \
                     help='')

optParser.add_option('-t', '--samples_translation', action='store', dest='samples_translation', type='string', \
                     help='')

optParser.add_option('--contigs_info', action='store', dest='contigs_info', type='string')

optParser.add_option('-f', '--output_format', action='store', dest='output_format', type='string', \
                     help='summary, detail or tabular')

optParser.add_option('--max_missing', action='store', dest='max_missing', type='float', help='')

optParser.add_option('--max_heteros', action='store', dest='max_heteros', type='float', help='')

optParser.add_option('--maf', action='store', dest='maf', type='float', help='')

optParser.add_option('--het_to_miss', action='store_true', dest='het_to_miss', help='')

optParser.add_option('-m', '--show_monomorph', action='store_true', dest='show_monomorph', help='')

optParser.add_option('-e', '--show_effects', action='store_true', dest='show_effects', help='')

optParser.add_option('-b', '--biallelic', action='store_true', dest='biallelic', help='')

optParser.add_option('-n', '--numeric', action='store_true', dest='numeric', help='')

optParser.add_option('-k', '--cluster_samples', action='store_true', dest='cluster_samples', help='')

(options, arguments) = optParser.parse_args()

if not arguments or len(arguments)==0:
    optParser.exit(0, "You may wish to run '-help' option.\n")

vcf_filename = arguments[0] # This is mandatory

if options.vcf_header: vcf_header = options.vcf_header
else: vcf_header = ""

if options.contigs_list:
    query_file = options.contigs_list
else: query_file = ""

if options.variants_list:
    variants_file = options.variants_list
else: variants_file = ""

#if query_file == "" and variants_file == "":
#    raise Exception("Either a list of contigs or of variants is required.")

if options.contigs_info: contigs_info = options.contigs_info
else: contigs_info = ""

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

if options.het_to_miss: het_to_miss = True
else: het_to_miss = False

if options.maf or options.maf == 0.0: maf = options.maf
else: maf = 0.0

if options.show_monomorph: show_monomorph = options.show_monomorph
else: show_monomorph = False

if options.show_effects: show_effects = options.show_effects
else: show_effects = False

if options.biallelic: biallelic = "bi"
else: biallelic = "mono"

if options.numeric: numeric = True
else: numeric = False

if options.cluster_samples: cluster_samples = True
else: cluster_samples = False

if cluster_samples and output_format != "tabular":
    sys.stderr.write("WARNING: samples won't be clustered. It is only compatible with output format \"tabular\".")

###############
###############

_print_parameters(options)

try:
    genotypes_dict = {}
    names_dict = {}
    variants_dict = {}
    
    #### Parse queries file
    ####
    if query_file != "":
        sys.stderr.write("Processing queries list...\n")
        query_list = parse_queries_file(query_file)
    
    #### Variants to show
    ####
    if variants_file != "":
        sys.stderr.write("Processing variants list...\n")
        variants_list = parse_queries_file(variants_file, keys=(1,2))
    
    #### Parse samples translation
    ####
    sys.stderr.write("Processing samples translation list...\n")
    samples_trans_dict = parse_samples_translation(samples_translation)
    
    #### Parse samples list
    ####
    sys.stderr.write("Parsing samples list...\n")
    samples_list = parse_samples_list(samples_filename)
    
    #### Parse headers file
    ####
    header_found = parse_vcf_header_file(vcf_header, genotypes_dict, names_dict, \
                                         samples_filename, samples_list, samples_translation, samples_trans_dict)
    
    #### Parse VCF file
    ####
    total_records = 0
    total_variants = 0
    total_output = 0
    sys.stderr.write("Parsing VCF file...\n")
    vcf_file = open(vcf_filename, 'r')
    for line in vcf_file:
        
        if line.startswith("##"): continue
        
        line_data = line.strip().split()
        
        # If a header is found, and if no header file was specified
        # record names of genotypes
        if line.startswith("#") and vcf_header == "":
            parse_header(line_data, genotypes_dict, names_dict, \
                         samples_filename, samples_list, samples_translation, samples_trans_dict)
            header_found = True
            continue
        
        if not header_found:
            raise Exception("No header found nor provided for VCF data.")
        
        total_records += 1
        
        contig = line_data[VCF_CONTIG_COL]
        if query_file != "" and not contig in query_list: continue
        
        pos = line_data[VCF_POS_COL]
        if variants_file != "" and not [contig, pos] in variants_list: continue
        
        total_variants+=1
        var_id = total_variants
        
        variant_dict = {'var_id':var_id, 'contig':contig, 'pos':pos, \
                        'ref':line_data[VCF_REF_COL], 'alt':line_data[VCF_ALT_COL], \
                        'alleles':{}, 'effs':{}, 'eff_x':[]}
        
        if show_effects: load_effects(line_data, variant_dict)
        
        alleles = parse_alleles(line_data, genotypes_dict, biallelic)
        
        variant_dict['alleles'] = alleles
        
        ok_variant = preprocess_variant(alleles, max_heteros, max_missing, show_monomorph, maf, \
                                        biallelic, het_to_miss)
        if ok_variant:
            variants_dict[var_id] = variant_dict
            total_output += 1
            for allele in alleles:
                for j in alleles[allele]:
                    genotypes_dict[j][var_id] = allele
    
    #### Output
    ####
    sys.stderr.write("Generating output...\n")
    print_variants_contigs(variants_dict, genotypes_dict, names_dict, samples_list, contigs_info, \
                           show_effects, output_format, \
                           biallelic, numeric, cluster_samples)
    
    sys.stderr.write("Total records read: "+str(total_records)+"\n")
    sys.stderr.write("Total variants parsed: "+str(total_variants)+"\n")
    sys.stderr.write("Total variants output: "+str(total_output)+"\n")
    sys.stderr.write("Finished.\n")

except Exception as e:
    print e
    traceback.print_exc()
finally:
    vcf_file.close()

## END