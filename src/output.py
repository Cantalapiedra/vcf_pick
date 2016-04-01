#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys
from util.vcf_util import *
from util.util import *

def f_get_samples_rows(outputs_list, samples_list, genotypes_dict, names_dict, biallelic):
    samples_rows = []
    
    sys.stderr.write("Preparing data of samples for clustering...\n")
    header_row = []
    header_row = samples_list
    #for sample in samples_list:
    #    header_row.append(sample)
    
    samples_rows.append(header_row)
    
    for output_line in outputs_list:
        var_id = output_line[0]
        var_row = []
        for sample in samples_list:
            genotype = names_dict[sample]
            #for genotype in genotypes_dict:
            #    if genotypes_dict[genotype]["good_name"] == sample:
            genotype_var_allele = get_numeric_allele(genotypes_dict[genotype][var_id], biallelic)
            var_row.append(genotype_var_allele)
        samples_rows.append(var_row)
    
    sys.stderr.write("Ready for clustering.\n")
    
    return samples_rows

#

def print_rows(outputs_list, samples_list, variants_dict, genotypes_dict, names_dict, output_fmt, biallelic, numeric):
    # For each variant/gene position to output
    sys.stderr.write("Printing genotypes...\n")
    
    for output_line in outputs_list:
        var_id = output_line[0]
        #isof_id = output_line[1]
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
                genotype = names_dict[sample]
                #for genotype in genotypes_dict:
                #    if genotypes_dict[genotype]["good_name"] == sample:
                if biallelic == "bi" or not numeric:
                    sys.stdout.write("\t"+genotypes_dict[genotype][var_id])
                else: # biallelic = "mono" and numeric:
                    numeric_allele = get_numeric_allele(genotypes_dict[genotype][var_id], biallelic)
                    sys.stdout.write("\t"+str(numeric_allele))
            sys.stdout.write("\n")
        
        else:
            raise Exception("Unrecognized output format "+output_fmt+".")
    
    return

#

def load_contigs_info(contigs_info):
    contigs_info_dict = {}
    
    ## TODO: check all rows have the same number of fields or raise an
    ## Exception("The info file should have the same number of fields/columns in all the rows.")
    if contigs_info != "":
        for i, line in enumerate(open(contigs_info, 'r')):
            line_data = line.strip().split("\t")
            if i == 0:
                contigs_info_dict["header"] = line_data[1:]
                contigs_info_dict["void"] = len(line_data[1:])*["-"]
                contigs_info_dict["fields"] = len(line_data[1:])
            else:
                if line_data[0] not in contigs_info_dict:
                    contigs_info_dict[line_data[0]] = line_data[1:]
                else:
                    sys.stderr.write("WARNING: duplicated info for "+str(line_data[0])+"\n")
    
    return contigs_info_dict

##

def get_contigs_header(show_effects, contigs_info_dict):
    header_list = ["#", "contig", "pos"]
    if len(contigs_info_dict)>0:
        header_list = header_list + contigs_info_dict["header"]
    header_list = header_list + ["ref", "alt"]
    
    if show_effects:
        header_list = header_list + ["effect", "isof", "change", "other"]
    
    return header_list

def get_contigs_variant(show_effects, variant, contigs_info_dict):
    contig = variant["contig"]
    output_list = [variant["var_id"], contig, variant["pos"]]
    if len(contigs_info_dict)>0:
        if contig in contigs_info_dict:
            output_list = output_list + contigs_info_dict[contig]
        else:
            output_list = output_list + contigs_info_dict["void"]
        
    output_list = output_list + [variant["ref"], variant["alt"]]
    
    eff_type = []
    eff_isof = []
    eff_aa = []
    
    if show_effects:
        if len(variant["effs"]) == 0:
            eff_type = ["-"]
            eff_isof = ["-"]
            eff_aa = ["-"]
        else:
            for eff in sorted(variant["effs"]):
                variant_eff = variant["effs"][eff]
                eff_type.append(eff)
                eff_isof.append("("+",".join(variant_eff["eff_isof"])+")")
                eff_aa.append("("+",".join(variant_eff["eff_aa"])+")")
        
        output_list = output_list + [",".join(eff_type), ",".join(eff_isof), ",".join(eff_aa), ",".join(variant["eff_x"])]
    
    return output_list

def print_variants_contigs(variants_dict, genotypes_dict, names_dict, samples_list, contigs_info = "", \
                           show_effects = False, output_fmt = "tabular", \
                           biallelic = "mono", numeric = False, cluster_samples = False):
    
    contigs_info_dict = load_contigs_info(contigs_info)
    
    # Prepare output lines
    outputs_list = []
    for variant_id in variants_dict:
        variant = variants_dict[variant_id]
        output_list = get_contigs_variant(show_effects, variant, contigs_info_dict)
        outputs_list.append(output_list)
    
    # Sort the list of variants to output
    # by isoform, contig and position
    outputs_list = sorted(outputs_list, key=lambda x: (x[1], int(x[2])))  
    
    # Header
    header_list = get_contigs_header(show_effects, contigs_info_dict)
    sys.stdout.write("\t".join(header_list))
    
    if output_fmt == "tabular":
        if cluster_samples:
            samples_rows = f_get_samples_rows(outputs_list, samples_list, genotypes_dict, names_dict, biallelic)
            samples_list = f_cluster_samples(samples_rows, biallelic)
        # else: samples_list = samples_list
        for sample in samples_list:
            sys.stdout.write("\t"+sample)
        
    sys.stdout.write("\n") # END of header row
    
    # Rows
    print_rows(outputs_list, samples_list, variants_dict, genotypes_dict, names_dict, output_fmt, biallelic, numeric)
    
    return

##

def get_genes_header(contigs_info_dict, genes_info_dict):
    header_list = ["#", "isof"]
    if len(genes_info_dict)>0:
        header_list = header_list + genes_info_dict["header"]
    
    header_list = header_list + ["effect", "eff_other", "ref", "alt", "change", "contig", "pos"]
    
    if len(contigs_info_dict)>0:
        header_list = header_list + contigs_info_dict["header"]
    
    return header_list

def get_genes_effect(variant, variant_eff_isof, variant_eff_aa, eff, contigs_info_dict, genes_info_dict):
    output_list = []
    
    contig = variant["contig"]
    
    output_list = [variant["var_id"], variant_eff_isof]
    
    if len(genes_info_dict)>0:
        if variant_eff_isof in genes_info_dict:
            output_list = output_list + genes_info_dict[variant_eff_isof]
        else:
            output_list = output_list + genes_info_dict["void"]
    
    output_list = output_list + [eff, ",".join(variant["eff_x"]), \
                    variant["ref"], variant["alt"], variant_eff_aa, \
                    contig, variant["pos"]]
    
    if len(contigs_info_dict)>0:
        if contig in contigs_info_dict:
            output_list = output_list + contigs_info_dict[contig]
        else:
            output_list = output_list + contigs_info_dict["void"]
    
    return output_list

def print_variants_genes(variants_dict, genotypes_dict, names_dict, samples_list, contigs_info, genes_info, \
                         output_fmt = "tabular", \
                         biallelic = "mono", numeric = False, cluster_samples = False):
    
    contigs_info_dict = load_contigs_info(contigs_info)
    genes_info_dict = load_contigs_info(genes_info) # Using the same function as for contigs
    
    # Prepare output lines
    outputs_list = []
    for variant_id in variants_dict:
        variant = variants_dict[variant_id]
        
        for eff in variant["effs"]:
            variant_effs = variant["effs"][eff]
            for i, variant_eff_isof in enumerate(variant_effs["eff_isof"]):
                variant_eff_aa = variant_effs["eff_aa"][i]
                output_list = get_genes_effect(variant, variant_eff_isof, variant_eff_aa, eff, \
                                               contigs_info_dict, genes_info_dict)
                outputs_list.append(output_list)    
        
    # Sort the list of variants to output
    # by isoform, contig and position
    if genes_info == "":
        outputs_list = sorted(outputs_list, key=lambda x: (x[1], x[7], int(x[8])))
    else:
        add_fields = genes_info_dict["fields"]
        
        outputs_list = sorted(outputs_list, key=lambda x: (x[1], x[7+add_fields], int(x[8+add_fields])))
    
    # Header
    header_list = get_genes_header(contigs_info_dict, genes_info_dict)
    sys.stdout.write("\t".join(header_list))
    
    if output_fmt == "tabular":
        if cluster_samples:
            samples_rows = f_get_samples_rows(outputs_list, samples_list, genotypes_dict, names_dict, biallelic)
            samples_list = f_cluster_samples(samples_rows, biallelic)
        # else: samples_list = samples_list
        for sample in samples_list:
            sys.stdout.write("\t"+sample)
    
    sys.stdout.write("\n")
    
    # Rows
    print_rows(outputs_list, samples_list, variants_dict, genotypes_dict, names_dict, output_fmt, biallelic, numeric)
    
    return

##

def print_filtered_samples(line_data, samples_fields):
    sys.stdout.write("\t".join(line_data[:VCF_SAMPLES_INI_COL]+\
                                [line_data[i] for i in samples_fields])+"\n")


## END