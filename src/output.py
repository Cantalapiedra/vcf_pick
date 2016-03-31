#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys
from util.vcf_util import *
from util.util import *

GENES = 0
ISOFS = 1

def f_get_samples_rows(outputs_list, samples_list, genotypes_dict, biallelic):
    samples_rows = []
    
    header_row = []
    for sample in samples_list:
        header_row.append(sample)
    
    samples_rows.append(header_row)
    
    for output_line in outputs_list:
        var_id = output_line[0]
        var_row = []
        for sample in samples_list:
            for genotype in genotypes_dict:
                if genotypes_dict[genotype]["good_name"] == sample:
                    genotype_var_allele = get_numeric_allele(genotypes_dict[genotype][var_id], biallelic)
                    var_row.append(genotype_var_allele)
        samples_rows.append(var_row)
    
    return samples_rows

#

def print_rows(outputs_list, samples_list, variants_dict, genotypes_dict, output_fmt, biallelic, numeric):
    # For each variant/gene position to output
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
                for genotype in genotypes_dict:
                    if genotypes_dict[genotype]["good_name"] == sample:
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

def print_variants_contigs(variants_dict, genotypes_dict, samples_list, show_effects, output_fmt = "tabular", \
                           biallelic = "mono", numeric = False, cluster_samples = False):
    
    # Prepare output lines
    outputs_list = []
    for variant_id in variants_dict:
        variant = variants_dict[variant_id]
        
        output_list = []
        
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
            
            output_list = [variant_id, variant["contig"], variant["pos"], variant["ref"], variant["alt"], \
                           ",".join(eff_type), ",".join(eff_isof), ",".join(eff_aa), ",".join(variant["eff_x"])]
        else:
            output_list = [variant_id, variant["contig"], variant["pos"], variant["ref"], variant["alt"]]
        
        outputs_list.append(output_list)
    
    # Sort the list of variants to output
    # by isoform, contig and position
    outputs_list = sorted(outputs_list, key=lambda x: (x[1], int(x[2])))  
    
    ### Print output
    genotypes_index_list = sorted(genotypes_dict, key=lambda x: genotypes_dict[x]["good_name"])
    
    # Header
    if show_effects:
        sys.stdout.write("#\tcontig\tpos\tref\talt\teffect\tisof\tchange\tother")
    else:
        sys.stdout.write("#\tcontig\tpos\tref\talt")
    
    if output_fmt == "tabular":
        if cluster_samples:
            samples_rows = f_get_samples_rows(outputs_list, samples_list, genotypes_dict, biallelic)
            samples_list = f_cluster_samples(samples_rows, biallelic)
        # else: samples_list = samples_list
        for sample in samples_list:
            sys.stdout.write("\t"+sample)
        
    sys.stdout.write("\n") # END of header row
    
    # Rows
    print_rows(outputs_list, samples_list, variants_dict, genotypes_dict, output_fmt, biallelic, numeric)
    
    return

##

def print_variants_genes(variants_dict, genotypes_dict, samples_list, query_type, output_fmt = "tabular", \
                         biallelic = "mono", numeric = False, cluster_samples = False):
    
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
        if cluster_samples:
            samples_rows = f_get_samples_rows(outputs_list, samples_list, genotypes_dict, biallelic)
            samples_list = f_cluster_samples(samples_rows, biallelic)
        # else: samples_list = samples_list
        for sample in samples_list:
            sys.stdout.write("\t"+sample)
            
    sys.stdout.write("\n")
    
    # Rows
    print_rows(outputs_list, samples_list, variants_dict, genotypes_dict, output_fmt, biallelic, numeric)
    
    return

## END