#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import sys, traceback

genes_dict = {}

vcf_filename = sys.argv[1]

sys.stderr.write("Input file: "+vcf_filename+"\n")

total_variants = 0
total_effects = 0
total_genes = 0
try:
    vcf_file = open(vcf_filename)
    for line in vcf_file:
            #code
            if line.startswith("#"): continue
            
            total_variants += 1
            line_data = line.strip().split()
            #print line_data
            info = line_data[7]
            
            # The split(";")[0] is to avoid fields beyond EFF= (like LOF=)
            effs_list = info.split("EFF=")[1].split(";")[0].split(",")
            
            for eff in effs_list:
                total_effects += 1
                eff_data = eff.split("(")
                eff_type = eff_data[0]
                eff_fields = eff_data[1].split("|")
                
                eff_gene = eff_fields[5]                
                eff_isof = eff_fields[8]
                
                if eff_gene == "": eff_gene = "-"
                if eff_isof == "": eff_isof = "-"
                
                if eff_gene in genes_dict:
                    gene_dict = genes_dict[eff_gene]
                    if eff_isof in gene_dict:
                        isof_dict = gene_dict[eff_isof]
                        if eff_type in isof_dict:
                            isof_dict[eff_type] += 1
                        else:
                            isof_dict = {eff_type:1}
                    else:
                        isof_dict = {eff_type:1}
                        gene_dict[eff_isof] = isof_dict
                else:
                    total_genes += 1
                    isof_dict = {eff_type:1}
                    gene_dict = {eff_isof:isof_dict}
                    genes_dict[eff_gene] = gene_dict
                    
except Exception as e:
    sys.stderr.write("An error ocurred while parsing VCF file.\n")
    #print line
    #traceback.print_exc()
    sys.exit(-1)
finally:
    vcf_file.close()

for gene in genes_dict:
    gene_dict = genes_dict[gene]
    for isof in gene_dict:
        isof_dict = gene_dict[isof]
        for eff_type in isof_dict:
            eff_type_count = isof_dict[eff_type]
            outlist = [gene,isof,eff_type,str(eff_type_count)]
            sys.stdout.write("\t".join(outlist)+"\n")
    
sys.stderr.write("Total variants read: "+str(total_variants)+"\n")
sys.stderr.write("Total effects read: "+str(total_effects)+"\n")
sys.stderr.write("Total genes found: "+str(total_genes)+"\n")

sys.stderr.write("Finished.\n")

## END