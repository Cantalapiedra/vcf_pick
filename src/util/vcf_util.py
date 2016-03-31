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
VCF_GENOTYPE_SEP = ":"
VCF_GENOTYPE_FIELD = 0
VCF_ALLELES_SEP = "/"
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

def get_variant_genes_isof(line_data, query_type):
    GENES = 0
    ISOFS = 1
    genes_isof_list = []
    
    info = line_data[VCF_INFO_COL]
    
    if SNPEFF_FIELD in info:
        effs_list = info.split(SNPEFF_FIELD)[1].split(";")[0].split(SNPEFF_EFF_SEP)
        for eff in effs_list:
            eff_data = eff.split("(")
            eff_fields = eff_data[1].split("|")
            
            if query_type == GENES: genes_isof_list.append(eff_fields[SNPEFF_GENE_POS])
            if query_type == ISOFS: genes_isof_list.append(eff_fields[SNPEFF_ISOF_POS])
        
    else:
        raise Exception("Could not find snpEff fields, which is required to filter by genes or isoforms.")
    
    return genes_isof_list

##

def get_numeric_allele(allele, biallelic):
    if biallelic == "bi":
        raise Exception("Biallelic alleles not yet allowed to yield numeric alleles.")
    # else: biallelic == "mono":
    if allele == MISS:
        allele = -1.0
    elif allele == HETERO:
        allele = 0.5
    else:
        allele = int(allele) * 1.0
    
    return allele

##

def load_effects_genes(line_data, variant_dict, query_type, query_list):
    ok_variant = False
    
    info = line_data[VCF_INFO_COL]
    
    if SNPEFF_FIELD in info:
        effs_list = info.split(SNPEFF_FIELD)[1].split(";")[0].split(SNPEFF_EFF_SEP)
        
        if parse_effects_genes(effs_list, query_type, query_list, variant_dict):
            ok_variant = True
            
            # The split(";")[0] is to avoid fields beyond EFF= (like LOF=)
            for snpeff_data in info.split(SNPEFF_FIELD)[1].split(";")[1:]:
                variant_dict["eff_x"].append(snpeff_data[:3])
            if len(variant_dict["eff_x"]) == 0: variant_dict["eff_x"] = ["-"]   
    else:
        sys.stderr.write("WARNING: VCF record withouth snpEff field ("+SNPEFF_FIELD+") found: "+line)
    
    return ok_variant


def load_effects(line_data, variant_dict):
    info = line_data[VCF_INFO_COL]
    
    if SNPEFF_FIELD in info:
        effs_list = info.split(SNPEFF_FIELD)[1].split(";")[0].split(SNPEFF_EFF_SEP)
        parse_effects_contig(effs_list, variant_dict)
        
        for snpeff_data in info.split(SNPEFF_FIELD)[1].split(";")[1:]:
            variant_dict["eff_x"].append(snpeff_data[:3])
        if len(variant_dict["eff_x"]) == 0: variant_dict["eff_x"] = ["-"]
    
    return

def parse_effects_contig(effs_list, variant_dict):
    effects_dict = variant_dict["effs"]
    
    for eff in effs_list:
        eff_data = eff.split("(")
        eff_fields = eff_data[1].split("|")
        
        eff_type = eff_data[0]
        eff_gene = eff_fields[SNPEFF_GENE_POS]
        eff_isof = eff_fields[SNPEFF_ISOF_POS]
        eff_aa = eff_fields[SNPEFF_AA_POS]
        
        if eff_type in effects_dict:
            eff_dict = effects_dict[eff_type]
            eff_dict["eff_gene"].append(eff_gene)
            eff_dict["eff_isof"].append(eff_isof)
            eff_dict["eff_aa"].append(eff_aa)
        else:
            eff_dict = {'eff_gene':[eff_gene], 'eff_type':eff_type, 'eff_isof':[eff_isof], 'eff_aa':[eff_aa]}
            effects_dict[eff_type] = eff_dict
    
    return

def parse_effects_genes(effs_list, query_type, query_list, variant_dict):
    GENES = 0
    ISOFS = 1
    
    retValue = False
    effects_dict = variant_dict["effs"]
    
    for eff in effs_list:
        eff_data = eff.split("(")
        eff_fields = eff_data[1].split("|")
        
        eff_type = eff_data[0]
        eff_gene = eff_fields[SNPEFF_GENE_POS]
        eff_isof = eff_fields[SNPEFF_ISOF_POS]
        eff_aa = eff_fields[SNPEFF_AA_POS]
        
        if query_type == GENES and (eff_gene not in query_list): continue
        if query_type == ISOFS and (eff_isof not in query_list): continue
        
        if not retValue: retValue = True
        
        if eff_type in effects_dict:
            eff_dict = effects_dict[eff_type]
            eff_dict["eff_gene"].append(eff_gene)
            eff_dict["eff_isof"].append(eff_isof)
            eff_dict["eff_aa"].append(eff_aa)
        else:
            eff_dict = {'eff_gene':[eff_gene], 'eff_type':eff_type, 'eff_isof':[eff_isof], 'eff_aa':[eff_aa]}
            effects_dict[eff_type] = eff_dict
    
    return retValue

##

def filter_header(line_data, samples_list = [], samples_translation = "", samples_trans_dict = {}):
    samples_fields = []
    
    for j, genotype in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
        if samples_translation == "":
            if genotype in samples_list:
                samples_fields.append(j+VCF_SAMPLES_INI_COL)
        else:
            if samples_trans_dict[genotype] in samples_list:
                samples_fields.append(j+VCF_SAMPLES_INI_COL)
    
    return samples_fields

##

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

def parse_alleles(line_data, genotypes_dict, alleles_output = "mono"):
    alleles = {}
    
    for j, line_line in enumerate(line_data[VCF_SAMPLES_INI_COL:]):
        
        if j not in genotypes_dict: continue
        
        genotype = line_line.split(VCF_GENOTYPE_SEP)[VCF_GENOTYPE_FIELD]
        
        if alleles_output == "mono":
            allele = parse_monoallelic(genotype)
        elif alleles_output == "bi":
            allele = parse_biallelic(genotype)
        else:
            raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
        
        if allele in alleles:
            alleles[allele].append(j)
        else:
            alleles[allele] = [j]
    
    return alleles

def parse_monoallelic(genotype):
    allele = ""
    
    alleles = genotype.split(VCF_ALLELES_SEP)
    if len(alleles)==2:
        if alleles[0] == alleles[1]:
            allele = alleles[0]
        else: # alleles[0] != alleles[1]
            allele = HETERO
    else:
        raise Exception("Genotype with more fields than expected: "+genotype+".")
    
    return allele

def parse_biallelic(genotype):
    return genotype

##

 # Heteros --> missing; missing --> rm; monomorph --> rm
def preprocess_variant(alleles, max_heteros, max_missing, show_monomorph, maf = 0.3, alleles_output = "mono"):
    retValue = True
    
    genotypes_keys = set(alleles)
    
    if len(genotypes_keys) <= 1 and not show_monomorph:
        retValue = False
    else:
        num_alleles = get_num_alleles(alleles, alleles_output)
        num_genotypes = get_num_genotypes(num_alleles, alleles_output)
        
        # Heteros --> missing
        percent_heteros = preprocess_heteros(alleles, genotypes_keys, num_genotypes, max_heteros, alleles_output)
        
        genotypes_keys = set(alleles)
        
        num_alleles = get_num_alleles(alleles, alleles_output)
        num_genotypes = get_num_genotypes(num_alleles, alleles_output)
        
        if len(genotypes_keys) <= 1 and not show_monomorph:
            retValue = False
        else:
            percent_miss = calculate_percent_miss(genotypes_keys, alleles, num_genotypes, alleles_output)
            
            if percent_miss > max_missing:
                retValue = False
            else:
                ## MAF
                alleles_keys = get_alleles_keys(genotypes_keys, alleles_output)
                alleles_maf = calculate_maf(alleles_keys, alleles, num_alleles, alleles_output)
                if alleles_maf < maf:
                    retValue = False
    
    #print "\t"+str(alleles)+"\t"+str(retValue)
    
    return retValue

##

def get_num_alleles(alleles, alleles_output):
    num_alleles = 0
    if alleles_output == "mono":
        num_alleles = sum([len(allele_genotypes) for allele_genotypes in alleles.values()])+1
        
    elif alleles_output == "bi":
        num_alleles = sum([len(allele_genotypes) for allele_genotypes in alleles.values()])+1
        num_alleles = num_alleles * 2
    else:
        raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
    
    return num_alleles

##

def get_num_genotypes(num_alleles, alleles_output):
    num_genotypes = 0
    
    if alleles_output == "mono":
        num_genotypes = num_alleles
    elif alleles_output == "bi":
        num_genotypes = num_alleles / 2
    else:
        raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
    
    return num_genotypes

##

def get_alleles_keys(genotypes_keys, alleles_output):
    alleles_keys = []
    
    if alleles_output == "mono":
        alleles_keys = genotypes_keys
    elif alleles_output == "bi":
        for genotype in genotypes_keys:
            alleles = genotype.split(VCF_ALLELES_SEP)
            for allele in alleles:
                if allele not in set(alleles_keys):
                    alleles_keys.append(allele)
    else:
        raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
    
    return alleles_keys

##

def preprocess_heteros(alleles, genotypes_keys, num_genotypes, max_heteros, alleles_output):
    percent_heteros = 0
    
    if alleles_output == "mono":
        percent_heteros = preprocess_heteros_mono(alleles, genotypes_keys, num_genotypes, max_heteros)
    elif alleles_output == "bi":
        percent_heteros = preprocess_heteros_bi(alleles, genotypes_keys, num_genotypes, max_heteros)
    else:
        raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
    
    return percent_heteros

def preprocess_heteros_mono(alleles, genotypes_keys, num_genotypes, max_heteros):
    percent_heteros = 0
    
    ## Heteros --> missing
    if HETERO in genotypes_keys:
        percent_heteros = len(alleles[HETERO]) * 1.0 / num_genotypes
    
    if percent_heteros > max_heteros:
        if MISS not in alleles:
            alleles[MISS] = []
        alleles[MISS] = alleles[MISS] + alleles[HETERO]
        del alleles[HETERO]
    
    return percent_heteros

def preprocess_heteros_bi(alleles, genotypes_keys, num_genotypes, max_heteros):
    percent_heteros = 0
    
    ## Heteros --> missing
    num_heteros = 0
    genotype_heteros = []
    for genotype in genotypes_keys:
        genotype_alleles = genotype.split(VCF_ALLELES_SEP) 
        if genotype_alleles[0] != genotype_alleles[1]:
            num_heteros += len(alleles[genotype])
            genotype_heteros.append(genotype)
    
    percent_heteros = num_heteros * 1.0 / num_genotypes
    
    if percent_heteros > max_heteros:
        missing = MISS+VCF_ALLELES_SEP+MISS
        if missing not in alleles:
            alleles[missing] = []
        for genotype_hetero in genotype_heteros:
            alleles[missing]= alleles[missing] + alleles[genotype_hetero]
            del alleles[genotype_hetero]
    
    return percent_heteros

##

def calculate_percent_miss(genotypes_keys, alleles, num_genotypes, alleles_output):
    percent_miss = 0
    
    if alleles_output == "mono":
        percent_miss = calculate_percent_miss_mono(genotypes_keys, alleles, num_genotypes)
    elif alleles_output == "bi":
        percent_miss = calculate_percent_miss_bi(genotypes_keys, alleles, num_genotypes)
    else:
        raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
    
    return percent_miss

def calculate_percent_miss_mono(genotypes_keys, alleles, num_genotypes):
    percent_miss = 0
    
    if MISS in genotypes_keys:
        percent_miss = len(alleles[MISS]) * 1.0 / num_genotypes
    
    return percent_miss

def calculate_percent_miss_bi(genotypes_keys, alleles, num_genotypes):
    percent_miss = 0
    
    for genotype in genotypes_keys:
        if MISS in genotype:
            percent_miss += len(alleles[genotype])
    
    percent_miss = percent_miss * 1.0 / num_genotypes
    
    return percent_miss

##

def calculate_maf(alleles_keys, alleles, num_alleles, alleles_output):
    maf = 1.0
    
    if alleles_output == "mono":
        maf = calculate_maf_mono(alleles_keys, alleles, num_alleles)
    elif alleles_output == "bi":
        maf = calculate_maf_bi(alleles_keys, alleles, num_alleles)
    else:
        raise Exception("Unrecognized allele output option "+str(alleles_output)+".")
    
    return maf

def calculate_maf_mono(alleles_keys, alleles, num_alleles):
    maf = 1.0
    
    for allele in alleles_keys:
        if allele == HETERO or allele == MISS:
            continue
        else:
            curr = len(alleles[allele])
            if allele == REF: curr += 1
            curr = curr * 1.0 / num_alleles
            
            if curr < maf: maf = curr
    
    return maf

def calculate_maf_bi(alleles_keys, alleles, num_alleles):
    maf = 1.0
    
    for allele in alleles_keys:
        if allele == MISS:
            continue
        else:
            curr = 0
            for genotype in alleles:
                for allele_genotype in genotype.split(VCF_ALLELES_SEP):
                    if allele == allele_genotype:
                        curr+=1
            if allele == REF: curr += 2
            curr = curr * 1.0 / num_alleles
            
            if curr < maf: maf = curr
    
    return maf

##