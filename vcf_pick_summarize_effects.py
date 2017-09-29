#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2017

import sys, traceback

SNPEFF_FIELD = "EFF"
SNPEFF_FIELD_EFF_FIELD_IMPACT = 0 # High, Moderate, Low, Modifier
SNPEFF_FIELD_EFF_FIELD_CLASS = 1 # NONE, SILENT, MISSENSE, NONSENSE
SNPEFF_FIELD_EFF_FIELD_CODON_CHANGE = 2
SNPEFF_FIELD_EFF_FIELD_AACHANGE = 3
SNPEFF_FIELD_EFF_FIELD_PROT_LENGTH = 4
SNPEFF_FIELD_EFF_FIELD_GENE = 5
SNPEFF_FIELD_EFF_FIELD_BIOTYPE = 6
SNPEFF_FIELD_EFF_FIELD_ISOF = 7
SNPEFF_FIELD_EFF_FIELD_GENOTYPE = 8 # 0: ref, 1: alt, ...(more for indels)
SNPEFF_FIELD_EFF_FIELD_ERRORS = 9

# deleterious order (from http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf)
del_order = {
"chromosome_number_variation":1,
"exon_loss_variant":2,
"exon_loss":2, # I added this one
"frameshift_variant":3,
"stop_gained":4,
"stop_lost":5,
"start_lost":6,
"splice_acceptor_variant":7,
"splice_donor_variant":8,
"rare_amino_acid_variant":9,
"missense_variant":10,
"inframe_insertion":11,
"disruptive_inframe_insertion":12,
"inframe_deletion":13,
"disruptive_inframe_deletion":14,
"5_prime_UTR_truncation+exon_loss_variant":15,
"3_prime_UTR_truncation+exon_loss":16,
"splice_branch_variant":17,
"splice_region_variant":18,
"splice_branch_variant":19,
"stop_retained_variant":20,
"initiator_codon_variant":21,
"synonymous_variant":22,
"initiator_codon_variant+non_canonical_start_codon":23,
"non_canonical_start_codon":23.5, # I added this one
"stop_retained_variant":24,
"coding_sequence_variant":25,
"5_prime_UTR_variant":26,
"3_prime_UTR_variant":27,
"5_prime_UTR_premature_start_codon_gain_variant":28,
"3_prime_UTR_truncation":28.5, # I added this one
"5_prime_UTR_truncation":28.6, # I added this one
"upstream_gene_variant":29,
"downstream_gene_variant":30,
"TF_binding_site_variant":31,
"regulatory_region_variant":32,
"miRNA":33,
"custom":34,
"sequence_feature":35,
"conserved_intron_variant":36,
"intron_variant":37,
"intragenic_variant":38,
"conserved_intergenic_variant":39,
"intergenic_region":40,
"coding_sequence_variant":41,
"non_coding_exon_variant":42,
"nc_transcript_variant":43,
"gene_variant":44,
"chromosome":45
}

LOF_FIELD = "LOF"
LOF_FIELDS_GENE = 0
LOF_FIELDS_ID = 1
LOF_FIELDS_NUM_TRANSCRIPTS = 2
LOF_FIELDS_PERCEN_TRANSCRIPTS = 3

NMD_FIELD = "NMD"
NMD_FIELDS_GENE = 0
NMD_FIELDS_ID = 1
NMD_FIELDS_NUM_TRANSCRIPTS = 2
NMD_FIELDS_PERCEN_TRANSCRIPTS = 3

indels_dict = {}
snps_dict = {}

vcf_filename = sys.argv[1]

sys.stderr.write("Input file: "+vcf_filename+"\n")

total_variants = 0
try:
    vcf_file = open(vcf_filename)
    for line in vcf_file:
            #code
            if line.startswith("#"): continue
            
            total_variants += 1
            line_data = line.strip().split()
            #print line_data
            
            # Check if variant is indel or snp
            ref_allele = line_data[3]
            alt_allele = line_data[4]
            
            is_indel = False
            curr_dict = snps_dict
            if len(ref_allele) > 1 or len(alt_allele) > 1:
                curr_dict = indels_dict
                is_indel = True
                
            # Check if variant is unique
            var_id = "-".join(line_data[:5]) # contig-pos-.-ref-alt
            if var_id in curr_dict:
                print var_id
                print curr_dict[var_id]
                raise Exception("Var "+var_id+" already exists in current dict.")
            else:
                var_dict = {"var_id":var_id, "genes":{}}
                curr_dict[var_id] = var_dict
            
            # Parse INFO field to info_fields_dict
            info = line_data[7]
            #print "INFO"
            #print info
            
            info_fields = info.split(";")
            info_fields_dict = {}
            for info_field in info_fields:
                info_field_split = info_field.split("=")
                if len(info_field_split) == 2:
                    info_fields_dict[info_field_split[0]] = info_field_split[1]
                elif len(info_field_split) == 1:
                    info_fields_dict[info_field_split[0]] = "novalue"
                else:
                    raise Exception("Field in info should have 'FIELD=VALUE' or 'FIELD' format")
            
            #print "DICT"
            #print info_fields_dict
            
            if SNPEFF_FIELD in info_fields_dict:
                snpeff_field = info_fields_dict[SNPEFF_FIELD]
                snpeff_field_effs = snpeff_field.split(",")
                for snpeff_eff in snpeff_field_effs:
                    snpeff_eff_parts = snpeff_eff.split("(")
                    snpeff_eff_type = snpeff_eff_parts[0]
                    # if there are several snpeff types (type1+type2)
                    # I take the del_order which is smaller
                    if "+" in snpeff_eff_type:
                        snpeff_eff_del_order = 99999
                        for snpeff_eff_subtype in snpeff_eff_type.split("+"):
                            snpeff_eff_subtype_del_order = del_order[snpeff_eff_subtype]
                            if snpeff_eff_subtype_del_order < snpeff_eff_del_order:
                                snpeff_eff_del_order = snpeff_eff_subtype_del_order
                    else:
                        snpeff_eff_del_order = del_order[snpeff_eff_type]
                    
                    #print str(snpeff_eff_type)+" - "+str(snpeff_eff_del_order)
                    
                    snpeff_eff_fields = snpeff_eff_parts[1].replace(")", "").split("|")
                    snpeff_eff_gene = snpeff_eff_fields[SNPEFF_FIELD_EFF_FIELD_GENE]
                    #print snpeff_eff_gene
                    if snpeff_eff_gene == "":
                        snpeff_eff_gene = snpeff_eff_type
                        
                    if snpeff_eff_gene in var_dict["genes"]:
                        snpeff_eff_gene_dict = var_dict["genes"][snpeff_eff_gene]
                        if snpeff_eff_del_order < snpeff_eff_gene_dict["del_order"]:
                            snpeff_eff_gene_dict["type"] = snpeff_eff_type
                            snpeff_eff_gene_dict["del_order"] = snpeff_eff_del_order
                    else:
                        snpeff_eff_gene_dict = {"type":snpeff_eff_type,
                                                "del_order":snpeff_eff_del_order,
                                                "lof":False,
                                                "nmd":False}
                        var_dict["genes"][snpeff_eff_gene] = snpeff_eff_gene_dict
                        
            
            if LOF_FIELD in info_fields_dict:
                lof_fields = info_fields_dict[LOF_FIELD]
            #    print ""
            #    print info
            #    print lof_fields
                for lof_fields_lof in lof_fields.split("),"):
                    lof_fields_lof_fields = lof_fields_lof.replace("(", "").replace(")", "").split("|")
                    #print lof_fields_lof_fields
                    lof_gene = lof_fields_lof_fields[LOF_FIELDS_GENE]
                    #print lof_gene
                    
                    if lof_gene in var_dict["genes"]:
                        var_dict["genes"][lof_gene]["lof"] = True
                    else:
                        print info
                        print var_dict
                        raise Exception("LOF gene "+lof_gene+" not in variants dict.")
            
            if NMD_FIELD in info_fields_dict:
                nmd_fields = info_fields_dict[NMD_FIELD]
                #print ""
                #print info
                #print nmd_fields
                for nmd_fields_nmd in nmd_fields.split("),"):
                    nmd_fields_nmd_fields = nmd_fields_nmd.replace("(", "").replace(")", "").split("|")
                    nmd_gene = nmd_fields_nmd_fields[NMD_FIELDS_GENE]
                    if nmd_gene in var_dict["genes"]:
                        var_dict["genes"][nmd_gene]["nmd"] = True
                    else:
                        raise Exception("NMD gene "+nmd_gene+" not in variants dict.")
                
            
            #print info
            #print var_dict
            
except Exception as e:
    sys.stderr.write("An error ocurred while parsing VCF file.\n")
    #print line
    traceback.print_exc()
    sys.exit(-1)
finally:
    vcf_file.close()

snps_types_dict = {}
snps_lof = 0
snps_nmd = 0
snps_total = 0
for var_id in snps_dict:
    var_dict = snps_dict[var_id]
    snps_total += 1
    for gene in var_dict["genes"]:
        gene_dict = var_dict["genes"][gene]
        snpeff_type = gene_dict["type"]
        lof = gene_dict["lof"]
        nmd = gene_dict["nmd"]
        
        if snpeff_type in snps_types_dict:
            snps_types_dict[snpeff_type] += 1
        else:
            snps_types_dict[snpeff_type] = 1
        
        if lof:
            snps_lof += 1
        
        if nmd:
            snps_nmd += 1

indels_types_dict = {}
indels_lof = 0
indels_nmd = 0
indels_total = 0
for var_id in indels_dict:
    var_dict = indels_dict[var_id]
    indels_total += 1
    for gene in var_dict["genes"]:
        gene_dict = var_dict["genes"][gene]
        snpeff_type = gene_dict["type"]
        lof = gene_dict["lof"]
        nmd = gene_dict["nmd"]
        
        if snpeff_type in indels_types_dict:
            indels_types_dict[snpeff_type] += 1
        else:
            indels_types_dict[snpeff_type] = 1
        
        if lof:
            indels_lof += 1
        
        if nmd:
            indels_nmd += 1

print "#\tSNPs\t"+str(snps_total)+"\t"+str(snps_lof)+"\t"+str(snps_nmd)
# Print SNPs effs
for snpeff_type in sorted(del_order, key=lambda x: del_order[x]):
    if snpeff_type in snps_types_dict:
        print ">\t"+snpeff_type+"\t"+str(snps_types_dict[snpeff_type])

print "#\tIndels\t"+str(indels_total)+"\t"+str(indels_lof)+"\t"+str(indels_nmd)
# Print Indels effs
for snpeff_type in sorted(del_order, key=lambda x: del_order[x]):
    if snpeff_type in indels_types_dict:
        print ">\t"+snpeff_type+"\t"+str(indels_types_dict[snpeff_type])


print "#\tALL"
# Print SNPs effs
for snpeff_type in sorted(del_order, key=lambda x: del_order[x]):
    show = False
    snps_count = 0
    indels_count = 0
    if snpeff_type in snps_types_dict:
        show = True
        snps_count = snps_types_dict[snpeff_type]
        
    if snpeff_type in indels_types_dict:
        show = True
        indels_count = indels_types_dict[snpeff_type]
    
    if show:
        print ">\t"+snpeff_type+"\t"+str(snps_count)+"\t"+str(indels_count)

## END