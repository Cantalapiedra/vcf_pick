#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

import os, sys, tempfile, csv
from subprocess import Popen, PIPE

def cluster_table(temp_file):
    curr_path = os.path.dirname(os.path.abspath(__file__))
    __command = " ".join([curr_path+"/cluster_samples.R ", temp_file])
    
    sys.stderr.write("\tcluster_samples.R: "+temp_file+"\n")
    
    p = Popen(__command, shell=True, stdout=PIPE, stderr=PIPE)
    com_list = p.communicate()
    output = com_list[0]
    output_err = com_list[1]
    retValue = p.returncode
    
    if retValue != 0: raise Exception("cluster_samples.R: return != 0. "+__command+"\n"+str(output_err)+"\n")
    
    #sys.stderr.write("cluster_samples.R: return value "+str(retValue)+"\n")
    
    return output

def f_cluster_samples(samples_rows, biallelic):
    new_samples_list = []
    
    sys.stderr.write("Clustering samples...\n")
    if biallelic == "bi":
        raise Exception("Clustering biallelic alleles not yet allowed. Use monoallelic option instead.")
    else: # biallelic == "mono":
        try:
            ##
            ## Create temp file
            file_tmp = tempfile.NamedTemporaryFile()
            file_name = file_tmp.name
            
            sys.stderr.write("\ttemp file created: "+file_name+"\n")
            
            writer = csv.writer(file_tmp, dialect='excel', delimiter="\t")
            
            for sample_row in samples_rows:
                writer.writerow(sample_row)
            
            file_tmp.flush()
            
            ##
            ## Cluster from temp file
            output = cluster_table(file_name)
            
            ##
            ## Recover order of samples
            for line in output.strip().split("\n"):
                new_samples_list.append(line.strip())
            
        except Exception:
            raise
        finally:
            file_tmp.close()
    
    sys.stderr.write("Samples clustered.\n")
    
    return new_samples_list

##

def parse_queries_file(query_file, keys=(0,)):
    query_list = []
    if query_file == "":
        pass
    else:
        try:
            queries_file = open(query_file, 'r')
            for line in queries_file:
                line_data = line.strip().split()
                if len(keys) == 1:
                    query = line_data[keys[0]]
                else:
                    query = []
                    for key in keys:
                        query.append(line_data[key - 1])
                
                query_list.append(query)
            
        except Exception:
            raise
        finally:
            queries_file.close()
    
    return query_list

##

def parse_samples_list(samples_filename):
    samples_list = []
    if samples_filename == "":
        pass
    else:
        try:
            samples_file = open(samples_filename, 'r')
            for line in samples_file:
                sample = line.strip().split()[0]
                if sample not in set(samples_list):
                    samples_list.append(sample)
            
        except Exception:
            raise
        finally:
            samples_file.close()
    
    return samples_list

##

def parse_samples_translation(samples_translation):
    samples_trans_dict = {}
    if samples_translation == "":
        pass
    else:
        try:
            samples_file = open(samples_translation, 'r')
            for sample in samples_file:
                sample_data = sample.strip().split()
                samples_trans_dict[sample_data[0]] = sample_data[1]
            
        except Exception:
            raise
        finally:
            samples_file.close()
    
    return samples_trans_dict

##