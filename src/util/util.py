#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra - EEAD - CSIC - 2016

def parse_queries_file(query_file, keys=(0,)):
    query_list = []
    if query_file == "":
        raise Exception("No queries especified.")
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