
import numpy as np
import pickle
import pandas as pd
import os
import re
import pickle

file_path = "C:/Users/Tal/OneDrive/Desktop/Igem/webTool/search/cds_from_genomic.fna"

with open(file_path, "r") as file:
    annotations = file.read()
    
transcriptome = annotations.split('>lcl')[1:]

protein_pattern = r'\[protein=(.*?)\]'
gene_pattern = r'\[gene=(.*?)\]'
sequence_pattern = re.compile(r'([ACGT]+)')

def build_suffix_array(elements):
    suffixes_dict = []
    for item in elements:
        try:
            sequence = item.split('[gbkey=CDS]')[1].replace('\n', '')
            protein = re.findall(protein_pattern, item)[0]
            gene = re.findall(gene_pattern, item)[0]
            suffixes = [(sequence[i:], i) for i in range(len(sequence))]
            suffixes.sort(key=lambda x: x[0])
            suffixes_dict.append({'gene': gene, 'protein': protein, 'array': [suffix[1] for suffix in suffixes], 'sequence': sequence})
        except:
            continue
    return suffixes_dict

suffixes = build_suffix_array(transcriptome)

with open('Escherichia_coli_ASM886v2.pkl', 'wb') as f:
    pickle.dump(suffixes, f)