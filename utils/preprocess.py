
import itertools
import json
import pickle
import random
import time
from datetime import datetime
from multiprocessing import Pool
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from fuzzysearch import find_near_matches
from utils.suffix_tree import Suffixtree
from utils.binary_search import get_similar_transcripts
from utils.suffix_build import build_suffix_array
import sys
from utils.send_mail import send_email, send_email_with_attachment
random.seed(42)




def find_similar_sequences(trigger: str, sequences_dict: dict[str, str]) -> Align.Alignment:
    """"function to find similar sequences to a trigger sequence
    :param trigger: trigger sequence
    :param sequences_dict: dictionary of sequences
    :return: alignment of the trigger sequence and the most similar sequence
    """""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    alignments = aligner.align(trigger, sequences_dict)
    return alignments[0]


def chunks(data_iter, size):
    """"function to split data into chunks for processing
    :param data_iter: data to be split
    :param size: size of the chunk
    :return: chunked data
    """""
    it = iter(data_iter)
    for i in range(0, len(data_iter), size):
        yield {k: data_iter[k] for k in itertools.islice(it, size)}


def run_chunks(data_path: str, func, args):
    """""function to run large data sets in chunks using multiprocessing
    :param data_path: path to the data file
    :param func: function to be applied to the data
    :param args: arguments to the function
    :return: None
    """""
    file_path = data_path
    with open(file_path, 'rb') as file:
        file_data = pickle.load(file)
    # Create enumeration for the data for reproducibility and file data indicator
    current_time = datetime.now().strftime("%H:%M:%S")
    enumerated_data = {key: i for i, key in enumerate(file_data.keys())}
    with open(f'enumerated_data_{current_time}.pkl', 'wb') as f:
        pickle.dump(enumerated_data, f)
    # Process the data in chunks
    n_process = 5
    for i, chunk in enumerate(chunks(file_data, n_process)):
        with Pool(n_process) as pool:
            args = []
            res = pool.map(func, args)
        # Combine the results and save of to intermediate file
        current_time = datetime.now().strftime("%H:%M:%S")
        if res:
            with open(f'intermediate_results_[{i * n_process}, {(i + 1) * n_process})_{current_time}.pkl', 'wb') as f:
                pickle.dump(res, f)


def find_matches_fuzzy(trigg: str, data_dict, edit_distance) -> dict[str, str]:
    """"function to find matches using fuzzy search
    :param edit_distance: maximum edit distance allowed
    :param trigg: trigger sequence
    :param data_dict: dictionary of sequences key=gene name, value= sequence
    :return: dictionary of sequences that match the trigger sequence
    """""
    sub_seqs = []
    for seq_name, seq in data_dict.items():
        seq_match_mapping = {"Protein": seq_name, "matched_seqs": [], "index's": []}
        try:
            matches = find_near_matches(trigg, seq, max_insertions=0, max_deletions=0, max_l_dist=edit_distance)
            for match_obj in matches:
                locus_start = match_obj.start
                end_locus = match_obj.end
                sub_seq = seq[locus_start:end_locus]
                seq_match_mapping["matched_seqs"].append(sub_seq)
                seq_match_mapping["index's"].append((locus_start, end_locus))

            if seq_match_mapping["matched_seqs"]:
                sub_seqs.append(seq_match_mapping)
        except Exception as e:
            print(f"Error in finding matches for {seq_name} with error {e}")
            return sub_seqs
    return sub_seqs

def find_with_suffix_tree(trigg: str, data_dict: dict[str, str]) -> dict[str, str]:
    """"function to find matches using suffix tree
    :param trigg: trigger sequence
    :param data_dict: dictionary of sequences key=gene name/rna name, value= rna/dna sequence
    :return: dictionary of sequences that match the trigger sequence
    """""
    sub_seqs_dict = {}
    for seq_name, seq in data_dict.items():
        tree = Suffixtree(seq, seq_name)
        try:
            matches = tree.search(trigg)
            for match in matches:
                sub_seq = (seq[match:match + len(trigg)], match)
                sub_seqs_dict[seq_name] = sub_seq
        except Exception as e:
            print(f"Error in finding matches for {seq_name} with error {e}")

    return sub_seqs_dict
def find_with_suffix_array(trigger, data_dict):
    """"function to find matches using suffix array"""
    suffixes_dict = []
    for gene_name, sequence in data_dict.items():
        try:
            sequence = sequence
            protein = "protein"
            gene = gene_name
            suffixes = [(sequence[i:], i) for i in range(len(sequence))]
            suffixes.sort(key=lambda x: x[0])
            suffixes_dict.append({'gene': gene, 'protein': protein, 'array': [suffix[1] for suffix in suffixes], 'sequence': sequence})
        except:
            continue

    print('getting similar transcripts')
    off_target = get_similar_transcripts(trigger, suffixes_dict, 2)
    return off_target




def main(gene, trigger, rna, cell_type, file_dict,email):
    print(f"Gene: {gene}, Reporter: {trigger}, RNA: {rna}, Cell Type: {cell_type}")
    """"main function to combine Peleg's code with user request
    output to be sent by email 
    """""
    if file_dict != "EMPTY":
        try:
            data_dict = json.loads(file_dict)
        except json.JSONDecodeError:
            data_dict = {}
    else:
        data_dict = {}
    print(f"Data dict: {data_dict}")
    print("processing off target sequences for the trigger sequence")
    regular = False
    if regular:
        off_target = find_matches_fuzzy(trigger, data_dict, 2)
    else:
        off_target = find_with_suffix_array(trigger, data_dict)
    body = "Thank you for using our service. We are processing your request and will send you the results shortly. \n\nPROtech Team"
    temp = ''
    email_msg = send_email_with_attachment(email, "Processed Files", body,[temp])
    send_email(email_msg)
    return



if __name__ == '__main__':
    # Get the arguments from the user form and send to main function
    # 1. no trigger and gene present-> find window in gene -> find similar sequences in the off target -> make toehold.
    # 2. trigger and gene non relevant -> find similar sequences in the off target with trigger -> make toehold.
        # if file exist -> use the file to off target
        # if file does not exist -> use the data base default off target

    """
    gene = sys.argv[1]  # pick window
    trigger = sys.argv[2] # toehold

    reporter = sys.argv[3]
    file_dict = sys.argv[6]
    
    cell_type = sys.argv[4] # cell type: how to deal with reporter on different organisem?
                            # ribozom starting site
                            # off ratio -> similar sequences
    email = sys.argv[5]
    main(gene, trigger, reporter, cell_type, file_dict, email)
    """
