import pickle
import random
import numpy as np
import time
from itertools import combinations
import matplotlib.pyplot as plt
from utils.preprocess import find_matches_fuzzy, find_near_matches
random.seed(42)

def homo_sapiance_search_plot():
    path = '/Users/netanelerlich/PycharmProjects/webTool/data/transcripts_data.pkl'
    with open(path, 'rb') as f:
        data = pickle.load(f)
    times = []
    for i in range(10, 25, 1):
        trigger = ''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(i)])
        print(f"Processing trigger of length {i}")
        genes = list(data.keys())
        start = time.time()
        matches_over_genes = {}
        for gene in genes:
            print(f"Processing gene {gene}")
            transcripts_list = data[gene]  # list ot dict
            input_transcripts = {}
            for transcript_dict in transcripts_list:
                protein_name = transcript_dict['Protein']  # protein name
                sequence = transcript_dict['Sequence']  # sequence
                input_transcripts[protein_name] = sequence
            res = find_matches_fuzzy(trigger, input_transcripts, 2)
            if res:
                matches_over_genes[gene] = res
        end = time.time()
        times.append(end - start)

    plt.plot(list(range(10, 25, 1)), times)
    plt.title('Transcriptome search time for different trigger lengths using fuzzy search')
    plt.xlabel('Length of trigger')
    plt.ylabel('Time in seconds')
    plt.show()

def gen_mutations(trigger: str, sub: int) -> dict[str, int]:
    trig_len = len(trigger)
    indexes_permutation = list(combinations(range(trig_len), sub))
    nucleotides = ['A', 'C', 'G', 'T']
    trigger_muts = {}

    for mutated_index in indexes_permutation:
        trig_char_list = list(trigger)
        for nuc_idx in mutated_index:
            choices = nucleotides.copy()
            nuc = trig_char_list[nuc_idx]
            choices.remove(nuc)
            trig_char_list[nuc_idx] = random.choice(choices)
        mutated_trigger = ''.join(trig_char_list)
        trigger_muts[mutated_trigger] = mutated_index

    return trigger_muts
def construct_dummy_seq(trigger_muts: dict[str, int]) -> tuple[str, dict[str, tuple[int, int, int]]]:
    trig_n = len(trigger_muts)
    n_list = random.sample(range(1, trig_n * 10), trig_n)
    fillers = []
    for x in n_list:
        filler_seq = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(x))
        fillers.append(filler_seq)

    dummy_seq = ''
    mut_trigger_mapping = {}
    for trig, mut_index in trigger_muts.items():
        dummy_seq += fillers.pop(0)
        start_index = len(dummy_seq)
        dummy_seq += trig
        end_index = len(dummy_seq)
        mut_trigger_mapping[trig] = (start_index, end_index, mut_index)

    return dummy_seq, mut_trigger_mapping


def test_fuzzy_search(trig_len: int, sub: int):
    trigger = ''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(trig_len)])
    mutated_triggers = gen_mutations(trigger, sub)
    test_seq, trig_muts = construct_dummy_seq(mutated_triggers)
    test = find_near_matches(trigger, test_seq, max_deletions=0, max_insertions=0, max_substitutions=sub, max_l_dist=sub)

    test_res = {}
    for match in test:
        start = match.start
        end = match.end
        sub_seq = test_seq[start:end]
        test_res[sub_seq] = (start, end)
    return test_res, trig_muts, test_seq



if __name__ == "__main__":
        subs = 5
        triggers_lens = 25
        res_array = []

        for trig_len in range(15, 25, 2):
            for subs in range(5):
                trig_len = trig_len
                test_res, mutations_mapping, dummy_seq = test_fuzzy_search(trig_len=trig_len, sub=subs)
                for mutated_trigger, mut_loc in mutations_mapping.items():
                    for search_res_trigger, res_loc in test_res.items():
                        if mut_loc[0] == res_loc[0] and mut_loc[1] == res_loc[1]:
                            mutations_mapping[mutated_trigger] = all(np.array(list(mutated_trigger)) == np.array(list(search_res_trigger)))

                print(f'for trigger len={trig_len} with {subs} subs, all possible mutated triggers found: {all(mutations_mapping.values())}')
                res_array.append(f'for trigger len={trig_len} with {subs} subs, all possible mutated triggers found: {all(mutations_mapping.values())}')









