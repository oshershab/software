from suffix_trees import STree
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class SuffixTree:

    def __init__(self, text, gene_name):

        self.name = gene_name
        self.tree = STree.STree(text)

    def search(self, pattern):
        return self.tree.find_all(pattern)


if __name__ == "__main__":
    text = "ACATGATAGCATCGCGATAGCTAGAT"
    tree = SuffixTree(text, 'dummy_gene_name')
    print(tree.search("ATAG"))
    tree= STree.STree("ACATGATAGCATCGCGATAGCTAGAT")





















