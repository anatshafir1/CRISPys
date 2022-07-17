"""Candidate class file"""
from typing import Dict


class Candidate:
    """
        a class representing a sgRNA candidate to target CRISPR cut sites in a gene.
        a toy example:
        ATCCAATGTCTCTGGTTTT, 0.5, 0.19999999999999996, {'AT3G21670.1': 0.19999999999999996},
        {'AT3G21670.1': [['AATCCAACGTCTCTGGTTTT', {7: ['T', 'C']}]]}.
        Initiated attributes are explained under __init__ function.

        Attributes
        ----------

        num_of_genes_above_thr : int

        cleave_all_above_thr : int

        off_targets : bool
        """

    def __init__(self, seq: str, cut_expectation: float = 0.0, genes_score_dict: Dict = None,
                 mismatch_site_dict: Dict = None):
        """

        :param seq: the DNA sequence of the candidate
        :param cut_expectation: the probability that all the genes that expected to be cleaved by this
        sgRNA with high probability will be cleaved

        :param genes_score_dict: dict of gene name -> cleaving probability of this gene
        :param mismatch_site_dict: dict of gene name -> list of lists, each match sites and a mismatches dictionary.

        :return: a Candidate object
        """
        if genes_score_dict is None:
            genes_score_dict = dict()
        if mismatch_site_dict is None:
            mismatch_site_dict = dict()
        self.seq = seq
        self.cut_expectation = cut_expectation
        self.genes_score_dict = genes_score_dict
        self.targets_dict = mismatch_site_dict
        # initiated without parameters:
        self.num_of_genes_above_thr = 0
        self.cleave_all_above_thr = 1.0
        self.off_targets = False

    def fill_default_fields(self, gene_names):
        """
        For use when the sgRNA is constructed in the leaf of the BU tree

        :param gene_names: a list of gene names with a perfect match to the given sgRNA
        """
        self.genes_score_dict = dict()
        self.cut_expectation = 1.0
        for gene_name in gene_names:
            self.genes_score_dict[gene_name] = 1.0
            self.targets_dict = {gene_name: [[self.seq, {}]]}

    def __str__(self):
        return self.seq + ", " + str(self.cut_expectation) + ", " + str(self.genes_score_dict) + ", " + str(
            self.targets_dict)

    def __repr__(self):
        return self.__str__()

    def total_num_of_mismatches(self):
        """
        Counts the number of mismatches between a given candidate and its potential targets
        :return: number of mismatches
        :rtype: int
        """
        num_of_mismatches = 0
        for val in self.targets_dict.values():  # here the value is a list of targets. Each represented by list. Each
            # list: a match sites and a mismatches dictionary.
            for target in val:
                num_of_mismatches += len(target[1])
        return num_of_mismatches
