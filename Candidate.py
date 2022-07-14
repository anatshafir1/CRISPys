
class Candidate:

    # a toy example:
    # AATCCAATGTCTCTGGTTTT, 0.5, 0.19999999999999996, {'AT3G21670.1': 0.19999999999999996}, {'AT3G21670.1': [['AATCCAACGTCTCTGGTTTT', {7: ['T', 'C']}]]}

    def __init__(self, seq, cut_expectation=0, genes_score_dict=dict(), targets_dict=dict()):

        '''
        :param seq:
        :param fraction_of_cut: the fraction of the genes in the current node, expected to be cleaved with high probability
        :param cut_expectation: the probability that all the genes that expected to be cleaved by this sgRNA with high probability will be cleaved
        :param genes_score_dict: key: gene_name. value: cleaving pobability of this gene
        :param targets_dict (old name: match_sites_dict: key: gene_name. value: list of lists. each list: a match sites and a mismatches dictionary. right now it is not implemented like this at the Naive code- have to made deabuging, and to deaside maybe it is better not to save the mm locations in here.
        :param lowest_cut_site: the probability of cutting the site with the lowest cut probability by this sgRNA, among the sites of the sub tree
        width: multipication of cleaving oll of the tergets of the node in the tergets tree.
        :return:
        '''
        self.seq = seq
        self.cut_expectation = cut_expectation
        self.genes_score_dict = genes_score_dict
        self.targets_dict = targets_dict  # change the name to missmatch site dict.
        self.score_2 = None  # the score for the objective function
        self.num_of_genes_above_thr = 0
        self.cleave_all_above_thr = 1.0
        self.off_targets = False

    def fill_default_fields(self, gene_names):
        '''
        for use when the sgRNA is constructed in the leaf of the BU tree
        :param gene_names: a list of gene names with a perfect match to the given sgRNA
        :return:
        '''
        self.genes_score_dict = dict()
        self.cut_expectation = 1.0
        self.lowest_cut_site_prob = 1.0
        for gene_name in gene_names:
            self.genes_score_dict[gene_name] = 1.0
            self.targets_dict = {gene_name: [[self.seq, {}]]}

    def __str__(self):
        return self.seq + ", " + str(self.cut_expectation) + ", " + str(self.genes_score_dict) + ", " + str(
            self.targets_dict)

    def __repr__(self):
        return self.__str__()

    def total_num_of_mismatches(self):
        num_of_mm = 0
        for val in self.targets_dict.values():  # here the value is a list of targets. Each represented by list. Each list: a match sites and a mismatches dictionary.
            for target in val:
                num_of_mm += len(target[1])
        return num_of_mm
