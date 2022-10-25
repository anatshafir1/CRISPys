"""MOFF implementation for CRISPys"""

# ######### Import all the packages ########## #
from itertools import combinations
from itertools import product
import numpy as np
from scipy.stats import gmean
import globals
#################################


def OneHotEncoding(sg_ls):
    """This function encodes the sgRNA sequences into (16,19) vector with 0,1 presentation of certain dinucleotide at
    certain position.

       o input:  sg_ls: A list of sgRNAs for one-hot encoding

       o Output: The function will return a numpy array of feature matrix for machine learning.
    """
    di_ls = [s[0] + s[1] for s in list(product('ATGC', repeat=2))]  # Get all possible di-nucleotide combinations
    di_dic = {}
    for i in range(len(di_ls)):
        di_dic[di_ls[i]] = i  # Record the index of each dinucleotide

    ls_all = []  # Initialize a list to add vector for different sgRNA
    for sg in sg_ls:
        vec_all = []
        for i in range(len(sg) - 1):
            vec = [0] * len(di_ls)  # Make all position to be 0
            di = sg[i:i + 2]
            vec[di_dic[di]] = 1  # Assign 1 if certain dinucleotide appear at certain position
            vec_all.append(vec)

        ls_all.append(np.array(vec_all).T)

    return np.array(ls_all)


def GetMutType(s1, s2):
    """This function is obtain the mismatches between sgRNA and the target

       o input:  1). s1: the sequence of sgRNA; 2). s2: the sequence of target DNA

       o Output: The function will return: 1). A list of positions where mismatch happen.
                                           2). A list of mismatch types at each position.
    """
    pos_ls = []
    mut_ls = []
    for i in range(20):  # Go through the index along the 20bp sequence
        d = ''
        if s1[i] != s2[i]:
            pos = 20 - i  # The index relative to PAM
            if s1[i] == 'T':
                r = 'U'  # Replace 'T' with 'U' in sgRNA.
            else:
                r = s1[i]
            # Get mutation type given the nt at sgRNA and target
            if s2[i] == 'A':
                d = 'T'
            elif s2[i] == 'T':
                d = 'A'
            elif s2[i] == 'C':
                d = 'G'
            elif s2[i] == 'G':
                d = 'C'
            elif s2[i] == 'N':
                d = s1[i]
            pos_ls.append(pos)
            mut_ls.append('p' + str(pos) + 'r' + r + 'd' + d)  # p3rAdC: mismatch A-G at index 3 to PAM
    return pos_ls, mut_ls


def Multiply(m1_dic, sg_ls, tg_ls):
    """This function calculates the off-target effect by multiplying the MDE at each position

       o input:1). m1_dic: Python dic contains MDE of all the possible nucleotide mismatches (12)
                   at all possible positions (20)
               2). sg_ls: A list of sgRNAs
               3). tg_ls: A list of DNA targets

       o Output: A list of calculated mismatch-dependent effect.
    """
    me_ls = []
    for i in range(len(sg_ls)):
        s1 = sg_ls[i][0:20].upper()
        s2 = tg_ls[i][0:20].upper()
        mut_ls = GetMutType(s1, s2)[1]
        score = 1
        for mut in mut_ls:  # Multiply all the 1-mismatch effects
            score = score * m1_dic[mut]  # m1_dic: dic file
        me_ls.append(score)

    return me_ls


def MisNum(sg_ls, tg_ls):
    """This function gets mismatch numbers of gRNA-target pairs"""
    num_ls = []
    for i in range(len(sg_ls)):
        s1 = sg_ls[i][0:20].upper()
        s2 = tg_ls[i][0:20].upper()

        num = len(GetMutType(s1, s2)[0])
        num_ls.append(num)

    return num_ls


def CombineGM(m2_dic, sg_ls, tg_ls):
    """This function is Calculate Combinatorial effect (CE) for given mismatch positions

       o input:1). m2_dic: Python dic contains CE of all the possible position combinations
               2). sg_ls: A list of sgRNAs
               3). tg_ls: A list of DNA targets

       o Output: A list of calculated combinatorial effects.
    """
    cm_ls = []
    for i in range(len(sg_ls)):
        s1 = sg_ls[i][0:20].upper()
        s2 = tg_ls[i][0:20].upper()
        pos_ls = sorted(GetMutType(s1, s2)[0])

        # Combinatorial effect at certain position combination.
        di_ls = list(combinations(pos_ls, 2))
        c_ls = [m2_dic[str(di[0]) + '&' + str(di[1])] for di in di_ls]

        if len(pos_ls) > 1:
            m = gmean(c_ls) ** (len(pos_ls) - 1)  # Geometric mean of all possible combinations
        else:
            m = 1
        cm_ls.append(m)
    return cm_ls


def MOFF_score(m1_dic, m2_dic, sg_ls, tg_ls):
    """This function is predict off-target MOFF score for given gRNA-target pairs

       o input:1). m2_dic: Python dic contains CE of all the possible position combinations (20*19)
               2). m1_dic: Python dic contains MDE of all the possible nucleotide mismatches (12)
                   at all possible positions (20)
               3). df: A panda dataframe with one column of sgRNA and another column of DNA targets

       o Output: A panda dataframe with off-target predictions using different models (factors)
    """

    np.random.seed(24)  # for reproducibility
    pred_test = list(globals.moff_loaded_model.predict(OneHotEncoding([s.upper()[0:20] for s in sg_ls])))

    gop = [g[0] for g in pred_test]
    mde = Multiply(m1_dic, sg_ls, tg_ls)
    ce = CombineGM(m2_dic, sg_ls, tg_ls)
    mms = MisNum(sg_ls, tg_ls)
    gmt = [a ** b for a, b in zip(gop, mms)]
    moff = [a * b * c for a, b, c in zip(mde, ce, gmt)]
    return moff
