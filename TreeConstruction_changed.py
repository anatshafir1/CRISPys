# taken from http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction-pysrc.html#DistanceCalculator
# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its 
# license. Please see the LICENSE file that should have been included 
# as part of this package. 

"""Classes and methods for tree construction"""
__docformat__ = "restructuredtext en"

import copy
from typing import List
from Bio.Phylo import BaseTree, TreeConstruction


class DistanceTreeConstructor(TreeConstruction.TreeConstructor):
    """Distance based tree constructor.

    :Parameters:
        method : str
            Distance tree construction method, 'upgma'.
        distance_calculator : DistanceCalculator
            The distance matrix calculator for multiple sequence alignment.
            It must be provided if `build_tree` will be called.

    Example
    --------

    >> from TreeConstruction import DistanceTreeConstructor
    >> constructor = DistanceTreeConstructor()

    UPGMA Tree:

    >> upgmatree = constructor.upgma(dm)
    >> print upgmatree
    Tree(rooted=True)
        Clade(name='Inner4')
            Clade(branch_length=0.171955155115, name='Inner1')
                Clade(branch_length=0.111111111111, name='Epsilon')
                Clade(branch_length=0.111111111111, name='Delta')
            Clade(branch_length=0.0673103855608, name='Inner3')
                Clade(branch_length=0.0907558806655, name='Inner2')
                    Clade(branch_length=0.125, name='Gamma')
                    Clade(branch_length=0.125, name='Beta')
                Clade(branch_length=0.215755880666, name='Alpha')
    """

    methods = ['upgma']

    def __init__(self, distance_calculator=None, method="upgma"):
        if (distance_calculator is None
                or isinstance(distance_calculator, TreeConstruction.DistanceCalculator)):
            self.distance_calculator = distance_calculator
        else:
            raise TypeError("Must provide a DistanceCalculator object.")
        if isinstance(method, str) and method in self.methods:
            self.method = method
        else:
            raise TypeError("Bad method: " + method +
                            ". Available methods: " + ", ".join(self.methods))

    def build_tree(self, msa):
        if self.distance_calculator:
            dm = self.distance_calculator.get_distance(msa)
            tree = None
            if self.method == 'upgma':
                tree = self.upgma(dm)
            return tree
        else:
            raise TypeError("Must provide a DistanceCalculator object.")

    def upgma(self, distance_matrix):
        """Construct and return an UPGMA tree.

        Constructs and returns an Unweighted Pair Group Method
        with Arithmetic mean (UPGMA) tree.

        :Parameters:
            distance_matrix : _DistanceMatrix
                The distance matrix for tree construction.
        """
        if not isinstance(distance_matrix, TreeConstruction._DistanceMatrix):
            raise TypeError("Must provide a _DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        weights = [1 for i in range(len(dm))]  # number of items at eah node
        # init terminal clades
        clades = [CladeNew(None, name) for name in dm.names]
        leaves_lst = []
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        inner_clade = None
        while len(dm) > 1:
            min_dist = dm[1, 0]
            # find minimum index
            for i in range(1, len(dm)):
                for j in range(0, i):
                    if min_dist >= dm[i, j]:
                        min_dist = dm[i, j]
                        min_i = i
                        min_j = j

            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = CladeNew(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            clade1.set_parent(inner_clade)
            inner_clade.clades.append(clade2)
            clade2.set_parent(inner_clade)
            # assign branch length
            if clade1.is_terminal():
                clade1.branch_length = min_dist * 1.0 / 2
            else:
                clade1.branch_length = min_dist * 1.0 / 2 - self._height_of(clade1)

            if clade2.is_terminal():
                clade2.branch_length = min_dist * 1.0 / 2
            else:
                clade2.branch_length = min_dist * 1.0 / 2 - self._height_of(clade2)

            # update nodes list
            if "Inner" not in clades[min_j].name:
                leaves_lst.append(clades[min_j])
            clades[min_j] = inner_clade
            if "Inner" not in clades[min_i].name:
                leaves_lst.append(clades[min_i])
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            new_weight = weights[min_i] + weights[min_j]
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] * weights[min_i] / new_weight + dm[min_j, k] * weights[
                        min_j] / new_weight)

            dm.names[min_j] = "Inner" + str(inner_count)
            weights[min_j] = new_weight
            del dm[min_i]
            del weights[min_i]

            dm.names[min_j] = "Inner" + str(inner_count)

        inner_clade.branch_length = 0
        return Tree_new(inner_clade, leaves=leaves_lst)

    def _height_of(self, clade):
        """calculate clade height -- the longest path to any terminal."""
        height = 0
        if clade.is_terminal():
            height = clade.branch_length
        else:
            height = height + max([self._height_of(c) for c in clade.clades])
        return height


class CladeNew(BaseTree.Clade):
    """
    Implementation of BaseTree.Clade for the purposes of CRISPys algorithm
    """

    def __init__(self, branch_length: int = None, name: str = None, clades: List = None, confidence=None,
                 color: BaseTree.BranchColor = 'w', width=None):
        """

        """
        super().__init__(branch_length, name, clades, confidence, color, width)
        self.parent = None
        self.node_targets = list()
        self.candidates = dict()
        self.polymorphic_sites = set()

    def set_parent(self, parent):
        self.parent = parent

    def add_node_target(self, leaf):
        self.node_targets.append(leaf)

    def fill_polymorphic_sites(self, polymorphic_sites: set = None):
        self.polymorphic_sites = polymorphic_sites


class Tree_new(BaseTree.Tree):
    def __init__(self, root=None, rooted=True, id=None, name=None, leaves=None):
        BaseTree.Tree.__init__(self, root, rooted, id, name)
        self.leaves = leaves
