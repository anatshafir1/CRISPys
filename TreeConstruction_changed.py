# taken from http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction-pysrc.html#DistanceCalculator
# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its 
# license. Please see the LICENSE file that should have been included 
# as part of this package. 

"""Classes and methods for tree construction"""
__docformat__ = "restructuredtext en"

import itertools
import copy
from Bio.Phylo import BaseTree, TreeConstruction
from Bio.Align import MultipleSeqAlignment
from Bio.SubsMat import MatrixInfo
from Bio import _py3k


def _is_numeric(x):
    return _py3k._is_int_or_long(x) or isinstance(x, (float, complex))


class _Matrix(object):
    """Base class for distance matrix or scoring matrix 
 
    Accepts a list of names and a lower triangular matrix.:: 
 
        matrix = [[0], 
                  [1, 0], 
                  [2, 3, 0], 
                  [4, 5, 6, 0]] 
        represents the symmetric matrix of 
        [0,1,2,4] 
        [1,0,3,5] 
        [2,3,0,6] 
        [4,5,6,0] 
 
    :Parameters: 
        names : list 
            names of elements, used for indexing 
        matrix : list 
            nested list of numerical lists in lower triangular format 
 
    Example 
    ------- 
 
    >> from Bio.Phylo.TreeConstruction import _Matrix 
    >> names = ['Alpha', 'Beta', 'Gamma', 'Delta'] 
    >> matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]] 
    >> m = _Matrix(names, matrix) 
    >> m 
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]) 
 
    You can use two indices to get or assign an element in the matrix. 
 
    >> m[1,2] 
    3 
    >> m['Beta','Gamma'] 
    3 
    >> m['Beta','Gamma'] = 4 
    >> m['Beta','Gamma'] 
    4 
 
    Further more, you can use one index to get or assign a list of elements related to that index. 
 
    >> m[0] 
    [0, 1, 2, 4] 
    >> m['Alpha'] 
    [0, 1, 2, 4] 
    >> m['Alpha'] = [0, 7, 8, 9] 
    >> m[0] 
    [0, 7, 8, 9] 
    >> m[0,1] 
    7 
 
    Also you can delete or insert a column&row of elemets by index. 
 
    >> m 
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]]) 
    >> del m['Alpha'] 
    >> m 
    _Matrix(names=['Beta', 'Gamma', 'Delta'], matrix=[[0], [4, 0], [5, 6, 0]]) 
    >> m.insert('Alpha', [0, 7, 8, 9] , 0) 
    >> m 
    _Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]]) 
 
    """

    def __init__(self, names, matrix=None):
        """Initialize matrix by a list of names and a list of 
        lower triangular matrix data"""
        # check names 
        if isinstance(names, list) and all(isinstance(s, str) for s in names):
            if len(set(names)) == len(names):
                self.names = names
            else:
                raise ValueError("Duplicate names found")
        else:
            raise TypeError("'names' should be a list of strings")

            # check matrix
        if matrix is None:
            # create a new one with 0 if matrix is not assigned 
            matrix = [[0] * i for i in range(1, len(self) + 1)]
            self.matrix = matrix
        else:
            # check if all elements are numbers 
            if (isinstance(matrix, list)
                    and all(isinstance(l, list) for l in matrix)
                    and all(_is_numeric(n) for n in [item for sublist in matrix
                                                     for item in sublist])):
                # check if the same length with names 
                if len(matrix) == len(names):
                    # check if is lower triangle format 
                    if [len(m) for m in matrix] == list(range(1, len(self) + 1)):
                        self.matrix = matrix
                    else:
                        raise ValueError(
                            "'matrix' should be in lower triangle format")
                else:
                    raise ValueError(
                        "'names' and 'matrix' should be the same size")
            else:
                raise TypeError("'matrix' should be a list of numerical lists")

    def __getitem__(self, item):
        """Access value(s) by the index(s) or name(s). 
 
        For a _Matrix object 'dm':: 
 
            dm[i]                   get a value list from the given 'i' to others; 
            dm[i, j]                get the value between 'i' and 'j'; 
            dm['name']              map name to index first 
            dm['name1', 'name2']    map name to index first 
        """
        # Handle single indexing 
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
                # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            return [self.matrix[index][i] for i in range(0, index)] + [self.matrix[i][index] for i in
                                                                       range(index, len(self))]
            # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
                # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            if row_index > col_index:
                return self.matrix[row_index][col_index]
            else:
                return self.matrix[col_index][row_index]
        else:
            raise TypeError("Invalid index type.")

    def __setitem__(self, item, value):
        """Set value by the index(s) or name(s). 
 
        Similar to __getitem__:: 
 
            dm[1] = [1, 0, 3, 4]    set values from '1' to others; 
            dm[i, j] = 2            set the value from 'i' to 'j' 
 
        """

        # Handle single indexing 
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
                # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
                # check and assign value
            if isinstance(value, list) and all(_is_numeric(n) for n in value):
                if len(value) == len(self):
                    for i in range(0, index):
                        self.matrix[index][i] = value[i]
                    for i in range(index, len(self)):
                        self.matrix[i][index] = value[i]
                else:
                    raise ValueError("Value not the same size.")
            else:
                raise TypeError("Invalid value type.")
                # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
                # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
                # check and assign value
            if _is_numeric(value):
                if row_index > col_index:
                    self.matrix[row_index][col_index] = value
                else:
                    self.matrix[col_index][row_index] = value
            else:
                raise TypeError("Invalid value type.")
        else:
            raise TypeError("Invalid index type.")

    def __delitem__(self, item):
        """Delete related distances by the index or name"""
        index = None
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self.names.index(item)
        else:
            raise TypeError("Invalid index type.")
            # remove distances related to index
        for i in range(index + 1, len(self)):
            del self.matrix[i][index]
        del self.matrix[index]
        # remove name 
        del self.names[index]

    def insert(self, name, value, index=None):
        """Insert distances given the name and value. 
 
        :Parameters: 
            name : str 
                name of a row/col to be inserted 
            value : list 
                a row/col of values to be inserted 
        """
        if isinstance(name, str):
            # insert at the given index or at the end 
            if index is None:
                index = len(self)
            if not isinstance(index, int):
                raise TypeError("Invalid index type.")
                # insert name
            self.names.insert(index, name)
            # insert elements of 0, to be assigned 
            self.matrix.insert(index, [0] * index)
            for i in range(index, len(self)):
                self.matrix[i].insert(index, 0)
                # assign value
            self[index] = value
        else:
            raise TypeError("Invalid name type.")

    def __len__(self):
        """Matrix length"""
        return len(self.names)

    def __repr__(self):
        return self.__class__.__name__ + "(names=%s, matrix=%s)" % tuple(map(repr, (self.names, self.matrix)))

    def __str__(self):
        """Get a lower triangular matrix string"""
        matrix_string = '\n'.join(
            [self.names[i] + "\t" + "\t".join([str(n) for n in self.matrix[i]])
             for i in range(0, len(self))])
        matrix_string = matrix_string + "\n\t" + "\t".join(self.names)
        return matrix_string


class _DistanceMatrix(_Matrix):
    """Distance matrix class that can be used for distance based tree algorithms. 
 
    All diagonal elements will be zero no matter what the users provide. 
    """

    def __init__(self, names, matrix=None):
        #    print("names: ",names)  #for debugging
        #    print("matrix: ",matrix)  #for debugging
        _Matrix.__init__(self, names, matrix)
        self._set_zero_diagonal()

    def __setitem__(self, item, value):
        _Matrix.__setitem__(self, item, value)
        self._set_zero_diagonal()

    def _set_zero_diagonal(self):
        """set all diagonal elements to zero"""
        for i in range(0, len(self)):
            self.matrix[i][i] = 0


class DistanceCalculator(object):
    """Class to calculate the distance matrix from a DNA or Protein 
 
    Multiple Sequence Alignment(MSA) and the given name of the 
    substitution model. 
 
    Currently only scoring matrices are used. 
 
    :Parameters: 
        model : str 
            Name of the model matrix to be used to calculate distance. 
            The attribute `dna_matrices` contains the available model 
            names for DNA sequences and `protein_matrices` for protein 
            sequences. 
 
    Example 
    ------- 
 
    >> from Bio.Phylo.TreeConstruction import DistanceCalculator 
    >> from Bio import AlignIO 
    >> aln = AlignIO.read(open('Tests/TreeConstruction/msa.phy'), 'phylip') 
    >> print aln 
    SingleLetterAlphabet() alignment with 5 rows and 13 columns 
    AACGTGGCCACAT Alpha 
    AAGGTCGCCACAC Beta 
    GAGATTTCCGCCT Delta 
    GAGATCTCCGCCC Epsilon 
    CAGTTCGCCACAA Gamma 
 
    DNA calculator with 'identity' model:: 
 
        >> calculator = DistanceCalculator('identity') 
        >> dm = calculator.get_distance(aln) 
        >> print dm 
        Alpha   0 
        Beta    0.230769230769  0 
        Gamma   0.384615384615  0.230769230769  0 
        Delta   0.538461538462  0.538461538462  0.538461538462  0 
        Epsilon 0.615384615385  0.384615384615  0.461538461538  0.153846153846  0 
                Alpha           Beta            Gamma           Delta           Epsilon 
 
    Protein calculator with 'blosum62' model:: 
 
        >> calculator = DistanceCalculator('blosum62') 
        >> dm = calculator.get_distance(aln) 
        >> print dm 
        Alpha   0 
        Beta    0.369047619048  0 
        Gamma   0.493975903614  0.25            0 
        Delta   0.585365853659  0.547619047619  0.566265060241  0 
        Epsilon 0.7             0.355555555556  0.488888888889  0.222222222222  0 
                Alpha           Beta            Gamma           Delta           Epsilon 
    """

    dna_alphabet = ['A', 'T', 'C', 'G']

    # BLAST nucleic acid scoring matrix 
    blastn = [[5],
              [-4, 5],
              [-4, -4, 5],
              [-4, -4, -4, 5]]

    # transition/transversion scoring matrix 
    trans = [[6],
             [-5, 6],
             [-5, -1, 6],
             [-1, -5, -5, 6]]

    protein_alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y',
                        'Z']

    # matrices available 
    dna_matrices = {'blastn': blastn, 'trans': trans}
    protein_models = MatrixInfo.available_matrices
    protein_matrices = dict((name, getattr(MatrixInfo, name))
                            for name in protein_models)

    dna_models = list(dna_matrices.keys())

    models = ['identity'] + dna_models + protein_models

    def __init__(self, model='identity'):
        """Initialize with a distance model"""

        if model == 'identity':
            self.scoring_matrix = None
        elif model in self.dna_models:
            self.scoring_matrix = _Matrix(self.dna_alphabet,
                                          self.dna_matrices[model])
        elif model in self.protein_models:
            self.scoring_matrix = self._build_protein_matrix(
                self.protein_matrices[model])
        else:
            raise ValueError("Model not supported. Available models: "
                             + ", ".join(self.models))

    def _pairwise(self, seq1, seq2):
        """Calculate pairwise distance from two sequences"""
        score = 0
        max_score = 0
        if self.scoring_matrix:
            max_score1 = 0
            max_score2 = 0
            skip_letters = ['-', '*']
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 in skip_letters or l2 in skip_letters:
                    continue
                if l1 not in self.scoring_matrix.names:
                    raise ValueError("Bad alphabet '%s' in sequence '%s' at position '%s'"
                                     % (l1, seq1.id, i))
                if l2 not in self.scoring_matrix.names:
                    raise ValueError("Bad alphabet '%s' in sequence '%s' at position '%s'"
                                     % (l2, seq2.id, i))
                max_score1 += self.scoring_matrix[l1, l1]
                max_score2 += self.scoring_matrix[l2, l2]
                score += self.scoring_matrix[l1, l2]

            max_score = max_score1 > max_score2 and max_score1 or max_score2
        else:
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 == l2:
                    score += 1
            max_score = len(seq1)

        return 1 - (score * 1.0 / max_score)

    def get_distance(self, msa):
        """Return a _DistanceMatrix for MSA object 
 
        :Parameters: 
            msa : MultipleSeqAlignment 
                DNA or Protein multiple sequence alignment. 
 
        """

        if not isinstance(msa, MultipleSeqAlignment):
            raise TypeError("Must provide a MultipleSeqAlignment object.")

        names = [s.id for s in msa]
        dm = _DistanceMatrix(names)
        for seq1, seq2 in itertools.combinations(msa, 2):
            dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)
        return dm

    def _build_protein_matrix(self, subsmat):
        """Convert matrix from SubsMat format to _Matrix object"""
        protein_matrix = _Matrix(self.protein_alphabet)
        for k, v in subsmat.items():
            aa1, aa2 = k
            protein_matrix[aa1, aa2] = v
        return protein_matrix


class TreeConstructor(object):
    """Base class for all tree constructor."""

    def build_tree(self, msa):
        """Caller to built the tree from a MultipleSeqAlignment object. 
 
        This should be implemented in subclass"""
        raise NotImplementedError("Method not implemented!")


class DistanceTreeConstructor(TreeConstructor):
    """Distance based tree constructor. 
 
    :Parameters: 
        method : str 
            Distance tree construction method, 'nj'(default) or 'upgma'. 
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
 
    NJ Tree: 
 
    >> njtree = constructor.nj(dm) 
    >> print njtree 
    Tree(rooted=False) 
        Clade(name='Inner3') 
            Clade(branch_length=0.0142054862889, name='Inner2') 
                Clade(branch_length=0.239265540676, name='Inner1') 
                    Clade(branch_length=0.0853101915988, name='Epsilon') 
                    Clade(branch_length=0.136912030623, name='Delta') 
                Clade(branch_length=0.292306275042, name='Alpha') 
            Clade(branch_length=0.0747705106139, name='Beta') 
            Clade(branch_length=0.175229489386, name='Gamma') 
 
 
    """

    methods = ['nj', 'upgma']

    def __init__(self, distance_calculator=None, method="nj"):
        if (distance_calculator is None
                or isinstance(distance_calculator, DistanceCalculator)):
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
            else:
                tree = self.nj(dm)
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

    def nj(self, distance_matrix):
        """Construct and return an Neighbor Joining tree. 
 
        :Parameters: 
            distance_matrix : _DistanceMatrix 
                The distance matrix for tree construction.
            usese Tree_new insted of tree and clade_new insted of clade
        """

        if not isinstance(distance_matrix, _DistanceMatrix):
            raise TypeError("Must provide a _DistanceMatrix object.")

            # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades 
        clades = [CladeNew(None, name) for name in dm.names]
        leaves_lst = []
        # init node distance 
        node_dist = [0] * len(dm)
        # init minimum index 
        min_i = 0
        min_j = 0
        inner_count = 0
        while len(dm) > 2:
            # calculate nodeDist 
            for i in range(0, len(dm)):
                node_dist[i] = 0
                for j in range(0, len(dm)):
                    node_dist[i] += dm[i, j]
                node_dist[i] = node_dist[i] / (len(dm) - 2)

            # find minimum distance pair 
            min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            min_i = 0
            min_j = 1
            for i in range(1, len(dm)):
                for j in range(0, i):
                    temp = dm[i, j] - node_dist[i] - node_dist[j]
                    if min_dist > temp:
                        min_dist = temp
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
            clade1.branch_length = (dm[min_i, min_j] + node_dist[min_i]
                                    - node_dist[min_j]) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            ##updating leaves DS
            if "Inner" not in clade1.name:
                leaves_lst.append(clade1)
            if "Inner" not in clade2.name:
                leaves_lst.append(clade2)

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix, 
            # set the distances of new node at the index of min_j 
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]
                                    - dm[min_i, min_j]) / 2.0

            dm.names[min_j] = "Inner" + str(inner_count)
            del dm[min_i]

            # set the last clade as one of the child of the inner_clade
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
            root = clades[0]
            if "inner_clade" not in clades[1].name:
                leaves_lst.append(clades[1])
            clades[1].set_parent(clades[0])

        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]
            if "inner_clade" not in clades[0].name:
                leaves_lst.append(clades[0])
            clades[0].set_parent(clades[1])

        return Tree_new(root, rooted=False, leaves=leaves_lst)

    def _height_of(self, clade):
        """calculate clade height -- the longest path to any terminal."""
        height = 0
        if clade.is_terminal():
            height = clade.branch_length
        else:
            height = height + max([self._height_of(c) for c in clade.clades])
        return height


class CladeNew(BaseTree.Clade):
    def __init__(self, branch_length=None, name=None, clades=None, confidence=None, colour='w', width=None, parent=None,
                 node_targets=None, candidates=None, polymorphic_sites=None, polymorphic_sites_set=None,
                 lowest_cut_site_prob=1):
        '''
        :param branch_length:
        :param name:
        :param clades:
        :param confidence:
        :param colour:
        :param width:
        :param parent:
        :param node_targets:
        :param candidates:
        :param polymorphic_sites:
        :param lowest_cut_site_prob:
        :param widest: the hiest probability to cleave all the targets by a single sgRNA
        :return:
        '''
        BaseTree.Clade.__init__(self, branch_length, name, clades, confidence, colour, width)
        self.parent = parent
        self.distance_from_leaf = None
        self.node_targets = node_targets or list()
        self.colour = colour
        self.candidates = candidates
        self.lowest_cut_site_prob = lowest_cut_site_prob
        self.set_polymorphic_sites(polymorphic_sites)
        self.set_polymorphic_sites_set(polymorphic_sites_set)

    def set_parent(self, parent):
        self.parent = parent

    def add_node_target(self, leaf):
        self.node_targets.append(leaf)

    def set_candidates_DS(self, candidates_DS=None):
        self.candidates = candidates_DS or dict()

    def set_polymorphic_sites(self, polymorphic_sites=None):
        self.polymorphic_sites = polymorphic_sites or dict()

    def set_polymorphic_sites_set(self, polymorphic_sites_set=None):
        self.polymorphic_sites_set = polymorphic_sites_set or set()


class Tree_new(BaseTree.Tree):
    def __init__(self, root=None, rooted=True, id=None, name=None, leaves=None):
        BaseTree.Tree.__init__(self, root, rooted, id, name)
        self.leaves = leaves
