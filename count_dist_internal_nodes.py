from ete3 import Tree
import pickle


NODE_DISTANCE = 3   # define here the desired node distance threshold


###################################################################################################################
"""
The main function to be used. The function removed amiRNAs that are designed to far target genes, and retains those 
whose targets are not too far.
:param dic_amiRNAs: dictionary of candidate AmiRNAs: key = AmiRNA, value = its attributes (should be different for sgRNAs).
:param subgroups: a list of subgroups of genes (representing internal nodes) sorted in a post order.
:param res_subgroups_path: A path to the results of amiRNAs per internal node: key = subgroup of genes (internal node), 
        value = AmiRNAs
:return:A dictionary where key = a subgroup of genes, representing an internal node that was empty before. value: a dictionary
        # of AmiRNAs and their attributes
"""
def removeFarGenesGroups(dic_amiRNAs, subgroups, res_subgroups_path):
    dic_subgroup_AmiRNA = {}    # A dictionary that will be returned  the function. It will contain internal nodes and
    # the new found sgRNAs with their attributes (target genes, score, etc).
    removed = set()
    with open(res_subgroups_path, 'rb') as handle:
        res_subgroups = pickle.load(handle)
    # We want to add AmiRNAs with subgroups that are concordant with the internal nodes in the tree
    # only for internal nodes for which there were not found any AmiRNAs

    #tree_subgroups - all internal nodes
    tree_subgroups = [subgroup for subgroup in res_subgroups if len(res_subgroups[subgroup]) == 0]
    for subgroup_in_tree in tree_subgroups:
        dic_subgroup_AmiRNA[subgroup_in_tree] = {}  # initializing the dictionary with new keys (internal nodes)
    sons_dic, fathers_dic = create_tree(subgroups)  # getting the tree structure as two dictionaries- dictionary of sons
                                                    # for each node, and dictionary for father for each node.
    dic_genes_distance = calculateDistBetweenGenes(sons_dic, fathers_dic, subgroups)    # calculate the number of separating
                                                                                        # nodes between each pair of genes

    subgroups_small_to_large = sorted(tree_subgroups, key= lambda x: len(x)) # sort subgroups of genes (internal nodes) from
                                                                                # smallest to largest
    subgroups_small_to_large = [set(tree_subgroup) for tree_subgroup in subgroups_small_to_large]
    for amiRNA in dic_amiRNAs:
        # find the target genes for a given AmiRNA
        subgroup_amiRNA = getSubgroupOfAmiRNA(dic_amiRNAs, amiRNA) #TODO replace with genes targeted 
        for tree_subgroup in subgroups_small_to_large:
            if subgroup_amiRNA.issubset(tree_subgroup): # finding the common ancestor of subgroup_amiRNA
                # if the target genes of a given AmiRNA are a subset of a particular
                                                        # internal node, it means that this is the closest subgroup of
                                                        # genes which contains this set of target genes
                tree_subgroup_tuple = tuple(sorted(list(tree_subgroup)))
                # cluster the genes to existing subgroups that represent internal nodes of the tree
                genes_clusters = getGeneClusters(tree_subgroup_tuple, subgroup_amiRNA, sons_dic)
                # check if the number of internal nodes separating between the clusters is not larger than three
                near = areAllClustersNear(genes_clusters, dic_genes_distance)
                if not near:
                    # set for removal
                    removed.add(amiRNA)
                    break
                # add the AmiRNA with all its attribute
                dic_subgroup_AmiRNA[tree_subgroup_tuple][amiRNA] = dic_amiRNAs[amiRNA]
                break
                # break the loop, because we are finished with this AmiRNA
    # Create set of the remained candidates
    remained_amiRNAs = set()
    for subgroup in dic_subgroup_AmiRNA:
        for amiRNA in dic_subgroup_AmiRNA[subgroup]:
            remained_amiRNAs.add(amiRNA)

    for miR in removed:
        del dic_amiRNAs[miR]
    return dic_subgroup_AmiRNA
######################################################################################


"""
Creates a structure of subgroups from gene tree, i.e., represents each internal nodes by its descendant leaves.
:param tree_path:  gene tree path.
:param pathForSubgroups: path for subgroups data structure: a list groups of genes sorted in postorder (each element in 
                        the list is a set of genes).
:return: No return value- the subgroups data structure is saved in a pickle file specified by @param2.
"""
def createSubgroupsOfGenes(tree_path, pathForSubgroups):
    setOfGenes = []
    tree = Tree(tree_path)
    dicNodesLeaves = tree.get_cached_content()
    for node in tree.traverse("postorder"):
        currSet = dicNodesLeaves[node]
        setOfNames = set()
        for n in currSet:
            setOfNames.add(n.name)
        if (len(setOfNames) == 1):
            continue
        setOfGenes.append(setOfNames)
    with open(pathForSubgroups, 'wb') as handle:
        pickle.dump(setOfGenes, handle)
    return

############################################################################################################
"""
The function gets as input a list of internal nodes (subgroups of genes), and returns dictionaries of sons and fathers
:param subgroups_not_ordered: unsorted list of subgroups.
:return: 1) sons is a dictionary where key = subgroup of genes (internal node), 
            and value = a list of direct sons' subgroups
         2) fathers is a dictionary where key = subgroup of genes (internal node), and value = subgroup of genes
            which represets the father internal node
"""
def create_tree(subgroups_not_ordered):
    subgroups = sorted(subgroups_not_ordered, key = lambda x: len(x))
    sons = {}
    fathers = {}
    for i in range(len(subgroups)-1):
        minor_group = tuple(sorted(list(subgroups[i])))
        if len(subgroups[i]) == 2:
            singeltons = [{sub} for sub in subgroups[i]]
            sons[minor_group] = [tuple(singlton) for singlton in singeltons]
            for singlton in singeltons:
                fathers[tuple(sorted(list(singlton)))] = minor_group
        for j in range(i+1, len(subgroups)):
            major_group = tuple(sorted(list(subgroups[j])))

            if subgroups[i].issubset(subgroups[j]):
                if major_group in sons:
                    sons[major_group].append(minor_group)
                else:
                    sons[major_group] = [minor_group]
                fathers[minor_group] = major_group
                break
    for node in sons:
        if len(sons[node]) < 2:
            son = set(list(sons[node][0]))
            node_set = set(list(node))
            gene = node_set.difference(son)
            sons[node].append(tuple(list(gene)))
            fathers[tuple(list(gene))] = node

    return sons, fathers
#######################################################################################################################
"""
gets the corresponding subgroup of genes for a given AmiRNA
:param dic_amiRNAs: dictionary of AmiRNAs, where key = AmiRNA, and value = AmiRNAs attributes.
:param amiRNA
:return: a set of target genes.
"""
def getSubgroupOfAmiRNA(dic_amiRNAs, amiRNA):
    subgroup = [target[0] for target in dic_amiRNAs[amiRNA]["targets"]] + [target[0] for target in
                                                                                 dic_amiRNAs[amiRNA][
                                                                                     "belowThresholdTargets"]]
    subgroup = set(subgroup)
    return subgroup

#######################################################################################################################
"""
For each pair of genes get the distance between the two genes (how many nodes separate between them)
:param sons_dic: dictionary of sons for each node in the tree: key = a subgroup of genes representing an internal node,
                value = list of subgroups of genes, representing the internal nodes of the direct node's sons
:param father_dic: dictionary of fathers- where key = subgroup of genes (internal node), and value = subgroup of genes
                    which represets the father internal node.
:param subgroups:  subgroups of genes sorted in a post order.
:return: dictionary, where key = tuple of pair of genes, value = the number of nodes separating between them
"""
def calculateDistBetweenGenes(sons_dic, father_dic, subgroups):
    dic_dist = {}
    genes = list(subgroups[-1]) # the last element refers to the root, so that all genes are contained within it
    for i in range(len(genes)-1):
        gene1 = genes[i]
        for j in range(i+1, len(genes)):
            subset = {gene1}
            subset_tuple = tuple([gene1])
            gene2 = genes[j]
            distance = 0
            while not (gene2 in subset):
                subset_tuple = father_dic[subset_tuple]
                subset = set(list(subset_tuple))
                distance += 1
            while subset != {gene2}:
                son1 = sons_dic[subset_tuple][0]
                setSon1 = set(list(son1))
                son2 = sons_dic[subset_tuple][1]
                setSon2 = set(list(son2))
                if gene2 in setSon1:
                    subset = setSon1
                    subset_tuple = son1
                elif gene2 in setSon2:
                    subset = setSon2
                    subset_tuple = son2
                else:
                    raise Exception("Error in calculate distance!\n")
                distance += 1

            dic_dist[tuple(sorted([gene1, gene2]))] = distance-1
    return dic_dist

###################################################################################################################
"""
make a partition to monophyletic targets
:param tree_subgroup_tuple: internal node in a tree (subgroup of genes)
:param subgroup_amiRNA: A given subgroup of genes (which does not represent an internal node in a tree)
:param sons_dic: A dictionary of sons, where key = subgroup of genes (internal node), and value = a list of direct 
                    sons' subgroups
:return: list of clusters of genes that are clustered according to subtrees in a tree.
"""
def getGeneClusters(tree_subgroup_tuple, subgroup_amiRNA, sons_dic):
    lst_of_clusters = []
    initial_subgroup = subgroup_amiRNA
    son1 = sons_dic[tree_subgroup_tuple][0]
    son2 = sons_dic[tree_subgroup_tuple][1]

    getGeneClustersRec(lst_of_clusters, initial_subgroup, son1, son2, sons_dic)
    return lst_of_clusters
###################################################################################################################
"""
Auxilary function for getGeneClusters. A recursive function which searches for the required clusters of subgroups of genes
:param lst_of_clusters: list of formed clusters
:param initial_subgroup: the subgroup that should be further partitioned
:param son1
:param son2
:param sons_dic
:return: No return value.
"""
def getGeneClustersRec(lst_of_clusters, initial_subgroup, son1, son2, sons_dic):
    monophyletic_subgrpups = [set(), set()]
    for gene in initial_subgroup:
        if gene in son1:
            monophyletic_subgrpups[0].add(gene)
        elif gene in son2:
            monophyletic_subgrpups[1].add(gene)
        else:
            print("gene", gene)
            print("son1:", son1)
            print("son2:", son2)
            raise Exception("ERROR!!! getGeneClustersRec(): gene is not found in any subgroup!")

    if len(initial_subgroup) == 0:
        return
    sons = [son1, son2]
    for i in range(2):
        if monophyletic_subgrpups[i] == set(sons[i]):
            lst_of_clusters.append(monophyletic_subgrpups[i])
        else:
            if len(monophyletic_subgrpups[i]) != 0:
                grand_son1 = sons_dic[sons[i]][0]
                grand_son2 = sons_dic[sons[i]][1]
                getGeneClustersRec(lst_of_clusters, monophyletic_subgrpups[i], grand_son1, grand_son2, sons_dic)


######################################################################################################################
"""
Test whether all clusters are not too far from each other
:param genes_clusters: list of clusters of genes that are clustered according to the tree.
:param dic_genes_distance: Dictionary, where key = tuple of pair of genes, value = the number of nodes separating them.
:return: True if the genes are near (at most three nodes apart). False- otherwise.
"""
def areAllClustersNear(genes_clusters, dic_genes_distance):
    for i in range(len(genes_clusters)-1):
        for j in range(i+1, len(genes_clusters)):
            near = False
            cluster1 = genes_clusters[i]
            cluster2 = genes_clusters[j]
            for gene_cluster1 in cluster1:
                for gene_cluster2 in cluster2:
                    if dic_genes_distance[tuple(sorted([gene_cluster1, gene_cluster2]))] <= NODE_DISTANCE:
                        near = True
                        break
                if near:
                    break
            if not near:
                return False
    return True