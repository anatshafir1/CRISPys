from time import perf_counter

def create_gene2distance_dict(node):
    """
    This function creates a dictionary that stores the distances between each pair of genes.
    Note that the distance between genes is defined as the number of internal nodes between the genes.
    :param node: A tree node
    :return: A dictionary (genei,genej) -> the distance between genei and gene j
    """
    gene_pair2distance_dict = {}
    gene2node_dict = {}  # Create a dictionary that stores the node
    for i, gene_i in enumerate(node.node_leaves):
        for gene_j in node.node_leaves[i + 1:]:
            if gene_i == gene_j or tuple(sorted([gene_i, gene_j])) in gene_pair2distance_dict:
                continue  # Continue if the genes are identical or if the distance was already calculated.

            dist_i = find_distance_from_lca(gene_i, gene_j, node,
                                            gene2node_dict)  # Find the distance of each gene from the least common ancestor
            dist_j = find_distance_from_lca(gene_j, gene_i, node, gene2node_dict)
            gene_pair2distance_dict[tuple(sorted([gene_i, gene_j]))] = dist_i + dist_j - 1
    return gene_pair2distance_dict


def find_distance_from_lca(gene_i, gene_j, node, gene2node_dict):
    """
    This function finds the distance of gene i from the least common ancestor of genes i and j
    :param gene2node_dict: A dictionary
    :param node: A tree node
    :param gene_i:
    :param gene_j:
    :return: the distance of gene i from the least common ancestor of genes i and j
    """
    node_i = find_leaf(node, gene_i, gene2node_dict)  # Find the leaf node containing the gene
    dist_i = 0
    while gene_j not in set(node_i.node_leaves):  # From each node, go to the least common-ancestor of genes i & j
        node_i = node_i.parent
        dist_i += 1
    return dist_i  # Return the distance


def find_leaf(node, gene, gene2node_dict):
    """
    This function finds the leaf node of a gene within a subtree.
    :param node: A tree node.
    :param gene: The name of the gene
    :param gene2node_dict: A dictionary of gene -> node
    :return: A CladeNew object of the node with the single gene. If it's not in gene2dist_from_root_dict, it is updated
    """
    if gene in gene2node_dict:  # If the gene is already in the dictionary, return the value.
        return gene2node_dict[gene]
    gene2node_dict[gene] = node
    while gene2node_dict[gene].node_leaves != [gene]:  # Go down in the tree until the gene-node is reached.
        if gene in gene2node_dict[gene].clades[0].node_leaves:
            gene2node_dict[gene] = gene2node_dict[gene].clades[0]
        else:
            gene2node_dict[gene] = gene2node_dict[gene].clades[1]
    return gene2node_dict[gene]


def get_gene_clusters(genes_targeted: set, node, list_of_clusters):
    """
    This recursive function divides the genes targeted by the sgRNA into monophyletic groups
    :param genes_targeted: The genes targeted by the sgRNA
    :param node: A tree node
    :param list_of_clusters: A list of monophyletic clusters
    :return: None
    """
    node_leaves_set = set(node.node_leaves)
    # If the genes targeted are the same as genes in the node, append the monophyletic group to the list
    if genes_targeted == node_leaves_set:
        list_of_clusters.append(node_leaves_set)
        return
    # If not, it means that the clusters are deeper in the tree.
    intersection_left = genes_targeted.intersection(set(node.clades[0].node_leaves))
    intersection_right = genes_targeted.intersection(set(node.clades[1].node_leaves))
    # If at least one gene appears in both the targeted genes set and the left/right clade, recursively call the function,
    # with the left/right child as the new node, and the intersection as the new set of genes targeted
    if intersection_left:
        get_gene_clusters(intersection_left, node.clades[0], list_of_clusters)
    if intersection_right:
        get_gene_clusters(intersection_right, node.clades[1], list_of_clusters)
    return


def decide_if_near(list_of_clusters, gene_pair2distance_dict, max_gap_distance=3):
    """
    This function determines if the gene clusters targeted by an sgRNA are close enough
    :param list_of_clusters: A list of monophyletic gene sets
    :param gene_pair2distance_dict: A dictionary (gene_i,gene_j) -> the number of internal nodes between them
    :param max_gap_distance: The maximal distance that is allowed between the genes targeted by the sgRNA
    :return: False if at least two monophyletic sets are too distant in the tree, else True
    """
    for i, cluster_i in enumerate(list_of_clusters):
        for cluster_j in list_of_clusters[i + 1:]:
            near = False
            for gene_i in cluster_i:
                for gene_j in cluster_j:
                    # If at least two genes from each cluster have a distance smaller or equal to the threshold,
                    # The two clusters will be considered as close enough
                    if gene_pair2distance_dict[tuple(sorted([gene_i, gene_j]))] <= max_gap_distance:
                        near = True
                        break
                if near:
                    break
            # If at least one pair of clusters is too spaced, the genes targeted by the sgRNA are not considered near.
            if not near:
                return False
    return True


def remove_candidates_with_distant_genes(node, list_of_candidates, max_gap_distance=3):
    """
    This function filters out genes that target genes that are too distanct in the tree.
    :param max_gap_distance: The maximal distance that is allowed between the genes targeted by the sgRNA
    :param node: The internal node in the gene tree
    :param list_of_candidates: A list of sgRNA Candidate objects
    :return:
    """
    # Calculate the distance between all genes in the internal node.
    # Note that the distance between genes is defined as the number of internal nodes between the genes
    gene_pair2distance_dict = create_gene2distance_dict(node)
    filtered_candidates_list = []
    for candidate in list_of_candidates:
        genes_targeted = set(candidate.genes_score_dict.keys())

        list_of_clusters = []
        get_gene_clusters(genes_targeted, node,
                          list_of_clusters)  # Divide the genes targeted by the sgRNA into monophyletic groups.

        is_near = decide_if_near(list_of_clusters, gene_pair2distance_dict,
                                 max_gap_distance)  # Determine if all clusters

        if is_near:
            filtered_candidates_list.append(candidate)
    return filtered_candidates_list
