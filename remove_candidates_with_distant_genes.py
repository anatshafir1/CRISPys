
def get_gene_cluster_mrca_nodes(genes_targeted: set, node, list_of_cluster_nodes):
    """
    This recursive function divides the genes targeted by the sgRNA into tree nodes containing only targeted genes
    :param genes_targeted: The genes targeted by the sgRNA
    :param node: A tree node
    :param list_of_cluster_nodes: A list of nodes, where each node is the most recent common ancestor (mrca) of targeted
    genes.
    :return: None
    """
    node_leaves_set = set(node.node_leaves)
    # If the genes targeted by the sgRNA are the same as the genes in the node, append the node to the list
    if genes_targeted == node_leaves_set:
        list_of_cluster_nodes.append(node)
        return
    # else, call recursively on the 'left' and 'right' child, if they contain any targeted genes
    genes_targeted_in_left_child = genes_targeted.intersection(set(node.clades[0].node_leaves))
    genes_targeted_in_right_child = genes_targeted.intersection(set(node.clades[1].node_leaves))
    if genes_targeted_in_left_child:
        get_gene_cluster_mrca_nodes(genes_targeted_in_left_child, node.clades[0], list_of_cluster_nodes)
    if genes_targeted_in_right_child:
        get_gene_cluster_mrca_nodes(genes_targeted_in_right_child, node.clades[1], list_of_cluster_nodes)
    return


def are_all_clusters_near(list_of_cluster_nodes, max_gap_distance=3):
    """
    This function determines if all targeted gene clusters are close enough. First, the distance, defined as the
    number of branches between the clusters, is calculated for each pair of clusters. Then, for each cluster, the min
    distance is compared to max_gap_distance. If the minimal distance of ALL clusters is lower than max_gap_distance,
    all gene clusters are near and the function will return True. If at least one cluster is too distant from all
    other clusters, return False.
    :param list_of_cluster_nodes: A list of nodes, where each node is the most recent common ancestor (mrca) of targeted
    :param max_gap_distance: The maximal distance that is allowed between the internal nodes of the genes
     targeted by the sgRNA
    :return: True if the gene clusters are close enough, False otherwise
    """
    cluster_pair2distance_dict = {}
    for i, node_i in enumerate(list_of_cluster_nodes):
        min_cluster_distance = -1
        for j, node_j in enumerate(list_of_cluster_nodes):
            if i == j:  # skip if node_i is the same as node_j
                continue
            node_names_pair = tuple(sorted([node_i.name, node_j.name]))
            cluster_distance = get_cluster_distance(node_names_pair, node_i, node_j, cluster_pair2distance_dict)
            if cluster_distance < min_cluster_distance or min_cluster_distance == -1:  # Update the min cluster distance
                min_cluster_distance = cluster_distance
        if min_cluster_distance > max_gap_distance:
            return False  # Return False if the min cluster distance of at least one cluster is above the threshold
    return True


def get_cluster_distance(node_names_pair, node_i, node_j, cluster_pair2distance_dict):
    """
    This function returns the distance between two clusters.
    :param node_names_pair: a tuple of two cluster names
    :param node_i: the node (most recent common ancestor) of cluster i
    :param node_j: the node (most recent common ancestor) of cluster j
    :param cluster_pair2distance_dict: a dictionary of node_names_pair -> the distance between the two clusters
    :return: the distance between the most recent common ancestors of clusters i and j.
    """
    if node_names_pair not in cluster_pair2distance_dict:
        dist_i = find_distance_from_mrca(node_i, node_j)
        dist_j = find_distance_from_mrca(node_j, node_i)
        cluster_pair2distance_dict[node_names_pair] = dist_i + dist_j  # calculate the distance between the clusters
    return cluster_pair2distance_dict[node_names_pair]


def find_distance_from_mrca(node_i, node_j):
    """
    This function finds the distance of gene i from the most common ancestor of clusters i and j
    :param node_i: the node (most recent common ancestor) of cluster i
    :param node_j: the node (most recent common ancestor) of cluster j
    :return: the distance of node_i from the least common ancestor of nodes i and j
    """
    dist_i = 0
    while not set(node_j.node_leaves).intersection(set(node_i.node_leaves)):  # Go up node i until the mrca is reached
        node_i = node_i.parent
        dist_i += 1
    return dist_i


def remove_candidates_with_distant_genes(node, list_of_candidates, max_gap_distance=3):
    """
    This function filters out sgRNAs that target genes that are too distant in the genes tree. :param
    sgRNA.
    Do note that this rule may not suffice for gene trees larger than 8.
    :param max_gap_distance: The maximal distance that is allowed between the internal nodes of the genes targeted by the
    :param node: The internal node in the gene tree
    :param list_of_candidates: A list of sgRNA Candidate
    objects
    :return: A filtered list of sgRNA candidates, where all Candidate object target genes that are close
    enough in the tree.
    """
    filtered_candidates_list = []
    for candidate in list_of_candidates:
        genes_targeted = set(candidate.genes_score_dict.keys())
        list_of_cluster_nodes = []
        # Divide the genes targeted by the sgRNA into cluster tree nodes.
        get_gene_cluster_mrca_nodes(genes_targeted, node,
                                    list_of_cluster_nodes)
        is_near = are_all_clusters_near(list_of_cluster_nodes, max_gap_distance)  # Determine if all clusters are near
        if is_near:
            filtered_candidates_list.append(candidate)
    return filtered_candidates_list
