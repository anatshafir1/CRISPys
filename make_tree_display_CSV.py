
import pickle
import CasSites


def sub_tree_display(path, candidates_lst, f, consider_homology, counter, genes_names, genes_lst):
    def change_mm_to_lowercase(target_str, mm_lst):
        """mm is mismatch"""
        target_in_lst = list(target_str)
        for place in mm_lst:
            target_in_lst[place] = target_in_lst[place].lower()
        return ''.join(target_in_lst)

    header_row = "sgRNA index,sgRNA,Score,Genes,Genes score,Target site,#mms,Position\n"  # new

    if consider_homology:
        # make a set of the genes in the group to write before each table. added by Udi 13/04/22
        genes = []
        for candidate in candidates_lst:
            genes += [g for g in candidate.genes_score_dict.keys()]
        genes = set(genes)
        f.write("Genes in group=," + str(genes) + "\n")

    f.write(header_row)
    sgRNA_index = 0
    for c in candidates_lst:
        sgRNA_index += 1

        num_of_targets = 0
        for targets in c.targets_dict.values():
            num_of_targets += len(targets)
        first_gene = 1
        l = list(c.targets_dict.items())
        l.sort(key=lambda item: c.genes_score_dict[item[0]], reverse=True)

        for gene, targets in l:
            targets.sort(key=lambda target: len(target[1]))
            gene_seq = genes_lst[genes_names.index(gene)]
            seen_sites = dict()
            first_target = 1
            for target in targets:
                pos = find_pos(target, gene_seq, seen_sites)

                if first_target == 1 and first_gene == 1:

                    f.write(str(sgRNA_index) + '.,' + c.seq + "," + str(c.cut_expectation)[:5])
                    f.write("," + gene)
                    score = str(c.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write("," + score)
                    f.write("," + change_mm_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    f.write("," + pos)
                    f.write("\n")
                    first_target = 0
                    continue
                if first_target != 1:
                    f.write(str(sgRNA_index) + ".,,,,,")

                    f.write(change_mm_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    f.write("," + pos)
                    f.write("\n")
                if first_target == 1 and first_gene != 1:
                    f.write(str(sgRNA_index) + ".,,,")
                    score = str(c.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write(gene)
                    f.write("," + score)
                    f.write("," + change_mm_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    f.write("," + pos)
                    f.write("\n")

                    first_target = 0
            first_gene = 0


def find_pos(target, gene_sequence, seen_sites):
    """sgRNA_targets is a list of target sites
	returns"""
    target_seq = target[0]
    if target_seq in seen_sites:
        directions_lst = seen_sites[target_seq]
    else:
        directions_lst = [0, 0]
    position = gene_sequence.find(target_seq, directions_lst[0])
    if position != -1:
        update_seen_sites_dict(seen_sites, target_seq, 0, position)
        position = str(position) + '+'
    else:
        position = CasSites.give_complementary(gene_sequence).find(target_seq, directions_lst[1])
        update_seen_sites_dict(seen_sites, target_seq, 1, position)
        position = str(position) + '-'
    if position == -1:
        position = ''
    return position


def update_seen_sites_dict(d, site_seq, direction, position):
    """
	d: dict: key: site_seq; val: [num of occurrences, directions_lst]
	direction: 0 for sense, 1 for antisense
	"""
    if site_seq in d:
        directions_lst = d[site_seq]
        directions_lst[direction] = position + 1
        d[site_seq] = directions_lst
    else:
        directions_lst = [0, 0]
        directions_lst[direction] = position + 1
        d[site_seq] = directions_lst


def tree_display(path, consider_homology=False, set_cover=False):
    if set_cover:
        candidates_lst = pickle.load(open(path + "/greedy_cover.p", "rb"))
    else:
        candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))
    genes_names = pickle.load(open(path + "/genes_names.p", "rb"))
    genes_list = pickle.load(open(path + '/genes_list.p', 'rb'))

    filepath = path + "/CRISPys_output.csv"
    f = open(filepath, 'w')
    counter = 0
    if not consider_homology:
        sub_tree_display(path, candidates_lst, f, consider_homology, counter, genes_names, genes_list)
    else:
        f.write(
            "The designed sgRNAs for the genes in your input are listed in the table below. Every section of the table corresponds to a homologous genes subgroup as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.\n")

        for subgroup_item in candidates_lst:
            # create the main table
            counter += 1
            sub_tree_display(path, subgroup_item.candidates_list, f, consider_homology, counter, genes_names, genes_list)
    f.close()


if __name__ == "__main__":
    tree_display("/groups/itay_mayrose/galhyams/1516893877")
