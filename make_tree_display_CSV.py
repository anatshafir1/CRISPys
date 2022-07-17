
import pickle
import make_tree_display


def sub_tree_display(path, candidates_lst, f, consider_homology, counter, genes_names, genes_lst):

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
    for candidate in candidates_lst:
        sgRNA_index += 1
        num_of_targets = 0
        for targets in candidate.targets_dict.values():
            num_of_targets += len(targets)
        first_gene = 1
        l = list(candidate.targets_dict.items())
        l.sort(key=lambda item: candidate.genes_score_dict[item[0]], reverse=True)

        for gene, targets in l:
            targets.sort(key=lambda target: len(target[1]))
            gene_seq = genes_lst[genes_names.index(gene)]
            seen_sites = dict()
            first_target = 1
            for target in targets:
                pos = make_tree_display.find_pos(target, gene_seq, seen_sites)

                if first_target == 1 and first_gene == 1:

                    f.write(str(sgRNA_index) + '.,' + candidate.seq + "," + str(candidate.cut_expectation)[:5])
                    f.write("," + gene)
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write("," + score)
                    f.write("," + change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    f.write("," + pos)
                    f.write("\n")
                    first_target = 0
                    continue
                if first_target != 1:
                    f.write(str(sgRNA_index) + ".,,,,,")

                    f.write(change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    f.write("," + pos)
                    f.write("\n")
                if first_target == 1 and first_gene != 1:
                    f.write(str(sgRNA_index) + ".,,,")
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write(gene)
                    f.write("," + score)
                    f.write("," + change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    f.write("," + pos)
                    f.write("\n")

                    first_target = 0
            first_gene = 0


def change_mismatch_to_lowercase(target_str, mm_lst):
    """

    :param target_str:
    :param mm_lst:
    :return:
    """
    target_in_lst = list(target_str)
    for place in mm_lst:
        target_in_lst[place] = target_in_lst[place].lower()
    return ''.join(target_in_lst)


def tree_display(path, consider_homology=False, set_cover=False):
    """

    :param path:
    :param consider_homology:
    :param set_cover:
    """
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
