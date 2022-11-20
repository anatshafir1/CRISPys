
import set_cover_greedy
import SubgroupRes
from typing import List, Dict


def sub_tree_display(candidates_lst, f, consider_homology):

    header_row = "sgRNA index,sgRNA,Score,Genes,Genes score,Target site,#mms,Position,strand,PAM\n"

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
            seen_sites = dict()
            first_target = 1
            for target in targets:
                 # pos = make_tree_display.find_pos(target, gene_seq, seen_sites) # this might be needed if running with the crispr website

                if first_target == 1 and first_gene == 1:

                    f.write(str(sgRNA_index) + '.,' + candidate.seq + "," + str(candidate.cut_expectation)[:5])
                    f.write("," + gene)
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write("," + score)
                    f.write("," + change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    # f.write("," + pos)
                    #write position (added by Udi 03/11/2022)
                    f.write("," + str(target[3]))
                    # write strand
                    f.write("," + target[4])
                    # write pam
                    f.write("," + target[2])
                    f.write("\n")
                    first_target = 0
                    continue
                if first_target != 1:
                    f.write(str(sgRNA_index) + ".,,,,,")

                    f.write(change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    # f.write("," + pos)
                    # write position
                    f.write("," + str(target[3]))
                    # write strand
                    f.write("," + target[4])
                    # write pam
                    f.write("," + target[2])
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
                    # f.write("," + pos)
                    # write position
                    f.write("," + str(target[3]))
                    # write strand
                    f.write("," + target[4])
                    # write pam
                    f.write("," + target[2])
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


def tree_display(path: str, subgroups_lst: list, genes_list: list, targets_genes_dict: dict,
                 omega: float, set_cover: bool = False, consider_homology: bool = False):
    """
    This function takes the results of crispys and write the crispys results in a CSV output to the output folder
    Args:
        path: path to output folder
        subgroups_lst: list os subgroup objects
        genes_list: a list of gene sequence
        targets_genes_dict: a dictionary with the all targets and the genes they capture
        omega: the value of the score threshold
        set_cover: boolean, if True the output will use the set cover function which output a minimum number of guides that will target all genes
        consider_homology: boolean, if the results were made with the 'consider homology' option

    Returns:

    """
    if set_cover:
        res = set_cover_greedy.find_set_cover(subgroups_lst, targets_genes_dict, omega)
        subgroups_lst = [SubgroupRes.SubgroupRes(genes_list, res, "set_group")]

    filepath = path + "/CRISPys_output.csv"
    f = open(filepath, 'w')
    f.write(
        "The designed sgRNAs for the genes in your input are listed in the table below. Every section of the table corresponds to a homologous genes subgroup as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.\n")

    for subgroup_item in subgroups_lst:
        # create the main table
        sub_tree_display(subgroup_item.candidates_list, f, consider_homology)
    f.close()


def create_output_multiplex(path: str, crispys_res: List, multiplex_dict: Dict, number_of_groups: int,
                            n_with_best_guide: int, n_sgrnas: int):
    """
    This function is used to write the output of multiplex
    Args:
        path: path to output folder
        crispys_res: the 'traditional' crispys output (list of subGroupRes)
        multiplex_dict: the output of crispys-chips. a dictionary of node_name:dictionary of best_sg_seq:BestSgGroup

    Returns:

    """
    filepath = path + "/CRISPys_output_multiplex.csv"
    f = open(filepath, 'w')
    f.write(f"Run multiplex with {number_of_groups} 'Best' groups each one with {n_with_best_guide} "
            f"gRNA and {n_sgrnas} in each multiplex\n")
    # go over each internal node
    for node in multiplex_dict:
        f.write(f"Node:,{node},")
        # get node genes
        for subgroup in crispys_res:
            if subgroup.name == node:
                genes_lst = subgroup.genes_lst
                # write the node name
                f.write(f"genes in node:,{str(genes_lst).strip('[]')}\n")
        # go over each 'best guide' group
        for bestseq in multiplex_dict[node].values():
            # write the sequence of best guide
            f.write(f"Group of:,{bestseq.best_candidate.seq}\n")
            # go over each pair (or more) of multiplex and write it to the file
            for subgroup in bestseq.subgroups:
                sub_tree_display(subgroup.candidates_list, f, consider_homology=False)
        f.write("\n")
    f.close()

if __name__ == "__main__":
    tree_display("/groups/itay_mayrose/galhyams/1516893877")
