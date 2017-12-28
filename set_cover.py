
import copy
import Candidate

def prob_cover_genes_lst(candidate, genes_lst):
    cover_all = 1
    for gene in genes_lst:
        cover_all *= candidate.genes_score_dict[gene]
    return cover_all

def find_set_cover(best_permutations_DS, sg_genes_dict, thr):
    '''for now, might won't work in a case when there is a gene that isn't covered by any of the permutations in the best_permutations_DS. not finished. can make it more readble'''
    temp_best_perm_DS = copy.copy(best_permutations_DS)
    res = list()#[temp_best_perm_DS[0]]
    uncovered_genes = set()
    for sg, genesLst in sg_genes_dict.items():
        for gene in genesLst:
            uncovered_genes.add(gene)
    while(len(uncovered_genes)) > 0 and len(temp_best_perm_DS) > 0:
    ##going over all the permutations, and return the permutation that cover the maximal amount of genes haven't been covered yet, in the highest probability among the maximal covered permutations
        best_current_perm, best_num_of_coverd, best_prob_of_covered = None, 0,0  #best_current_perm is the hole tuple
        i = 0
        while i < (len(temp_best_perm_DS)):
            new_genes_coverd = 0
            for gene, score in temp_best_perm_DS[i].genes_score_dict.items():
                if score < thr:
                    continue
                if gene in uncovered_genes:
                    new_genes_coverd.append(gene)
            if len(new_genes_coverd) == 0:
                del temp_best_perm_DS[i]
            elif len(new_genes_coverd) >= best_num_of_coverd:## and temp_best_perm_DS[i][2] > best_prob_of_covered:  ##need to check if 2 is the right index, and not 1.
                prob_cover = prob_cover_genes_lst(temp_best_perm_DS[i], new_genes_coverd)
                if prob_cover > best_prob_of_covered:
                    best_num_of_coverd, best_prob_of_covered = len(new_genes_coverd), prob_cover
                best_current_perm = temp_best_perm_DS[i]
            i+=1
        if(best_current_perm):
            res.append(best_current_perm)
            for gene, score in best_current_perm.genes_score_dict.items():
                if gene in uncovered_genes and score < thr: #there is a probability that this gene had already been covered bya prevuis sgRNA
                    uncovered_genes.remove(gene)
    return res
