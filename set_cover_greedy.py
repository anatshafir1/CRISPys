import copy
import Candidate
import pickle

def prob_cover_genes_lst(candidate, genes_lst):
    cover_all = 1
    for gene in genes_lst:
        cover_all *= candidate.genes_score_dict[gene]
    return cover_all

def find_set_cover(best_permutations_DS, sg_genes_dict, thr, genes_sg_dict = None):
    '''for now, might won't work in a case when there is a gene that isn't covered by any of the permutations in the best_permutations_DS. not finished. can make it more readble'''
    temp_best_perm_DS = copy.copy(best_permutations_DS)
    res = list()#[temp_best_perm_DS[0]]
    if genes_sg_dict:
        for gene, targets in genes_sg_dict.items():
            if len(targets) == 0:
                print("no targets for gene " + gene)
                genes_name_lst.remove(gene)
                continue
            c = Candidate.Candidate(targets[0])
            c.fill_default_fildes(sg_genes_dict[targets[0]])
            temp_best_perm_DS.append(c)

    uncovered_genes = set()
    for sg, genesLst in sg_genes_dict.items():
        for gene in genesLst:
            uncovered_genes.add(gene)
    while(len(uncovered_genes)) > 0 and len(temp_best_perm_DS) > 0:
        #print(uncovered_genes)
        #for gene in uncovered_genes:
            #print(gene)    
    ##going over all the permutations, and return the permutation that cover the maximal amount of genes haven't been covered yet, in the highest probability among the maximal covered permutations
        #print('len uncovered genes', len(uncovered_genes))
        best_current_perm, best_num_of_coverd, best_prob_of_covered = None, 0,0  #best_current_perm is the hole tuple
        i = 0
        while i < (len(temp_best_perm_DS)):
            new_genes_coverd = list()#0
            for gene, score in temp_best_perm_DS[i].genes_score_dict.items():
                if gene in uncovered_genes and score >= thr:
                    new_genes_coverd.append(gene)
                    #uncovered_genes.remove(gene)
            if len(new_genes_coverd) == 0:
                i+=1
                continue
                #del temp_best_perm_DS[i]
            elif len(new_genes_coverd) >= best_num_of_coverd:## and temp_best_perm_DS[i][2] > best_prob_of_covered:  ##need to check if 2 is the right index, and not 1.
                #print(new_genes_coverd)
                if len(new_genes_coverd) > best_num_of_coverd or prob_cover > best_prob_of_covered: # cover more gene or cover the same amount with greater prob.
                    prob_cover = prob_cover_genes_lst(temp_best_perm_DS[i], new_genes_coverd)
				#if prob_cover > best_prob_of_covered:
                    best_num_of_coverd, best_prob_of_covered = len(new_genes_coverd), prob_cover
                    best_current_perm = temp_best_perm_DS[i]
            i+=1
        if(best_current_perm):
            res.append(best_current_perm)
            for gene, score in best_current_perm.genes_score_dict.items():
                if gene in uncovered_genes and score >= thr: #there is a probability that this gene had already been covered bya prevuis sgRNA
                    uncovered_genes.remove(gene)
    return res

def test1(path = '/bioseq/data/results/multicrispr/1516613143' , thr = 0.45):
	best_permutations_DS, sg_genes_dict = pickle.load(open("/".join([path, "res_in_lst.p"]),'rb')), pickle.load(open("/".join([path, "sg_genes_dict.p"]), 'rb'))
	find_set_cover(best_permutations_DS, sg_genes_dict, thr)