import re
from itertools import combinations
from typing import Tuple, Dict, List

from Amplicon_construction.FindTargets import give_complementary, valid_distance_for_multiplex
from Amplicon_construction.Target_Obj import Family_Target_Obj, FamilyMultiplexTarget
from CRISPys_master.Candidate import Candidate
from CRISPys_master.SubgroupRes import SubgroupRes
from globals import SCORING_DELTA, ALLELES_EXPONENT


def get_cut_alleles_dict(targets_dict):
    sgrna_gene_cut_alleles_dict = {}
    for seq_id in targets_dict:
        gene_name = seq_id.split("_")[0]
        allele1_name = seq_id.split("_")[1]
        if gene_name not in sgrna_gene_cut_alleles_dict:
            sgrna_gene_cut_alleles_dict[gene_name] = [allele1_name]
        else:
            sgrna_gene_cut_alleles_dict[gene_name].append(allele1_name)
    return sgrna_gene_cut_alleles_dict


def calc_family_multiplex_score(sgrna_pair, genes_alleles_dict, family_targeting: int):
    sgrnas_scores_dict = {gene: {allele: 1 for allele in genes_alleles_dict[gene]} for gene in genes_alleles_dict}
    genes_num = 0
    alleles_num = 0
    sgrna_cut_genes_alleles_dict = {}
    if family_targeting == 1:
        sgrna_pair = [sgrna_pair]
    for sgrna in sgrna_pair:
        cut_genes_alleles_dict = {}
        for seq_id in sgrna.genes_score_dict:
            gene_name = seq_id.split("_")[0]
            allele_name = seq_id.split("_")[1]
            if gene_name not in cut_genes_alleles_dict:
                cut_genes_alleles_dict[gene_name] = [allele_name]
            else:
                cut_genes_alleles_dict[gene_name].append(allele_name)
            if sgrnas_scores_dict[gene_name][allele_name] == 1:
                sgrnas_scores_dict[gene_name][allele_name] = (1 - SCORING_DELTA * sgrna.genes_score_dict[seq_id])
            else:
                sgrnas_scores_dict[gene_name][allele_name] *= (1 - SCORING_DELTA * sgrna.genes_score_dict[seq_id])
        sgrna_cut_genes_alleles_dict[sgrna.seq] = cut_genes_alleles_dict
    multiplex_score = 0
    for gene in sgrnas_scores_dict:
        alleles_sum = 0
        for allele in sgrnas_scores_dict[gene]:
            alleles_sum += (1 - sgrnas_scores_dict[gene][allele])
            if sgrnas_scores_dict[gene][allele] > 0:
                alleles_num += 1
        multiplex_score += alleles_sum ** ALLELES_EXPONENT
        if alleles_sum > 0:
            genes_num += 1
    return alleles_num, genes_num, multiplex_score, sgrna_cut_genes_alleles_dict


def create_family_multiplex_target(target1: Family_Target_Obj, target2: Family_Target_Obj,
                                   exon_num: int, rank: int) -> FamilyMultiplexTarget:
    if target1.start_idx < target2.start_idx:
        upstream_target = target1
        downstream_target = target2
    else:
        upstream_target = target2
        downstream_target = target1
    family_multiplex_target = FamilyMultiplexTarget(upstream_target.start_idx, upstream_target.end_idx,
                                                    upstream_target.seq,
                                                    downstream_target.start_idx, downstream_target.end_idx,
                                                    downstream_target.seq,
                                                    0, upstream_target.cut_alleles_lst,
                                                    downstream_target.cut_alleles_lst,
                                                    exon_num, rank, upstream_target.strand, downstream_target.strand,
                                                    upstream_target.sgrna,
                                                    downstream_target.sgrna)
    return family_multiplex_target


def create_genes_single_targets_dict_permutations(genes_targets_dict: Dict[str, Dict[int, List[Family_Target_Obj]]]):
    updated_genes_targets_dict = {}
    for gene in genes_targets_dict:
        exon_targets_dicts_list = []
        for curr_exon in genes_targets_dict[gene]:
            for target in genes_targets_dict[gene][curr_exon]:
                curr_target_dict = {exon: ([] if exon != curr_exon else [target]) for exon in genes_targets_dict[gene]}
                exon_targets_dicts_list.append(curr_target_dict)
        updated_genes_targets_dict[gene] = exon_targets_dicts_list
    return updated_genes_targets_dict


def create_genes_multiplex_targets_dict_permutations(genes_targets_dict: Dict[str, Dict[int, List[Family_Target_Obj]]],
                                                     rank: int, max_amplicon_len: int, primer_length: int,
                                                     target_surrounding_region: int):
    updated_genes_targets_dict = {}
    for gene in genes_targets_dict:
        exon_targets_dicts_list = []
        exons_completed_dict = {exon: False for exon in genes_targets_dict[gene]}
        for exon1 in genes_targets_dict[gene]:
            for target1 in genes_targets_dict[gene][exon1]:
                for exon2 in genes_targets_dict[gene]:
                    if not exons_completed_dict[exon2]:
                        for target2 in genes_targets_dict[gene][exon2]:
                            if target1.seq != target2.seq:
                                if exon1 == exon2:
                                    if valid_distance_for_multiplex(target1, target2, max_amplicon_len, primer_length, target_surrounding_region):
                                        curr_target_dict = {exon: ([] if exon != exon1 else [
                                            create_family_multiplex_target(target1, target2, exon1, rank)]) for exon in
                                                            genes_targets_dict[gene]}
                                        curr_target_dict["multiplex"] = 1
                                        exon_targets_dicts_list.append(curr_target_dict)
                                    else:
                                        curr_target_dict = {exon: [] for exon in genes_targets_dict[gene]}
                                        target1.rank = rank + 0.1
                                        target2.rank = rank + 0.2
                                        curr_target_dict[exon1] = [target1, target2]
                                        curr_target_dict["multiplex"] = 0
                                        exon_targets_dicts_list.append(curr_target_dict)
                                else:
                                    curr_target_dict = {exon: [] for exon in genes_targets_dict[gene]}
                                    target1.rank = rank + 0.1
                                    target2.rank = rank + 0.2
                                    curr_target_dict[exon1] = [target1]
                                    curr_target_dict[exon2] = [target2]
                                    curr_target_dict["multiplex"] = 0
                                    exon_targets_dicts_list.append(curr_target_dict)

            exons_completed_dict[exon1] = True
        if exon_targets_dicts_list:
            updated_genes_targets_dict[gene] = exon_targets_dicts_list
        else:
            return {}
    return updated_genes_targets_dict


def get_single_sg_targets(genes_targets_dict, sgrna, genes_exons_seq_dict, intron_region_added, failed_targets):
    processed_genes_dict = {gene: False for gene in genes_exons_seq_dict}
    sgrna_gene_cut_alleles_dict = get_cut_alleles_dict(sgrna.targets_dict)
    for seq_id in sgrna.targets_dict:
        gene_name = seq_id.split("_")[0]
        if processed_genes_dict[gene_name]:
            continue
        allele1_name = seq_id.split("_")[1]
        target_seq = sgrna.targets_dict[seq_id][0][0] + sgrna.targets_dict[seq_id][0][2]  # target and PAM
        target_strand = sgrna.targets_dict[seq_id][0][4]
        gene_exons_seq_dict = genes_exons_seq_dict[gene_name][0]
        if gene_name not in genes_targets_dict:
            genes_targets_dict[gene_name] = {}
        for exon in gene_exons_seq_dict:
            if exon not in genes_targets_dict[gene_name]:
                genes_targets_dict[gene_name][exon] = []
            target = None
            target_found = False
            for allele in gene_exons_seq_dict[exon]:
                allele2_name = allele[0].split(":")[0][1:].upper()
                if allele1_name == allele2_name:
                    if target_strand == "+":
                        matches = list(re.finditer(target_seq, allele[1]))
                        if matches:
                            seq = matches[0]
                            target = Family_Target_Obj(seq.group(0), intron_region_added + seq.start(),
                                                       intron_region_added + seq.end() - 1, "+", allele2_name,
                                                       seq.group(0),
                                                       0, 0, sgrna_gene_cut_alleles_dict[gene_name], exon)
                            target_found = True
                            break
                    else:  # target_strand == "-"
                        complementary_target_seq = give_complementary(target_seq)
                        matches = list(re.finditer(complementary_target_seq, allele[1]))
                        if matches:
                            seq = matches[0]
                            target = Family_Target_Obj(seq.group(0), intron_region_added + seq.start(),
                                                       intron_region_added + seq.end() - 1, "-", allele2_name,
                                                       seq.group(0), 0, 0, sgrna_gene_cut_alleles_dict[gene_name],
                                                       exon)
                            target_found = True
                            break
            if target_found:
                if target.seq not in failed_targets:
                    genes_targets_dict[gene_name][exon].append(target)
                else:
                    if target not in failed_targets[target.seq]:
                        genes_targets_dict[gene_name][exon].append(target)

        processed_genes_dict[gene_name] = True


def get_genes_single_targets_dict(sgrna: Candidate, genes_exons_seq_dict: Dict[
    str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]], max_amplicon_len: int, primer_length: int,
                                  cut_location: int, target_surrounding_region: int):
    genes_targets_dict = {}
    intron_region_added = max_amplicon_len - primer_length - cut_location - target_surrounding_region  # 255 by default
    failed_targets = {}
    get_single_sg_targets(genes_targets_dict, sgrna, genes_exons_seq_dict, intron_region_added, failed_targets)
    updated_genes_targets_dict = create_genes_single_targets_dict_permutations(genes_targets_dict)
    return updated_genes_targets_dict


def get_genes_multiplex_targets_dict(sgrna_pair: Tuple[Candidate, Candidate], genes_exons_seq_dict: Dict[
                                       str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]],
                                     max_amplicon_len: int, primer_length: int, cut_location: int,
                                     target_surrounding_region: int, failed_targets: Dict[str, List], rank: int):
    genes_targets_dict = {}  # dictionary of {gene: {exon: targets list}}
    intron_region_added = max_amplicon_len - primer_length - cut_location - target_surrounding_region  # 255 by default
    for sgrna in sgrna_pair:
        get_single_sg_targets(genes_targets_dict, sgrna, genes_exons_seq_dict, intron_region_added, failed_targets)
    updated_genes_targets_dict = create_genes_multiplex_targets_dict_permutations(genes_targets_dict, rank, max_amplicon_len,
                                                                                  primer_length, target_surrounding_region)

    return updated_genes_targets_dict


def get_sgrna_key(sgrna, family_targeting: int):
    sgrna_key = f"{sgrna[0].seq}_{sgrna[1].seq}" if family_targeting == 2 else f"{sgrna.seq}"
    return sgrna_key


def get_sgrnas_from_crispys_output(res_grnas: SubgroupRes,
                                   genes_exons_seq_dict: Dict[
                                       str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]],
                                   family_targeting: int):
    sgrnas_dict = {}
    sgrnas_lst = res_grnas.candidates_list
    sgrnas = combinations(sgrnas_lst, 2) if family_targeting == 2 else sgrnas_lst
    sgrna_pairs_list = []
    genes_alleles_dict = {gene: [scaffold.upper() for scaffold in genes_exons_seq_dict[gene][1]] for gene in
                          genes_exons_seq_dict}
    for sgrna in sgrnas:
        sgrna_key = get_sgrna_key(sgrna, family_targeting)
        targ_alleles_num, targ_genes_num, multiplex_score, cut_genes_alleles_dict = calc_family_multiplex_score(sgrna, genes_alleles_dict, family_targeting)
        sgrnas_dict[sgrna_key] = targ_alleles_num, targ_genes_num, multiplex_score, sgrna, cut_genes_alleles_dict
        sgrna_pairs_list.append(sgrna)
    sorted_sgrnas_list = sorted(sgrna_pairs_list,
                                key=lambda sgrna: (-sgrnas_dict[get_sgrna_key(sgrna, family_targeting)][2],
                                                   -sgrnas_dict[get_sgrna_key(sgrna, family_targeting)][1],
                                                   -sgrnas_dict[get_sgrna_key(sgrna, family_targeting)][0]))
    return sgrnas_dict, sorted_sgrnas_list
