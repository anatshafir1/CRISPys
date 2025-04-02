from FindTargets import create_sgrna_permutations, calculate_off_scores, create_multiplex_targets
from Target_Obj import Combined_Target_Obj, Target_Obj


trg1 = "tgaaccacaacaaaattcattgg".upper()
trg2 = "tgaactacaacagaattcattgg".upper()
trg3 = "tgaaccacaacagaattcattgg".upper()

trg1_ungapped = trg1
trg2_ungapped = trg2
trg3_ungapped = trg3

trg_obj1 = Target_Obj(trg1, 10, 33, "+", "scf_1", trg1_ungapped)
trg_obj2 = Target_Obj(trg2, 10, 33, "+", "scf_2", trg2_ungapped)
trg_obj3 = Target_Obj(trg3, 10, 33, "+", "scf_3", trg3_ungapped)

comb_trg_obj1 = Combined_Target_Obj(10, 33, [trg_obj1, trg_obj2, trg_obj3])

trg4 = "tgaaccacaacagaattcattgg".upper()
trg5 = "tgaaccacaacagaattcattgg".upper()
trg6 = "tgaaccacaacagaattcattgg".upper()

trg4_ungapped = trg4
trg5_ungapped = trg5
trg6_ungapped = trg6

trg_obj1 = Target_Obj(trg1, 60, 83, "+", "scf_1", trg4_ungapped)
trg_obj2 = Target_Obj(trg2, 60, 83, "+", "scf_2", trg5_ungapped)
trg_obj3 = Target_Obj(trg3, 60, 83, "+", "scf_3", trg6_ungapped)

comb_trg_obj2 = Combined_Target_Obj(60, 83, [trg_obj1, trg_obj2, trg_obj3])

relevant_targets_dict = {1: [comb_trg_obj1, comb_trg_obj2]}
allele_ids_lst = ["scf_1", "scf_2", "scf_3"]

create_sgrna_permutations(relevant_targets_dict)
calculate_off_scores(relevant_targets_dict)
new_relevant_targets_dict = create_multiplex_targets(relevant_targets_dict, allele_ids_lst)

print(new_relevant_targets_dict)
