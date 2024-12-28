from Amplicon_construction.Target_Obj import Combined_Target_Obj, Target_Obj, sgRNA
from Amplicon_construction.FindTargets import calc_multiplex_score

trg1 = "tgaaccacaacaaaattcattgg".upper()
trg2 = "tgaactacaacagaattcattgg".upper()
trg3 = "tgaaccacaacagaattcattgg".upper()

trg1_ungapped = trg1
trg2_ungapped = trg2
trg3_ungapped = trg3

trg_obj1 = Target_Obj(trg1, 10, 33, "+", "scf_1", trg1_ungapped)
trg_obj2 = Target_Obj(trg2, 10, 33, "+", "scf_2", trg2_ungapped)
trg_obj3 = Target_Obj(trg3, 10, 33, "+", "scf_3", trg3_ungapped)

comb_trg_obj = Combined_Target_Obj(10, 23, [trg_obj1, trg_obj2, trg_obj3])

relevant_targets_dict = {1: [comb_trg_obj]}

allele1 = "scf_1"
allele2 = "scf_2"
allele3 = "scf_3"

sg1 = sgRNA(10, 32, trg1, {allele1: 1.0, allele2: 0.4, allele3: 1.0}, [trg_obj1, trg_obj2, trg_obj3])
sg2 = sgRNA(60, 82, trg2, {allele1: 0.3, allele2: 1.0, allele3: 0.3}, [trg_obj1, trg_obj2, trg_obj3])

multiplex_score = calc_multiplex_score(sg1, sg2, [allele1, allele2, allele3])
print(multiplex_score)
