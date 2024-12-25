from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.Primers_Obj import Primers_Obj
from Amplicon_construction.SNP_Obj import SNP_Obj
from Amplicon_construction.Target_Obj import Target_Obj


primers = Primers_Obj(0.7, "ACTGACTGACTGACTGAC", "ACTGACTGACTGACTGAC", 0, 18, 260, 18)
snps = [SNP_Obj(25, {1}), SNP_Obj(250, {2})]
target = Target_Obj("ACTGACTGACTGACTGACTGCGG", 100, 122, "+")

amplicon = Amplicon_Obj("", 0, 278, 1.0, 1.0, target, snps, primers)
amplicon_dict = amplicon.to_dict()
print(amplicon_dict)
