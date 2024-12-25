from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict


def genome_transformation(input_genome: str, num_comb: int, output_genome: str):

    # Step 1: Read the original FASTA file
    sequences = list(SeqIO.parse(input_genome, "fasta"))

    # Step 2: Sort sequences by length (longest to shortest)
    sequences.sort(key=lambda x: len(x.seq), reverse=True)

    # Step 3: Distribute sequences among the new scaffolds
    scaffold_bins = defaultdict(list)

    for i, seq in enumerate(sequences):
        scaffold_bins[i % num_comb].append(seq)

    # Step 4: Combine sequences, create the mapping, and write to a new FASTA file
    new_scaffolds = []
    mapping = {}

    for i in range(num_comb):
        combined_seq = []
        current_length = 0

        for seq in scaffold_bins[i]:
            combined_seq.append(str(seq.seq))
            start_index = current_length + 1
            end_index = current_length + len(seq.seq)
            mapping[seq.id] = {
                "new_scaffold": f"scaffold_{i + 1}",
                "start": start_index,
                "end": end_index
            }
            current_length = end_index

        new_record = SeqRecord(Seq("".join(combined_seq)), id=f"scaffold_{i + 1}", description="")
        new_scaffolds.append(new_record)

    SeqIO.write(new_scaffolds, output_genome, "fasta")

    print(f"New genome FASTA file with {num_comb} scaffolds created: {output_genome}")
    print(mapping)


genome_transformation("/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Banana_GAL.Phased_Scaffolds.fasta", 30, "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Banana_GAL.Phased_Scaffolds_comb.fasta")
