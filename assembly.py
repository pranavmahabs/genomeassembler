from contamination import get_seeds, decontaminate, read_text_file
from correction import correct
from debruijn import inference, extract_reads

import sys

if __name__ == "__main__":
    reads = extract_reads(sys.argv[1])
    vector = read_text_file(sys.argv[2])[0]

    # Decontaminate Reads
    contam_k = 10
    seeds = get_seeds(vector, contam_k)
    idxs, new_reads = decontaminate(seeds, reads, contam_k, vector)

    # Remove Any Reads that were Found to be Contaminated.
    num_dropped = len(idxs)
    # print(f"Dropping {100 * num_dropped / len(reads)}% of the given reads.")
    drop = set(idxs)
    reads = [reads[i] for i in range(len(new_reads)) if i not in drop]

    # Perform Correction of Reads
    corr_k, corr_t, corr_d = 15, 2, 2
    corrected = correct(reads, corr_k, corr_d, corr_t)

    # Assemble the Sequences and Print Optimal Alignments
    asm_k = 38
    targets = inference(corrected, asm_k)
    assembly = print(targets[0].strip())
