import sys


def read_text_file(infile: str) -> list:
    seqs = []
    with open(infile, "r") as f:
        matrixLines = f.readlines()
        for line in matrixLines:
            seqs.append(line.strip())
    return seqs


def get_kmer_counts(seqs: list, k: int) -> dict:
    kmers = {}
    for seq in seqs:
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i : i + k]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
    return kmers


def distance(a: str, b: str) -> int:
    d = 0
    for a_i, b_i in zip(a, b):
        d = d + 1 if a_i != b_i else d
    return d


def optimal_kmer(common_kmers: list, replace: str, d: int) -> str:
    best_d = d + 1
    best_str = replace
    for common in common_kmers:
        this_d = distance(common, replace)
        if this_d > d:
            continue
        if this_d < best_d:
            best_d = this_d
            best_str = common
            continue
    return best_str


def correct(seqs: list, k: int, d: int, t: int) -> list:
    kmer_dict = get_kmer_counts(seqs, k)

    common_kmers = []
    for kmer, count in kmer_dict.items():
        if count >= t:
            common_kmers.append(kmer)

    new_seqs = []
    for num, seq in enumerate(seqs):
        corrected = set()
        # print(f"SEQUENCE {num}")
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i : i + k]
            if kmer_dict.get(kmer, 0) < t:
                # print(f"correcting {kmer}, {i}")
                if i in corrected:
                    # print(f"skipping correction of {kmer}, {i}\n")
                    continue
                correction = optimal_kmer(common_kmers, kmer, d)
                for j in range(i, i + k):
                    corrected.add(j)
                seq = seq[:i] + correction + seq[i + k :]
        new_seqs.append(seq)
    return new_seqs


if __name__ == "__main__":
    seqs = read_text_file(sys.argv[1])
    k = 15
    t = 2
    d = 2

    corrected = correct(seqs, k, d, t)
    for seq in corrected:
        print(seq)
