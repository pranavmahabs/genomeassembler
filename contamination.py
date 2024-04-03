import sys


def read_text_file(infile: str) -> list:
    seqs = []
    with open(infile, "r") as f:
        matrixLines = f.readlines()
        for line in matrixLines:
            seqs.append(line.strip())
    return seqs


def get_seeds(vector: str, k: int) -> list:
    kmers = {}
    for i in range(0, len(vector) - k + 1):
        kmer = vector[i : i + k]
        if kmer in kmers:
            kmers[kmer].append(i)
        else:
            kmers[kmer] = [i]
    return kmers


def extend(r: str, v: str, index: int, fwd: bool, k: int):
    extension_length = 0
    r_i = 0 if fwd else len(r) - 1
    index = index if fwd else index + k - 1
    while True:
        if fwd and (r_i == len(r) or index == len(v)):
            break
        if not fwd and (r_i < 0 or index < 0):
            break
        if r[r_i] == v[index]:
            extension_length += 1
            if fwd:
                r_i += 1
                index += 1
            else:
                r_i -= 1
                index -= 1
        else:
            break

    return extension_length


def decontaminate(seeds: dict, reads: list, k: int, v: str) -> list:
    idxs = []
    new_reads = []
    for i, read in enumerate(reads):
        original_length = len(read)
        start_kmer = read[0:k]
        end_kmer = read[-k:]

        max_ef = 0
        if start_kmer in seeds:
            indices = seeds[start_kmer]
            for start in indices:
                e_len = extend(read, v, start, True, k)
                max_ef = max_ef if e_len < max_ef else e_len
        read = read[max_ef:]

        max_ef = 0
        if end_kmer in seeds:
            indices = seeds[end_kmer]
            for start in indices:
                e_len = extend(read, v, start, False, k)
                max_ef = max_ef if e_len < max_ef else e_len
        if max_ef > 0:
            read = read[:-max_ef]

        if len(read) < original_length:
            idxs.append(i)
        new_reads.append(read)

    return idxs, new_reads


if __name__ == "__main__":
    reads = read_text_file(sys.argv[1])
    vector = read_text_file(sys.argv[2])[0]
    k = int(sys.argv[3])
    seeds = get_seeds(vector, k)
    idxs, new_reads = decontaminate(seeds, reads, k, vector)

    p = ""
    for i in idxs:
        p += str(i) + ","
    print(p[:-1])
    print("--------------------")
    for read in new_reads:
        print(read)
