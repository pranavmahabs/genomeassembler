import networkx as nx
import sys


def extract_reads(file_path: str):
    reads = []
    with open(file_path, "r") as f:
        for line in f:
            reads.append(line.strip())
    return reads


def kmerize_reads(reads: list, k: int):
    kmers = set()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i : i + k]
            kmers.add(kmer)
    return list(kmers)


def build_debruijn(kmers: list, k: int):
    graph = nx.MultiDiGraph()
    nodes = set()
    # Build the graph.
    for kmer in kmers:
        n1, n2 = kmer[0 : k - 1], kmer[1:]
        nodes.add(n1)
        nodes.add(n2)
        graph.add_edge(n1, n2, label=kmer)
    nx.drawing.nx_pydot.write_dot(graph, "predebruijn.dot")
    # Simplify the Singletons
    for node in graph.nodes():
        in_e = list(graph.in_edges(node, data=True))
        ot_e = list(graph.out_edges(node, data=True))
        if len(in_e) == 1 and len(ot_e) == 1:
            in_edge = in_e[0]
            ot_edge = ot_e[0]
            in_mer = in_edge[2]["label"]
            ot_mer = ot_edge[2]["label"]
            graph.remove_edge(in_edge[0], in_edge[1])
            graph.remove_edge(ot_edge[0], ot_edge[1])
            graph.add_edge(in_edge[0], ot_edge[1], label=in_mer + ot_mer[k - 1 :])
    # Remove the Singleton Nodes from the Graph
    singeltons = [node for node in graph.nodes() if graph.degree(node) == 0]
    graph.remove_nodes_from(singeltons)
    nx.drawing.nx_pydot.write_dot(graph, "debruijn.dot")

    return graph


def dfs(graph: nx.MultiDiGraph, node: str, path: list):
    new_path = path + [node]
    visited = set(new_path)

    out_neighbors = list(graph.out_edges(node))
    out_nodes = [edge[1] for edge in out_neighbors]
    extend_dfs = []
    for neighbor in out_nodes:
        if neighbor not in visited:
            extend_dfs.append(neighbor)

    # End the recursive case.
    if len(extend_dfs) == 0:
        return [[node]]

    # Find maximal path from this position
    results = []
    for neighbor in extend_dfs:
        this_dfs = dfs(graph, neighbor, new_path)
        for pathway in this_dfs:
            results.append([node] + pathway)
    return results


def path_searching(graph: nx.MultiDiGraph, k: int):
    """
    Assumes that the created deBruijn graph is acyclic.
    """
    # Identify the starting node: OutDegree/Degree = 1
    starting_nodes = []
    for node in graph.nodes():
        if graph.out_degree(node) == 1 and graph.in_degree(node) == 0:
            starting_nodes.append(node)

    # Perform DFS on each starting node.
    paths = []
    for starting in starting_nodes:
        paths.extend(dfs(graph, starting, []))

    # Build the Potential Consensus Sequences
    sequences = []
    for path in paths:
        seq = ""
        for _, ld in graph.get_edge_data(path[0], path[1]).items():
            seq = ld["label"] if len(ld["label"]) > len(seq) else seq
        for i in range(1, len(path) - 1):
            edge = ""
            for _, ld in graph.get_edge_data(path[i], path[i + 1]).items():
                edge = ld["label"] if len(ld["label"]) > len(edge) else edge
            seq = seq + edge[k - 1 :]
        sequences.append(seq)

    # Return All Possible Sequences from the Graph
    return sequences


def inference(reads: list, k: int):
    kmers = kmerize_reads(reads, k)
    graph = build_debruijn(kmers, k)
    sequences = path_searching(graph, k)

    max_length = 0
    for seq in sequences:
        if len(seq) > max_length:
            max_length = len(seq)
    longest = []
    for seq in sequences:
        if len(seq) == max_length:
            longest.append(seq)
    return sorted(list(set(longest)))


if __name__ == "__main__":
    reads = extract_reads(sys.argv[1])
    k = len(reads[0])
    targets = inference(reads, k)
    for target in targets:
        print(target)
