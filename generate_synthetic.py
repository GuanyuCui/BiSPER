import networkx as nx

def save_edges(edge_list, filename):
    with open(filename, 'w') as f:
        for edge in edge_list:
            f.write(f"{edge[0]} {edge[1]}\n")

if __name__ == "__main__":
    n = 5000
    p = 0.005
    er_graph = nx.erdos_renyi_graph(n, p)

    # Find largest connected component induced graph.
    largest_component_nodes = max(nx.connected_components(er_graph), key = len)
    largest_component_subgraph = er_graph.subgraph(largest_component_nodes)
    largest_component_edges = list(largest_component_subgraph.edges())

    # Save.
    save_edges(largest_component_edges, "synthetic.txt")
