import networkx as nx


# This code is based on output of Chat-GPT.

def tryConstruction(G):
    g = nx.girth(G)

    for edge in G.edges:
        e = edge
        u,v = e
        break

    # Step 1: Identify the edge e and extract the endpoints u and v
    u, v = e

    # Step 2: Calculate (g+1)/4 and (g-1)/2 as we will use them multiple times
    dist_V1 = (g + 1) // 4 - 1
    dist_V2 = (g + 1) // 4
    exclusion_dist = (g - 1) // 2

    # Step 3: Find all nodes that belong to V1 and V2
    # Use shortest_path_length to get distances from u and v
    distance_u = nx.shortest_path_length(G, source=u)
    distance_v = nx.shortest_path_length(G, source=v)
    
    V1 = {node for node in G if min(distance_u.get(node, float('inf')), distance_v.get(node, float('inf'))) <= dist_V1}
    V2 = {node for node in G if min(distance_u.get(node, float('inf')), distance_v.get(node, float('inf'))) == dist_V2}

    # Step 4: Define a function to check the lower distance avoiding nodes at exclusion_dist
    def lower_distance(w1, w2):
        # Get nodes that are at exclusion_dist from the edge endpoints (u and v)
        exclusion_nodes = {node for node in G if min(distance_u.get(node, float('inf')), distance_v.get(node, float('inf'))) == exclusion_dist}
        # Create a subgraph excluding these nodes
        subgraph = G.subgraph(set(G.nodes()) - exclusion_nodes)
        try:
            # Compute shortest path in the subgraph
            return nx.shortest_path_length(subgraph, source=w1, target=w2)
        except nx.NetworkXNoPath:
            return float('inf')  # Return infinity if there's no path

    # Step 5: Add edges between pairs in V2 if their lower distance is exactly 2
    V2copy = V2.copy()
    for w1 in V2:
        for w2 in V2:
            if w1 < w2 and w1 in V2copy and w2 in V2copy and lower_distance(w1, w2) == 2:
                G.add_edge(w1, w2)
                V2copy.remove(w1)
                V2copy.remove(w2)
                break

    assert(len(V2copy) == 0)

    # Step 4: Remove all vertices in V1 from the graph
    G.remove_nodes_from(V1)

    return G



def printInfo(G):
    print(f"g={nx.girth(G)}")
    print(f"n={nx.number_of_nodes(G)}")
    print(f"degrees {nx.degree_histogram(G)}")



if __name__ == "__main__":
    file_path = "data/5_12.s6"
    G = nx.read_sparse6(file_path)
    # printInfo(G)
    newG = tryConstruction(G)
    # printInfo(newG)
    print(nx.to_sparse6_bytes(newG, header=False).decode("utf-8"))