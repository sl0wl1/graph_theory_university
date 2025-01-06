import networkx as nx
from typing import List, Set


def find_unequal_order_edges(G: nx.Graph) -> List[int]:
    """
    Identifies nodes that are part of edges with unequal orders in a graph, which are
    indicative of reaction centers.

    Parameters:
    - G (nx.Graph): The graph in which to identify nodes associated with reaction centers.

    Returns:
    - List[int]: List of node indices that are part of edges with unequal order, indicating
    reaction centers.
    """
    reaction_center_nodes: Set[int] = set()
    for u, v, data in G.edges(data=True):
        order = data.get("order", (1, 1))
        if isinstance(order, tuple) and order[0] != order[1]:
            if data.get("standard_order") != 0:
                reaction_center_nodes.add(u)
                reaction_center_nodes.add(v)
    return list(reaction_center_nodes)


def extract_subgraph(G: nx.Graph, node_indices: List[int]) -> nx.Graph:
    """
    Extracts a subgraph based on specified node indices from the given graph.

    Parameters:
    - G (nx.Graph): The original graph.
    - node_indices (List[int]): Indices of nodes to include in the subgraph.

    Returns:
    - nx.Graph: The extracted subgraph containing only the specified nodes.
    """
    return G.subgraph(node_indices).copy()


def get_rc(G: nx.Graph) -> nx.Graph:
    """
    Generates a subgraph of reaction centers from the input graph based on unequal edge orders.

    Parameters:
    - G (nx.Graph): The graph from which reaction centers are to be identified.

    Returns:
    - nx.Graph: A subgraph containing only the nodes involved in reaction centers.
    """
    rc_indices = find_unequal_order_edges(G)
    return extract_subgraph(G, rc_indices)
