import networkx as nx


def find_l_neighborhood_of_rc(
    graph: nx.Graph,
    reaction_centre: nx.Graph,
    l_neighborhood: int = 0,
) -> nx.Graph:
    """Recursive function to find the subgraph in your l neighborhood of your reaction centre.

    Args:
        graph (nx.Graph): This is your ITS graph.
        reaction_centre (nx.Graph): Reaction centre (subgraph) of your ITS graph (graph).
        l_neighborhood (int) = Defines the size of your neighborhood. (Defaults to 0)
    Returns:
        nx.Graph: A subgraph of your graph with the size of l_neighborhood starting from your reaction_centre.
    """

    edge_set = [edge for edge in reaction_centre.edges]

    # Base case
    if l_neighborhood <= 0:
        return reaction_centre

    # recursive part
    for node in reaction_centre.nodes:
        for n_node in graph.neighbors(node):
            edge_set.append((node, n_node))
            new_reaction_centre = nx.edge_subgraph(G=graph, edges=edge_set)
    return find_l_neighborhood_of_rc(
        graph=graph,
        reaction_centre=new_reaction_centre,
        l_neighborhood=l_neighborhood - 1,
    )
