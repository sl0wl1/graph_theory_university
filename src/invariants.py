from typing import List, Tuple
import networkx as nx


def vertex_degree_invariant(
    group_centre: nx.Graph, reaction_centre: nx.Graph
) -> Tuple[List[int], List[int]]:
    """Implementation for comparing vertex degree of two graphs, used for invariant comparison. Used in group_after_invariant

    Args:

        group_centre (nx.Graph): group centre from your already grouped reactions list
        reaction_centre (nx.Graph): reaction centre from your reactions list

    Returns:
        group_centre_invariant, reaction_centre_invariant Tuple[List[int], List[int]]: vertex degree list of inputs
    """

    group_centre_invariant = list(dict(group_centre.degree).values())
    group_centre_invariant.sort()

    reaction_centre_invariant = list(dict(reaction_centre.degree).values())
    reaction_centre_invariant.sort()

    return group_centre_invariant, reaction_centre_invariant
