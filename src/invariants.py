from typing import List, Tuple
import networkx as nx
from networkx import algebraic_connectivity
import numpy as np


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


def algebraic_connectivity_invariant(
    group_centre: nx.Graph, reaction_centre: nx.Graph
) -> Tuple[float, float]:
    def compute_algebraic_connectivity(graph: nx.Graph):
        try:
            return algebraic_connectivity(graph, normalized=True, tol=1e-6)
        except nx.NetworkXError:
            return 0  # disconnected
        except nx.NetworkXNotImplemented:  # when G is directed
            return 0  # probably better to handle this differently

    group_centre_connectivity = compute_algebraic_connectivity(group_centre)
    reaction_centre_connectivity = compute_algebraic_connectivity(reaction_centre)

    return group_centre_connectivity, reaction_centre_connectivity


def vertex_count_invariant(
    group_centre: nx.Graph, reaction_centre: nx.Graph
) -> Tuple[List[int], List[int]]:
    """Implementation for comparing vertex counts of two graphs, used for invariant comparison. Used in group_after_invariant

    Args:

        group_centre (nx.Graph): group centre from your already grouped reactions list
        reaction_centre (nx.Graph): reaction centre from your reactions list

    Returns:
        group_centre_invariant, reaction_centre_invariant Tuple[List[int], List[int]]: vertex count list of inputs
    """
    group_centre_invariant = [group_centre.number_of_nodes()]

    reaction_centre_invariant = [reaction_centre.number_of_nodes()]

    return group_centre_invariant, reaction_centre_invariant


def edge_count_invariant(
    group_centre: nx.Graph, reaction_centre: nx.Graph
) -> Tuple[List[int], List[int]]:
    """Implementation for comparing edge counts of two graphs, used for invariant comparison. Used in group_after_invariant

    Args:

        group_centre (nx.Graph): group centre from your already grouped reactions list
        reaction_centre (nx.Graph): reaction centre from your reactions list

    Returns:
        group_centre_invariant, reaction_centre_invariant Tuple[List[int], List[int]]: edge count list of inputs
    """
    group_centre_invariant = [group_centre.number_of_edges()]

    reaction_centre_invariant = [reaction_centre.number_of_edges()]

    return group_centre_invariant, reaction_centre_invariant


def rank_invariant(
    group_centre: nx.Graph, reaction_centre: nx.Graph
) -> Tuple[List[int], List[int]]:
    """Implementation for comparing ranks of two graphs (adjacency matrix), used for invariant comparison. Used in group_after_invariant

    Args:

        group_centre (nx.Graph): group centre from your already grouped reactions list
        reaction_centre (nx.Graph): reaction centre from your reactions list

    Returns:
        group_centre_invariant, reaction_centre_invariant Tuple[List[int], List[int]]: edge count list of inputs
    """
    group_centre_invariant = int(np.linalg.matrix_rank(nx.to_numpy_array(group_centre)))

    reaction_centre_invariant = int(
        np.linalg.matrix_rank(nx.to_numpy_array(reaction_centre))
    )

    return group_centre_invariant, reaction_centre_invariant


def histogram_invariant_check(histogram_1, histogram_2):
    list_1 = list(histogram_1.values())
    list_1.sort()
    list_2 = list(histogram_2.values())
    list_2.sort()
    if list_1 == list_2:
        return True
    else:
        return False
