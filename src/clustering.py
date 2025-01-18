from typing import Dict, List, Any
import networkx as nx
import networkx.algorithms.isomorphism as iso


from src.rc_extract import get_rc_updated
from src.invariants import (
    edge_count_invariant,
    vertex_degree_invariant,
    vertex_count_invariant,
    algebraic_connectivity_invariant,
    rank_invariant,
    histogram_invariant_check,
)
from src.add_combined_node_attributes import combine_charge_element_to_node
from src.weisfeiler_lehman_si import weisfeiler_lehman_isomorhpic_test, SharedHashTable


def cluster_reactions(list_reactions: List[Dict[Any, Any]]) -> Dict[str, Any]:
    """Simple function for clusterting chemical reactions

    Args:
        list_reactions (List[Dict[Any, Any]]): A list of reactions

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the cluster. Values are the isomorphic reactions.
    """

    # Create callables for is_isomorphic check
    nm_charge = iso.numerical_node_match("charge", 0)
    nm_element = iso.categorical_node_match("element", "C")
    em_order = iso.categorical_edge_match("order", 0)

    # Create an empty dict for storing reaction clusters and cluster counter for key naming
    cluster_dict = {}
    cluster_counter = 0

    # Check if the list of reactions is empty
    if list_reactions:
        for idx, reaction in enumerate(list_reactions):
            reaction_centre = get_rc_updated(reaction["ITS"])

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1

            else:
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    cluster_centre = get_rc_updated(value[0]["ITS"])

                    # TODO Maybe this could be optimized... and I do not know if it is correct -.-'
                    # Only need to be checked once
                    if (
                        nx.is_isomorphic(
                            cluster_centre, reaction_centre, node_match=nm_charge
                        )
                        and nx.is_isomorphic(
                            cluster_centre, reaction_centre, node_match=nm_element
                        )
                        and nx.is_isomorphic(
                            cluster_centre, reaction_centre, edge_match=em_order
                        )
                    ):
                        value.append(reaction)
                        cluster_dict[key] = value
                        break

                else:
                    # If no isomorphic reaction centre can be found, create a new entry (cluster)
                    cluster_dict[f"cluster_{cluster_counter}"] = [reaction]
                    cluster_counter += 1

    return cluster_dict


def group_after_invariant(
    list_reactions: List[Dict[Any, Any]], invariant: str
) -> Dict[str, Any]:
    """Simple function for grouping chemical reactions

    Args:
        list_reactions (List[Dict[Any, Any]]): A list of reactions
        invariant (str): Select invariant for clusterting

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the groups. Values are the isomorphic reactions.
    """

    invariants = [
        "vertex_counts",
        "edge_counts",
        "vertex_degrees",
        "algebraic_connectivity",
        "rank",
    ]

    if invariant not in invariants:
        raise ValueError("Not a valid invariant")

    # Create an empty dict for storing reaction groups and groups counter for key naming
    group_dict = {}
    group_counter = 0

    # Check if the list of reactions is empty
    if list_reactions:
        for idx, reaction in enumerate(list_reactions):
            reaction_centre = get_rc_updated(reaction["ITS"])

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                group_dict[f"group_{group_counter}"] = [list_reactions[idx]]
                group_counter += 1

            else:
                # Checks if invariants of the reaction centre already exist in a group
                for key, value in group_dict.items():
                    group_centre = get_rc_updated(value[0]["ITS"])

                    # TODO Maybe this could be optimized... and I do not know if it is correct -.-'
                    # Only need to be checked once
                    # Using switch cases for invariants
                    match invariant:
                        case "vertex_degrees":
                            group_centre_invariant, reaction_centre_invariant = (
                                vertex_degree_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )

                        # TODO Implement
                        case "vertex_counts":
                            group_centre_invariant, reaction_centre_invariant = (
                                vertex_count_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )

                        # TODO implement
                        case "edge_counts":
                            group_centre_invariant, reaction_centre_invariant = (
                                edge_count_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )

                        # TODO implement
                        case "algebraic_connectivity":
                            (
                                group_centre_invariant,
                                reaction_centre_invariant,
                            ) = (  # invariant doesn't exactly fit here (should be connectivity)
                                algebraic_connectivity_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )

                        # TODO check
                        case "rank":
                            group_centre_invariant, reaction_centre_invariant = (
                                rank_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )

                            # group_centre_invariant = np.linalg.eigvals(group_centre)
                            # reaction_centre_invariant = np.linalg.eigvals(reaction_centre)

                    if reaction_centre_invariant == group_centre_invariant:
                        value.append(reaction)
                        group_dict[key] = value
                        break

                else:
                    # If invariants are not equal, create a new entry (group)
                    group_dict[f"group_{group_counter}"] = [reaction]
                    group_counter += 1

    return group_dict


def cluster_after_invariant_grouping(
    group_dict: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """Function for clustering after grouping. The result dict of group_after_invariant() is clustered using isomorphism check.

    Args:
    group_dict Dict[str, Any]: The result dict of group_after_invariant()

    Returns:
        Dict[Dict[str, Any]]: Returns a dicts in a dict. Outer dicts keys are the groups number. In this group, the keys of the inner dicts are the cluster numbers.
    """
    cluster_after_group_dict: Dict[str, Dict[str, Any]] = {}

    if not group_dict:
        return cluster_after_group_dict

    for key, values in group_dict.items():
        temporary_cluster_dict = cluster_reactions(values)
        cluster_after_group_dict[key] = temporary_cluster_dict

    return cluster_after_group_dict


def cluster_weisfeiler_lehman_nx(
    list_reactions: List[Dict[Any, Any]],
    iterations: int = 3,
    use_edge_node_attr: bool = False,
) -> Dict[str, Any]:
    """Simple function for clusterting chemical reactions using Weisfeiler-Lehman from NetworkX

    Args:
        list_reactions (List[Dict[Any, Any]]): A list of reactions
        iterations (int): Number of neighbor aggregations to perform. Defaults to 3
        use_edge_node_attr (bool): Set to True for using edge and node attributes (order, charge, element). Defaults to False.

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the cluster. Values are the isomorphic reactions.
    """

    # Create an empty dict for storing reaction clusters and cluster counter for key naming
    cluster_dict = {}
    cluster_counter = 0

    # Check if the list of reactions is empty
    if list_reactions:
        for idx, reaction in enumerate(list_reactions):
            reaction_centre = get_rc_updated(reaction["ITS"])

            if use_edge_node_attr:
                if (
                    reaction_centre.graph.get(
                        "reaction_centre_wl_hash_node_edge_attributes"
                    )
                    is None
                ):
                    combine_charge_element_to_node(reaction_centre)
                    reaction_centre_hash = nx.weisfeiler_lehman_graph_hash(
                        reaction_centre,
                        iterations=iterations,
                        edge_attr="order",
                        node_attr="element_charge",
                    )
                    reaction_centre.graph[
                        "reaction_centre_wl_hash_node_edge_attributes"
                    ] = reaction_centre_hash

            else:
                if reaction_centre.graph.get("reaction_centre_wl_hash_vanilla") is None:
                    reaction_centre_hash = nx.weisfeiler_lehman_graph_hash(
                        reaction_centre, iterations=iterations
                    )
                    reaction_centre.graph["reaction_centre_wl_hash_vanilla"] = (
                        reaction_centre_hash
                    )

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1

            else:
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    cluster_centre = get_rc_updated(value[0]["ITS"])

                    cluster_centre_hash = (
                        cluster_centre.graph.get(
                            "reaction_centre_wl_hash_node_edge_attributes"
                        )
                        if use_edge_node_attr
                        else cluster_centre.graph.get("reaction_centre_wl_hash_vanilla")
                    )
                    reaction_centre_hash = (
                        reaction_centre.graph.get(
                            "reaction_centre_wl_hash_node_edge_attributes"
                        )
                        if use_edge_node_attr
                        else reaction_centre.graph.get(
                            "reaction_centre_wl_hash_vanilla"
                        )
                    )

                    # TODO Maybe this could be optimized... and I do not know if it is correct -.-'
                    # Only need to be checked once
                    if cluster_centre_hash == reaction_centre_hash:
                        value.append(reaction)
                        cluster_dict[key] = value
                        break

                else:
                    # If no isomorphic reaction centre can be found, create a new entry (cluster)
                    cluster_dict[f"cluster_{cluster_counter}"] = [reaction]
                    cluster_counter += 1

    return cluster_dict


def cluster_weisfeiler_lehman_si(
    list_reactions: List[Dict[Any, Any]],
    extract_reaction_centre: bool = True,
    reset: bool = True,
) -> Dict[str, Any]:
    """Simple function for clusterting chemical reactions using Weisfeiler-Lehman from NetworkX

    Args:
        list_reactions (List[Dict[Any, Any]]): A list of reactions
        iterations (int): Number of neighbor aggregations to perform. Defaults to 3
        use_edge_node_attr (bool): Set to True for using edge and node attributes (order, charge, element). Defaults to False.

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the cluster. Values are the isomorphic reactions.
    """

    # Create an empty dict for storing reaction clusters and cluster counter for key naming
    cluster_dict = {}
    cluster_counter = 0

    # Check if the list of reactions is empty
    if list_reactions:
        for idx, reaction in enumerate(list_reactions):
            reaction_centre = get_rc_updated(reaction["ITS"])

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1

            else:
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    cluster_centre = get_rc_updated(value[0]["ITS"])

                    # TODO Maybe this could be optimized... and I do not know if it is correct -.-'
                    # Only need to be checked once

                    shared_hash_table = SharedHashTable()
                    if weisfeiler_lehman_isomorhpic_test(
                        cluster_centre,
                        reaction_centre,
                        shared_hash_table,
                        extract_reaction_centre=extract_reaction_centre,
                        reset=reset,
                    ):
                        value.append(reaction)
                        cluster_dict[key] = value
                        break

                else:
                    # If no isomorphic reaction centre can be found, create a new entry (cluster)
                    cluster_dict[f"cluster_{cluster_counter}"] = [reaction]
                    cluster_counter += 1

    return cluster_dict


def cluster_histograms(
    list_reactions: List[Dict[Any, Any]],
) -> Dict[str, Any]:
    """Simple function for clusterting chemical reactions using Weisfeiler-Lehman from NetworkX

    Args:
        list_reactions (List[Dict[Any, Any]]): A list of reactions
        iterations (int): Number of neighbor aggregations to perform. Defaults to 3
        use_edge_node_attr (bool): Set to True for using edge and node attributes (order, charge, element). Defaults to False.

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the cluster. Values are the isomorphic reactions.
    """

    # Create an empty dict for storing reaction clusters and cluster counter for key naming
    cluster_dict = {}
    cluster_counter = 0

    # Check if the list of reactions is empty
    if list_reactions:
        for idx, reaction in enumerate(list_reactions):
            reaction_histogram = reaction["histogram"]

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1

            else:
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    cluster_histogram = value[0]["histogram"]

                    # TODO Maybe this could be optimized... and I do not know if it is correct -.-'
                    # Only need to be checked once

                    if histogram_invariant_check(cluster_histogram, reaction_histogram):
                        value.append(reaction)
                        cluster_dict[key] = value
                        break

                else:
                    # If no isomorphic reaction centre can be found, create a new entry (cluster)
                    cluster_dict[f"cluster_{cluster_counter}"] = [reaction]
                    cluster_counter += 1

    return cluster_dict
