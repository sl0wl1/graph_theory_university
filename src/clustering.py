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
)
from src.add_combined_node_attributes import combine_charge_element_to_node
from src.types_used import Config
from src.weisfeiler_lehman_si import weisfeiler_lehman_isomorphic_test


def cluster_by_isomorphism_test(list_reactions: List[Dict[Any, Any]]) -> Dict[str, Any]:
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
            reaction_centre = reaction[
                "reaction_centre"
            ]  # change to expect get_rc_updated to have already run

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1

            else:
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    cluster_centre = get_rc_updated(value[0]["ITS"])
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

# TODO: not yet implemented
def cluster_by_weisfeiler_lehman_si(
    list_reactions: List[Dict[Any, Any]],
) -> Dict[str, Any]:
   pass


def cluster_by_weisfeiler_lehman_nx(
    list_reactions: List[Dict[Any, Any]],
    **kwargs,
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
    weisfeiler_lehman_params = kwargs.get(
        "weisfeiler_lehman_params", {"iterations": 3, "use_edge_node_attr": False}
    )
    iterations = weisfeiler_lehman_params.get("iterations")
    use_edge_node_attr = weisfeiler_lehman_params.get("use_edge_node_attr")

    cluster_dict = {}
    cluster_counter = 0

    # Check if the list of reactions is empty
    if list_reactions:
        for idx, reaction in enumerate(list_reactions):
            reaction_centre = reaction["reaction_centre"]

            if use_edge_node_attr:
                combine_charge_element_to_node(reaction_centre)
                reaction_centre_hash = nx.weisfeiler_lehman_graph_hash(
                    reaction_centre,
                    iterations=iterations,
                    edge_attr="order",
                    node_attr="element_charge",
                )
                reaction_centre.graph["reaction_centre_wl_hash"] = reaction_centre_hash

            else:
                reaction_centre_hash = nx.weisfeiler_lehman_graph_hash(
                    reaction_centre, iterations=iterations
                )
                reaction_centre.graph["reaction_centre_wl_hash"] = reaction_centre_hash

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1

            else:
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    cluster_centre = value[0]["reaction_centre"]
                    if (
                        cluster_centre.graph["reaction_centre_wl_hash"]
                        == reaction_centre.graph["reaction_centre_wl_hash"]
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
        "edge_count",
        "vertex_degree",
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
            reaction_centre = reaction["reaction_centre"]

            # Create first entry in dict. For the first reaction there is nothing to compare
            if idx == 0:
                group_dict[f"group_{group_counter}"] = [list_reactions[idx]]
                group_counter += 1

            else:
                # Checks if invariants of the reaction centre already exist in a group
                for key, value in group_dict.items():
                    group_centre = get_rc_updated(value[0]["ITS"])

                    match invariant:
                        case "vertex_degree":
                            group_centre_invariant, reaction_centre_invariant = (
                                vertex_degree_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )
                        case "vertex_count":
                            group_centre_invariant, reaction_centre_invariant = (
                                vertex_count_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )
                        case "edge_count":
                            group_centre_invariant, reaction_centre_invariant = (
                                edge_count_invariant(
                                    group_centre=group_centre,
                                    reaction_centre=reaction_centre,
                                )
                            )
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
    config: Config,
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
        match config.algorithm:
            case "none":
                temporary_cluster_dict = values
            case "isomorphism_test":
                cluster_by_isomorphism_test(list_reactions=values)
            case "weisfeiler_lehmann_nx":
                temporary_cluster_dict = cluster_by_weisfeiler_lehman_nx(
                    list_reactions=values, **config.weisfeiler_lehman_params
                )
            case "weisfeiler_lehmann_si":
                temporary_cluster_dict = cluster_by_weisfeiler_lehman_si(
                    list_reactions=values
                )
        cluster_after_group_dict[key] = temporary_cluster_dict

    return cluster_after_group_dict

def cluster_without_invariant_grouping(
    list_reactions: List[Dict[Any, Any]],
    config: Config,
) -> Dict[str, Any]:
    """Function for clustering without grouping. The result dict of group_after_invariant() is clustered using isomorphism check.

    Args:
    list_reactions List[Dict[Any, Any]]: A list of reactions

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the cluster. Values are the isomorphic reactions.
    """

    match config.algorithm:
        case "none":
            cluster_dict = list_reactions
        case "isomorphism_test":
            cluster_dict = cluster_by_isomorphism_test(list_reactions=list_reactions)
        case "weisfeiler_lehmann_nx":
            cluster_dict = cluster_by_weisfeiler_lehman_nx(
                list_reactions=list_reactions, **config.weisfeiler_lehman_params
            )
        case "weisfeiler_lehmann_si":
            cluster_dict = cluster_by_weisfeiler_lehman_si(list_reactions=list_reactions)

    return cluster_dict