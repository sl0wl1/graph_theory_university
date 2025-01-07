import networkx as nx


def combine_charge_element_to_node(graph: nx.Graph) -> None:
    """Simple function for combining element and charge attribute for Weisfeiler-Lehman hashing. Be aware that this is in-place. Data will be changed.

    Args:
        graph (nx.Graph): The graph for adding

    Returns:
        None
    """

    for node in graph.nodes:
        graph.nodes[node]["element_charge"] = (
            f"{graph.nodes[node]["element"]}, {graph.nodes[node]["charge"]}"
        )
