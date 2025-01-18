from typing import Any, Callable, Dict, Tuple
import networkx as nx
from src.rc_extract import get_rc_updated
from collections import Counter


class SharedHashTable:
    def __init__(self, hash_function: Callable | None = None) -> None:
        if hash_function:
            self.hash_function = hash_function
        else:
            self._increment_hash = 1
        self.hash_function_exists = True if hash_function else False
        self.shared_hash_table: Dict[Any, Any] = dict()

    def set(self, key) -> None:
        if key not in self.shared_hash_table.keys():
            self.shared_hash_table[key] = (
                hash(key) if self.hash_function_exists else self._increment_hash
            )

            if not self.hash_function_exists:
                self._increment_hash += 1

    def get(self, key):
        return self.shared_hash_table.get(key, None)

    def __str__(self):
        return f"{self.shared_hash_table}"

    def __repr__(self):
        return f"{self.shared_hash_table}"


def weisfeiler_lehman_step(
    graph: nx.Graph, shared_hash_table: SharedHashTable, reset: bool = False
) -> Tuple[Tuple[int], Dict[int, int]]:
    # Check if you are in iteration step 0; if so, set compressed label of nodes to 1
    list_data_entries = tuple(
        (graph.nodes[node].get("compressed_label", None) for node in graph)
    )

    if None in list_data_entries or reset:
        shared_hash_table.set("initial")
        initial_hash = shared_hash_table.get("initial")
        nx.set_node_attributes(graph, initial_hash, "compressed_label")
        shared_hash_table.set("initial")

        return tuple(
            (value for value in dict(graph.nodes.data("compressed_label")).values())
        ), {}

    list_updated_compressed_labels = []

    for node in graph.nodes:
        compressed_label = graph.nodes[node]["compressed_label"]
        temporary_multiset_of_compressed_label_from_neighbours = [
            graph.nodes[neighbor].get("compressed_label", None)
            for neighbor in graph.neighbors(node)
        ]
        temporary_multiset_of_compressed_label_from_neighbours.sort()
        new_key = (
            compressed_label,
            tuple(temporary_multiset_of_compressed_label_from_neighbours),
        )
        shared_hash_table.set(new_key)
        new_hash = shared_hash_table.get(new_key)
        list_updated_compressed_labels.append(new_hash)

    # Update compressed node labels
    for idx, node in enumerate(graph.nodes):
        graph.nodes[node]["compressed_label"] = list_updated_compressed_labels[idx]

    list_updated_compressed_labels.sort()
    counter = Counter(list_updated_compressed_labels)
    histogram = dict(counter.items())

    return tuple(list_updated_compressed_labels), histogram


def weisfeiler_lehman_isomorhpic_test(
    graph_1: nx.Graph,
    graph_2: nx.Graph,
    shared_hash_table,
    extract_reaction_centre: bool = False,
    reset: bool = True,
) -> bool:
    if extract_reaction_centre:
        graph_1 = get_rc_updated(graph_1)
        graph_2 = get_rc_updated(graph_2)

    temporary_histogram_1: Dict[int, int] = {}
    temporary_histogram_2: Dict[int, int] = {}

    for i in range(len(graph_1.nodes) + 1):
        # Initial step
        if i == 0:
            temporary_compressed_labels_1, temporary_histogram_1 = (
                weisfeiler_lehman_step(
                    graph=graph_1, shared_hash_table=shared_hash_table, reset=reset
                )
            )
            temporary_compressed_labels_2, temporary_histogram_2 = (
                weisfeiler_lehman_step(
                    graph=graph_2, shared_hash_table=shared_hash_table, reset=reset
                )
            )

        if i >= 1:
            last_histogram_1_values = list(temporary_histogram_1.values())
            last_histogram_2_values = list(temporary_histogram_2.values())

        temporary_compressed_labels_1, temporary_histogram_1 = weisfeiler_lehman_step(
            graph=graph_1, shared_hash_table=shared_hash_table
        )
        temporary_compressed_labels_2, temporary_histogram_2 = weisfeiler_lehman_step(
            graph=graph_2, shared_hash_table=shared_hash_table
        )

        temporary_histogram_1_values = list(temporary_histogram_1.values())
        temporary_histogram_2_values = list(temporary_histogram_2.values())
        # Check if the ordered compressed label multiset is the same. If not -> not isomorphic
        if temporary_compressed_labels_1 != temporary_compressed_labels_2:
            return False

        # Check if partitions did not change; if so check if compressed labels are the same -> True, else -> False. This can save iterations.
        if i > 1:
            if (last_histogram_1_values == temporary_histogram_1_values) and (
                last_histogram_2_values == temporary_histogram_2_values
            ):
                if temporary_compressed_labels_1 != temporary_compressed_labels_2:
                    return False
                else:
                    return True
    else:
        return True
