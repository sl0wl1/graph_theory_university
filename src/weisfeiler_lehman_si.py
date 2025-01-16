from typing import Any, Callable, Dict, Tuple
import networkx as nx
import numpy as np
from collections import Counter
from src.rc_extract import get_rc_updated


class SharedHashTable:
    def __init__(self, hash_function: Callable | None = None) -> None:
        self.hash_function = hash_function
        self._increment_hash = 1
        self.hash_function_exists = hash_function is not None
        self.shared_hash_table: Dict[Any, Any] = dict()

    def set(self, key) -> None:
        if key not in self.shared_hash_table:
            if self.hash_function_exists:
                self.shared_hash_table[key] = self.hash_function(key)
            else:
                self.shared_hash_table[key] = self._increment_hash
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
        return tuple(
            (graph.nodes[node]["compressed_label"] for node in graph)
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


def weisfeiler_lehman_isomorphic_test(
    graph_1: nx.Graph,
    graph_2: nx.Graph,
    shared_hash_table: SharedHashTable = SharedHashTable(),
    extract_reaction_center: bool = False,
    reset: bool = False,
) -> bool:
    if extract_reaction_center:
        graph_1 = get_rc_updated(graph_1)
        graph_2 = get_rc_updated(graph_2)

    temporary_histogram_1: Dict[int, int] = {}
    temporary_histogram_2: Dict[int, int] = {}

    for _ in range(len(graph_1.nodes)):
        temporary_compressed_labels_1, temporary_histogram_1 = weisfeiler_lehman_step(
            graph=graph_1, shared_hash_table=shared_hash_table, reset=reset
        )
        temporary_compressed_labels_2, temporary_histogram_2 = weisfeiler_lehman_step(
            graph=graph_2, shared_hash_table=shared_hash_table, reset=reset
        )

        if temporary_compressed_labels_1 != temporary_compressed_labels_2:
            return False

        if temporary_histogram_1 == temporary_histogram_2:
            return True

    return True
