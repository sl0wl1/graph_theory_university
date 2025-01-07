from src.clustering import cluster_weisfeiler_lehman_nx
from synutility.SynIO.data_type import load_from_pickle


def test_clustering_wl_nx():
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]
    result = cluster_weisfeiler_lehman_nx(data)
    result_flattened = [entry for entries in result.values() for entry in entries]
    assert len(result_flattened) == len(data)


def test_clustering_wl_nx_with_attr():
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]
    result = cluster_weisfeiler_lehman_nx(data, use_edge_node_attr=True)
    result_flattened = [entry for entries in result.values() for entry in entries]
    assert len(result_flattened) == len(data)
