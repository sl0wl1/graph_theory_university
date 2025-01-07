from src.clustering import cluster_reactions, group_after_invariant, cluster_after_invariant_grouping
from synutility.SynIO.data_type import load_from_pickle

def test_clustering():
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]
    result = group_after_invariant(data, invariant="vertex_degrees")
    result = cluster_after_invariant_grouping(result)
    result_flattened = [entry_entry for entries in result.values() for entry in entries.values() for entry_entry in entry]
    
    assert len(result_flattened) == len(data)