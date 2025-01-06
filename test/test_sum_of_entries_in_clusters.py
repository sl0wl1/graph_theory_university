from src.clustering import cluster_reactions
from synutility.SynIO.data_type import load_from_pickle

def test_clustering():
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]
    result = cluster_reactions(data)
    result_flattened = [entries for entry in result.values() for entries in entry]
    assert len(result_flattened) == len(data)
