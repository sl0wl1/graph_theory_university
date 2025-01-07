from src.clustering import group_after_invariant
from synutility.SynIO.data_type import load_from_pickle
import pytest

def test_fail_invariant():
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]
    with pytest.raises(ValueError):
        group_after_invariant(list_reactions=data, invariant="I am not an invariant")
        
def test_sum_invariants_vertex_degrees():
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]
    result = group_after_invariant(data, invariant="vertex_degrees")
    result_flattened = [entry for entries in result.values() for entry in entries]
    assert len(result_flattened) == len(data)
