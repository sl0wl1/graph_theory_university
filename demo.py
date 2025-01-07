from src.clustering import group_after_invariant
from synutility.SynIO.data_type import load_from_pickle

# Load first 1000 entries
data = load_from_pickle("data/ITS_graphs.pkl.gz")  # [:1000]

# cluster_dict = cluster_reactions(data)

invariant_dict = group_after_invariant(data, invariant="vertex_degrees")
# cluster_after_invariant_grouping_dict = cluster_after_invariant_grouping(invariant_dict)

for key, values in invariant_dict.items():
    print(len(values))
