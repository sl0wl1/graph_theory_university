from src.clustering import cluster_weisfeiler_lehman_nx
from synutility.SynIO.data_type import load_from_pickle

# Load first 1000 entries
data = load_from_pickle("data/ITS_graphs.pkl.gz")  # [:1000]

# cluster_dict = cluster_reactions(data)

invariant_dict = cluster_weisfeiler_lehman_nx(data, use_edge_node_attr=True)
# cluster_after_invariant_grouping_dict = cluster_after_invariant_grouping(invariant_dict)

for key, values in invariant_dict.items():
    print(len(values))
