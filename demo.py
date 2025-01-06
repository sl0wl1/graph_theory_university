from src.clustering import cluster_reactions
from synutility.SynIO.data_type import load_from_pickle

# Load first 1000 entries
data = load_from_pickle("data/ITS_graphs.pkl.gz")[:1000]

cluster_dict = cluster_reactions(data)

print(cluster_dict.keys())