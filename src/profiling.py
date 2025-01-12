import pickle


# Type Definitions:
# - ReactionDict: Dict[str, Any]
# 	- "ITS": nx.Graph of Imaginary Transition State
# 	- "R-id": Any
# 	- "class": Optional[Any]
# - DataList: List[ReactionDict]
# This script performs the following tasks:
# 1. Loads reaction data from a pickle file.
# 2. Processes the reaction data to add reaction centres.
# 3. Clusters the reactions based on their properties.
# 4. Saves the clustered reactions to a new pickle file.

# Functions:
# - benchmark(data: DataList) -> List[dict]:
# 	Benchmark a selected clustering algorithm.
from typing import Dict, List, Any, Union

import networkx as nx
from src.rc_extract import get_rc_updated
from src.clustering import cluster_reactions, cluster_weisfeiler_lehman_nx
import copy
from enum import Enum

with open("ITS_graphs.pkl.gz", "rb") as file:
    data = pickle.load(file)


# generate a list of reaction graphs
# find a way to recover the reaction id from the reaction centre
# consider adding key to data dict with rc object
# pickle for quick


# retrieve reaction centre
# in next iteration consider making it independent of reaction centre method.
# i.e. to include WP6 neighbourhood function

# cluster method should receive the list of dicts and return a copy of this
# good so that we can still relate the id of the reaction to the reaction centre and the cluster
# to which it belongs
# Generate a list of reaction graphs

ReactionDict = Dict[str, Any]
ReactionDataList = List[ReactionDict]
ClusterDict = Dict[str, List[ReactionDict]]


class ComparisonMethod(Enum):
    VD = "vertex_degree"
    VC = "vertex_count"
    AC = "algebraic_connectivity"
    EC = "edge_count"
    RA = "rank"
    IS = "isomorphism"
    WL = "weisfeiler_lehman"


# run cprofile on this method and maybe even find a way to add the log to the data pickle


def run_clustering(
    raw_data: ReactionDataList, comparison_method: ComparisonMethod
) -> ClusterDict:
    """
    Benchmark a selected clustering algorithm
    """
    reactions = copy.deepcopy(raw_data)
    for reaction in reactions:
        reaction.pop("class", None)  # don't need the class property

    for reaction in reactions:
        reaction_centre = get_rc_updated(reaction["ITS"])
        reaction["reaction_centre"] = reaction_centre
        # save this as a pickle file for reuseability or visualisation

    match comparison_method:
        case (
            ComparisonMethod.VD
            | ComparisonMethod.VC
            | ComparisonMethod.AC
            | ComparisonMethod.EC
            | ComparisonMethod.RA
        ):
            # TODO: refactor invariant clustering to return the same format as the rest of the clustering
            pass
        case ComparisonMethod.IS:
            cluster_result = cluster_reactions(list_reactions=reactions)
        case ComparisonMethod.WL:
            cluster_result = cluster_weisfeiler_lehman_nx(
                list_reactions=reactions,
                iterations=3,  # TODO: add property to enum
                use_edge_node_attr=False,  # TODO: add property to enum
            )
        case _: 
            raise NotImplementedError('Unknown comparison method.')
        
    return cluster_result

def save_results(cluster_result: ClusterDict) -> None:
    with open("clustered_reactions.pkl.gz", "wb") as file:
        pickle.dump(cluster_result, file)

    return None







