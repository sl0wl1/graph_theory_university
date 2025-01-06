from typing import Dict, List, Any
import networkx as nx
import networkx.algorithms.isomorphism as iso

from src.rc_extract import get_rc

def cluster_reactions(list_reactions: List[Dict[Any, Any]]) -> Dict[str, Any]:
    """Simple function for clusterting chemical reactions

    Args:
        list_reactions (List[Dict[Any, Any]]): A list of reactions

    Returns:
        Dict[str, Any]: Returns a dict. Keys are the number of the cluster. Values are the isomorphic reactions. 
    """
    
    # Create callables for is_isomorphic check
    nm_charge = iso.numerical_node_match("charge", 0)
    nm_element = iso.categorical_node_match("element", "C")
    em_order = iso.categorical_edge_match("order", 0)
    
    # Create an empty dict for storing reaction clusters and cluster counter for key naming
    cluster_dict = {}
    cluster_counter = 0
    
    # Check if the list of reactions is empty
    if list_reactions:

        for idx, reaction in enumerate(list_reactions):
            
            reaction_centre = get_rc(reaction["ITS"])
            
            # Create first entry in dict. For the first reaction there is nothing to compare 
            if idx == 0:
                cluster_dict[f"cluster_{cluster_counter}"] = [list_reactions[idx]]
                cluster_counter += 1
                
            else:  
                
                # Checks if isomorphs of the reaction centre already exist in a cluster
                for key, value in cluster_dict.items():
                    
                    cluster_centre = get_rc(value[0]["ITS"])
                    
                    #TODO Maybe this could be optimized... and I do not know if it is correct -.-'
                    # Only need to be checked once
                    if nx.is_isomorphic(cluster_centre, reaction_centre, node_match=nm_charge) and \
                    nx.is_isomorphic(cluster_centre, reaction_centre, node_match=nm_element) and \
                    nx.is_isomorphic(cluster_centre, reaction_centre, edge_match=em_order):
                        
                        value.append(reaction)
                        cluster_dict[key] = value
                        break
                    
                else:
                    
                    # If no isomorphic reaction centre can be found, create a new entry (cluster)
                    cluster_dict[f"cluster_{cluster_counter}"] = [reaction]
                    cluster_counter += 1                
                 
    return cluster_dict
        