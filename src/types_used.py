from dataclasses import dataclass
from typing import Dict, List, Any, Literal

@dataclass
class Config:
    invariant: Literal["none", "rank", "connectivity", "edge_count", "vertex_count", "vertex_degree"]
    algorithm: Literal["none", "isomorphism_test", "weisfeiler_lehmann_nx", "weisfeiler_lehmann_si"]
    weisfeiler_lehman_params: Dict[str, Any] = None # expects a dictionary with keys "iterations": Int and "use_node_edge_att": Boolean

    def perform_benchmark(self) -> bool:
        return self.algorithm != "none" and self.invariant != "none"
        

# TYPE DEFINITIONS 
ReactionDict = Dict[str, Any]
ReactionDataList = List[ReactionDict]
ClusterDict = Dict[str, List[ReactionDict]]
BenchmarkingResult = Dict[str, Any]
ConfigList = List[Config]


