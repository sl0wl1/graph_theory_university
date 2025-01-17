from dataclasses import dataclass, field
from typing import Dict, List, Any, Literal

@dataclass
class Config:
    invariant: Literal["none", "rank", "algebraic_connectivity", "edge_count", "vertex_count", "vertex_degree"]
    algorithm: Literal["none", "isomorphism_test", "weisfeiler_lehmann_nx", "weisfeiler_lehmann_si", "neighbourhood"]
    weisfeiler_lehman_params: Dict[str, Any]
    neighbourhood_size: int = 3,

    def perform_benchmark(self) -> bool:
        if self.invariant == "algebraic_connectivity":
            return False
        if self.invariant == "none" and self.algorithm == "weisfeiler_lehmann_si":
            return False
        if self.invariant == "none" and self.algorithm == "none":
            return False
        return True 

# TYPE DEFINITIONS 
ReactionDict = Dict[str, Any]
ReactionDataList = List[ReactionDict]
ClusterDict = Dict[str, List[ReactionDict]]
BenchmarkingResult = Dict[str, Any]
ConfigList = List[Config]


