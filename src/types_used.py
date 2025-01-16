from dataclasses import dataclass, field
from typing import Dict, List, Any, Literal

@dataclass
class Config:
    invariant: Literal["none", "rank", "connectivity", "edge_count", "vertex_count", "vertex_degree"]
    algorithm: Literal["none", "isomorphism_test", "weisfeiler_lehmann_nx", "weisfeiler_lehmann_si"]
    weisfeiler_lehman_params: Dict[str, Any] = field(default_factory=lambda: {"iterations": 3, "use_node_edge_attr": False})

    def perform_benchmark(self) -> bool:
        return not (self.algorithm == "none" and self.invariant == "none")
        

# TYPE DEFINITIONS 
ReactionDict = Dict[str, Any]
ReactionDataList = List[ReactionDict]
ClusterDict = Dict[str, List[ReactionDict]]
BenchmarkingResult = Dict[str, Any]
ConfigList = List[Config]