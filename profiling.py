import pickle
from itertools import product
import cProfile
import pstats
import io
import copy

# type definitions
from src.types_used import (
    ReactionDataList,
    ClusterDict,
    Config,
    BenchmarkingResult,
)

# custom imports
from src.rc_extract import get_rc_updated
from src.clustering import (
    cluster_without_invariant_grouping,
    cluster_after_invariant_grouping,
    group_after_invariant
)

def _all_invariants():
    return [
        "none",
        "rank",
        "algebraic_connectivity",
        "edge_count",
        # "vertex_count",
        # "vertex_degree",
    ]


def _all_algorithms():
    return ["none", "isomorphism_test", "weisfeiler_lehmann_nx", "weisfeiler_lehmann_si"]


def all_configurations():
    invariants = _all_invariants()
    algorithms = _all_algorithms()

    all_combinations = product(invariants, algorithms)
    configurations = [
        Config(invariant, algorithm) for invariant, algorithm in all_combinations
    ]

    return configurations


def load_data() -> ReactionDataList:
    with open("data/ITS_graphs.pkl.gz", "rb") as file:
        data = pickle.load(file)
    return data


def run_benchmarking():
    configurations = all_configurations()

    try:
        with open("data_with_reaction_centres.pkl.gz", "rb") as file:
            refined_data = pickle.load(file)
    except Exception as e:
        print(f"no saved data found: {e}")
        raw_data = load_data()
        refined_data = refine_data(raw_data)
        save_refined_data(refined_data)
    
    

    for config in configurations:
        if config.perform_benchmark():
            result = benchmark_configuration(config, refined_data)
            save_result(result)


def refine_data(raw_data: ReactionDataList) -> ReactionDataList:
    reactions = copy.deepcopy(raw_data)
    for reaction in reactions:
        reaction.pop("class", None)  # don't need the class property
        reaction_centre = get_rc_updated(reaction["ITS"])
        reaction["reaction_centre"] = reaction_centre
    return reactions


def save_refined_data(reactions: ReactionDataList):
    try:
        with open("data_with_reaction_centres.pkl.gz", "wb") as file:
            pickle.dump(reactions, file)
        return True
    except Exception as e:
        print(f"Error saving reaction centres: {e}")
        return False


def benchmark_configuration(
    config: Config, refined_data: ReactionDataList
) -> BenchmarkingResult:
    pr = cProfile.Profile()
    pr.enable()

    cluster_result = run_clustering(refined_data, config)

    pr.disable()
    pr.dump_stats(f"{config.invariant}_{config.algorithm}_profile.prof")

    # Print profiling stats
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats()
    print(s.getvalue())

    return {
        "clusters": cluster_result,
        "configuration": config,
        "time": pr.total_tt,
        "cluster_count": len(cluster_result),
    }


def run_clustering(
    config: Config,
    refined_data: ReactionDataList,
) -> ClusterDict:
    match config.invariant:
        case "none":
            cluster_result = cluster_without_invariant_grouping(
                list_reactions=refined_data, config=config
            )
        case (
            "rank"
            | "algebraic_connectivity"
            | "edge_count"
            | "vertex_count"
            | "vertex_degree"
        ):
            group_dict = group_after_invariant(
                list_reactions=refined_data, invariant=config.invariant
            )
            cluster_result = cluster_after_invariant_grouping(
                group_dict=group_dict, config=config
            )
        case _:
            raise ValueError(f"Invalid config: {config}")
        
    return cluster_result


def save_result(benchmark_result: BenchmarkingResult):
    try:
        config = benchmark_result["configuration"]
        namestr = f"{config.invariant}_{config.algorithm}"
        with open(f"{namestr}_results.pkl.gz", "wb") as file:
            pickle.dump(benchmark_result, file)

        return True
    except Exception as e:
        print(f"Error saving benchmarking results: {e}")
        return False

if __name__ == "__main__":
    run_benchmarking()