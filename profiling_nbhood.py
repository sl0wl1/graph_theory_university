import pickle
from itertools import product
import cProfile
import pstats
import io
import copy
from synutility.SynIO.data_type import load_from_pickle


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
    group_after_invariant,
)


def _all_invariants():
    return [
        "none",
        # "rank",
        # "algebraic_connectivity",
        # "edge_count",
        # "vertex_count",
        # "vertex_degree",
    ]


def _all_algorithms():
    return [
        # "none",
        "isomorphism_test",
        # "weisfeiler_lehmann_nx",
        # "weisfeiler_lehmann_si",
    ]


def _weisfeiler_lehmann_params():
    return {"iterations": 3, "use_node_edge_attr": True}


def all_configurations():
    invariants = _all_invariants()
    algorithms = _all_algorithms()

    all_combinations = product(invariants, algorithms)
    configurations = [
        Config(invariant, algorithm, _weisfeiler_lehmann_params())
        for invariant, algorithm in all_combinations
    ]

    return configurations


def load_data() -> ReactionDataList:
    data = load_from_pickle("data/ITS_graphs_reduced.pkl.gz")
    return data


def run_benchmarking():
    configurations = all_configurations()

    try:
        refined_data = load_from_pickle("data/data_with_reaction_centres.pkl.gz")
    except Exception as e:
        print(f"no saved data found: {e}")
        raw_data = load_data()
        refined_data = refine_data(raw_data)
        save_refined_data(refined_data)

    for config in configurations:
        if config.perform_benchmark():
            print(f"Running benchmark for {config}")
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
        with open("data/data_with_reaction_centres.pkl.gz", "wb") as file:
            pickle.dump(reactions, file)
        return True
    except Exception as e:
        print(f"Error saving reaction centres: {e}")
        return False


def benchmark_configuration(
    config: Config, refined_data: ReactionDataList
) -> BenchmarkingResult:

    refined_data = copy.deepcopy(refined_data)
    pr = cProfile.Profile()
    pr.enable()
    cluster_result = run_clustering(config=config, refined_data=refined_data)

    pr.disable()
    pr.dump_stats(f"data/profiles/{config.invariant}_{config.algorithm}_profile.prof")
    stats = pstats.Stats(pr)

    # Get total time
    total_time = stats.total_tt  # This gives you the total time in seconds

    benchmark_result = {
        "clusters": trim_cluster_result(cluster_result),
        "configuration": config,
        "time": total_time,
        "cluster_count": len(cluster_result),
    }

    print(
        f"Config: {config}, Cluster Count: {benchmark_result['cluster_count']}, Time: {benchmark_result['time']}"
    )
    print("\n\n")
    return benchmark_result


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
        with open(f"data/results/{namestr}_results.pkl.gz", "wb") as file:
            pickle.dump(benchmark_result, file)

        return True
    except Exception as e:
        print(f"Error saving benchmarking results: {e}")
        return False


def trim_cluster_result(cluster_result: ClusterDict) -> ReactionDataList:
    """
    Removes 'reaction_centre' and 'ITS' keys from each reaction in the list.
    """
    return {
        cluster_name: trim_reaction_data(cluster)
        for cluster_name, cluster in cluster_result.items()
    }


def trim_reaction_data(list_reactions: ReactionDataList) -> ReactionDataList:
    """
    Removes 'reaction_centre' and 'ITS' keys from each reaction in the list.
    """
    return [
        {
            key: value
            for key, value in reaction.items()
            if key not in ["reaction_centre", "ITS"]
        }
        for reaction in list_reactions
    ]


if __name__ == "__main__":
    run_benchmarking()
