/***************************************************************************
 *  Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/
#include <boost/program_options.hpp>

#include <chrono>
#include <memory>
#include <string>

#include "lib/definitions.h"

#include "lib/core/Registrar.h"
#include "lib/io/HypergraphIO.h"
#include "lib/io/PartitioningOutput.h"
#include "lib/macros.h"
#include "partition/Configuration.h"
#include "partition/Factories.h"
#include "partition/Metrics.h"
#include "partition/Partitioner.h"
#include "tools/RandomFunctions.h"
#include "partition/initial_partitioning/RandomInitialPartitioner.h"
#include "partition/initial_partitioning/BFSInitialPartitioner.h"
#include "partition/initial_partitioning/GreedyHypergraphGrowingSequentialInitialPartitioner.h"
#include "partition/initial_partitioning/GreedyHypergraphGrowingGlobalInitialPartitioner.h"
#include "partition/initial_partitioning/GreedyHypergraphGrowingRoundRobinInitialPartitioner.h"
#include "partition/initial_partitioning/LabelPropagationInitialPartitioner.h"
#include "partition/initial_partitioning/nLevelInitialPartitioner.h"
#include "partition/initial_partitioning/PoolInitialPartitioner.h"
#include "partition/initial_partitioning/policies/StartNodeSelectionPolicy.h"
#include "partition/initial_partitioning/policies/GainComputationPolicy.h"
#include "partition/initial_partitioning/RecursiveBisection.h"
#include "partition/initial_partitioning/IInitialPartitioner.h"

#include "lib/serializer/SQLPlotToolsSerializer.h"

namespace po = boost::program_options;

using core::Registrar;
using partition::Partitioner;
using partition::InitialPartitioner;
using partition::Configuration;
using partition::Mode;
using partition::CoarseningAlgorithm;
using partition::RefinementAlgorithm;
using partition::InitialPartitionerAlgorithm;
using partition::RefinementStoppingRule;
using partition::CoarsenerFactory;
using partition::RefinerFactory;
using partition::InitialPartitioningFactory;
using partition::RandomWinsFullCoarsener;
using partition::RandomWinsLazyUpdateCoarsener;
using partition::RandomWinsHeuristicCoarsener;
using partition::HyperedgeCoarsener2;
using partition::TwoWayFMFactoryDispatcher;
using partition::HyperedgeFMFactoryDispatcher;
using partition::KWayFMFactoryDispatcher;
using partition::MaxGainNodeKWayFMFactoryDispatcher;
using partition::TwoWayFMFactoryExecutor;
using partition::HyperedgeFMFactoryExecutor;
using partition::KWayFMFactoryExecutor;
using partition::MaxGainNodeKWayFMFactoryExecutor;
using partition::LPRefiner;
using partition::IInitialPartitioner;
using partition::BFSInitialPartitioner;
using partition::LabelPropagationInitialPartitioner;
using partition::RandomInitialPartitioner;
using partition::GreedyHypergraphGrowingGlobalInitialPartitioner;
using partition::GreedyHypergraphGrowingRoundRobinInitialPartitioner;
using partition::GreedyHypergraphGrowingSequentialInitialPartitioner;
using partition::PoolInitialPartitioner;
using partition::nLevelInitialPartitioner;
using partition::RecursiveBisection;
using partition::FMGainComputationPolicy;
using partition::MaxPinGainComputationPolicy;
using partition::MaxNetGainComputationPolicy;
using partition::BFSStartNodeSelectionPolicy;
using serializer::SQLPlotToolsSerializer;
using defs::Hypergraph;
using defs::HypernodeID;
using defs::HypernodeWeight;
using defs::HyperedgeID;
using defs::HyperedgeIndexVector;
using defs::HyperedgeVector;
using defs::HyperedgeWeightVector;
using defs::HypernodeWeightVector;
using defs::HighResClockTimepoint;

InitialPartitionerAlgorithm stringToInitialPartitionerAlgorithm(
		std::string mode) {
	if (mode.compare("greedy") == 0) {
		return InitialPartitionerAlgorithm::greedy;
	} else if (mode.compare("greedy-global") == 0) {
		return InitialPartitionerAlgorithm::greedy_global;
	} else if (mode.compare("greedy-round") == 0) {
		return InitialPartitionerAlgorithm::greedy_round;
	} else if (mode.compare("greedy-maxpin") == 0) {
		return InitialPartitionerAlgorithm::greedy_maxpin;
	} else if (mode.compare("greedy-global-maxpin") == 0) {
		return InitialPartitionerAlgorithm::greedy_global_maxpin;
	} else if (mode.compare("greedy-round-maxpin") == 0) {
		return InitialPartitionerAlgorithm::greedy_round_maxpin;
	} else if (mode.compare("greedy-maxnet") == 0) {
		return InitialPartitionerAlgorithm::greedy_maxnet;
	} else if (mode.compare("greedy-global-maxnet") == 0) {
		return InitialPartitionerAlgorithm::greedy_global_maxnet;
	} else if (mode.compare("greedy-round-maxnet") == 0) {
		return InitialPartitionerAlgorithm::greedy_round_maxnet;
	} else if (mode.compare("recursive-greedy") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy;
	} else if (mode.compare("recursive-greedy-global") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_global;
	} else if (mode.compare("recursive-greedy-round") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_round;
	} else if (mode.compare("recursive-greedy-maxpin") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_maxpin;
	} else if (mode.compare("recursive-greedy-global-maxpin") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_global_maxpin;
	} else if (mode.compare("recursive-greedy-round-maxpin") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_round_maxpin;
	} else if (mode.compare("recursive-greedy-maxnet") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_maxnet;
	} else if (mode.compare("recursive-greedy-global-maxnet") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_global_maxnet;
	} else if (mode.compare("recursive-greedy-round-maxnet") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_round_maxnet;
	} else if (mode.compare("greedy-global-modify") == 0) {
		return InitialPartitionerAlgorithm::greedy_global_modify;
	} else if (mode.compare("recursive-greedy-global-modify") == 0) {
		return InitialPartitionerAlgorithm::rb_greedy_global_modify;
	} else if (mode.compare("lp") == 0) {
		return InitialPartitionerAlgorithm::lp;
	} else if (mode.compare("bfs") == 0) {
		return InitialPartitionerAlgorithm::bfs;
	} else if (mode.compare("random") == 0) {
		return InitialPartitionerAlgorithm::random;
	} else if (mode.compare("recursive-lp") == 0) {
		return InitialPartitionerAlgorithm::rb_lp;
	} else if (mode.compare("recursive-bfs") == 0) {
		return InitialPartitionerAlgorithm::rb_bfs;
	} else if (mode.compare("recursive-random") == 0) {
		return InitialPartitionerAlgorithm::rb_random;
	} else if (mode.compare("hMetis") == 0) {
		return InitialPartitionerAlgorithm::hMetis;
	} else if (mode.compare("PaToH") == 0) {
		return InitialPartitionerAlgorithm::PaToH;
	} else if (mode.compare("nLevel") == 0) {
		return InitialPartitionerAlgorithm::nLevel;
	} else if (mode.compare("direct_nLevel") == 0) {
		return InitialPartitionerAlgorithm::direct_nLevel;
	} else if (mode.compare("pool") == 0) {
		return InitialPartitionerAlgorithm::pool;
	} else if (mode.compare("ils") == 0) {
		return InitialPartitionerAlgorithm::ils;
	} else if (mode.compare("sa") == 0) {
		return InitialPartitionerAlgorithm::sa;
	}
	return InitialPartitionerAlgorithm::rb_greedy_global;
}

void configurePartitionerFromCommandLineInput(Configuration& config, const po::variables_map& vm) {
  if (vm.count("hgr") && vm.count("e") && vm.count("k")) {
    config.partition.graph_filename = vm["hgr"].as<std::string>();
    config.partition.k = vm["k"].as<PartitionID>();

    config.partition.coarse_graph_filename = std::string("/tmp/PID_")
                                             + std::to_string(getpid())
                                             + std::string("_coarse_")
                                             + config.partition.graph_filename.substr(
      config.partition.graph_filename.find_last_of("/")
      + 1);
    config.partition.graph_partition_filename = config.partition.graph_filename + ".part."
                                                + std::to_string(config.partition.k)
                                                + ".KaHyPar";
    config.partition.coarse_graph_partition_filename = config.partition.coarse_graph_filename
                                                       + ".part."
                                                       + std::to_string(config.partition.k);
    config.partition.epsilon = vm["e"].as<double>();

    if (vm.count("seed")) {
      config.partition.seed = vm["seed"].as<int>();
    }
    if (vm.count("nruns")) {
      config.partition.initial_partitioning_attempts = vm["nruns"].as<int>();
    }
    if (vm.count("part")) {
      if (vm["part"].as<std::string>() == "hMetis") {
        config.partition.initial_partitioner = InitialPartitioner::hMetis;
      } else if (vm["part"].as<std::string>() == "PaToH") {
        config.partition.initial_partitioner = InitialPartitioner::PaToH;
      } else if(vm["part"].as<std::string>() == "KaHyPar") {
    	  config.partition.initial_partitioner = InitialPartitioner::KaHyPar;
    	  if(vm.count("init-mode")) {
    		  config.initial_partitioning.mode = vm["init-mode"].as<std::string>();
    	  }
    	  if(vm.count("init-algo")) {
    		  config.initial_partitioning.algorithm = vm["init-algo"].as<std::string>();
    		  config.initial_partitioning.algo = stringToInitialPartitionerAlgorithm(config.initial_partitioning.algorithm);
    	  }
    	  if(vm.count("init-alpha")) {
    	     		  config.initial_partitioning.init_alpha = vm["init-alpha"].as<double>();
    	  }
      }
    }
    if (vm.count("mode")) {
      if (vm["mode"].as<std::string>() == "rb") {
        config.partition.mode = Mode::recursive_bisection;
      } else if (vm["mode"].as<std::string>() == "direct") {
        config.partition.mode = Mode::direct_kway;
      }
    }
    if (vm.count("part-path")) {
      config.partition.initial_partitioner_path = vm["part-path"].as<std::string>();
    } else {
      switch (config.partition.initial_partitioner) {
        case InitialPartitioner::hMetis:
          config.partition.initial_partitioner_path = "/software/hmetis-2.0pre1/Linux-x86_64/hmetis2.0pre1";
          break;
        case InitialPartitioner::PaToH:
          config.partition.initial_partitioner_path = "/software/patoh-Linux-x86_64/Linux-x86_64/patoh";
          break;
      }
    }
    if (vm.count("vcycles")) {
      config.partition.global_search_iterations = vm["vcycles"].as<int>();
    }
    if (vm.count("cmaxnet")) {
      config.partition.hyperedge_size_threshold = vm["cmaxnet"].as<HyperedgeID>();
      if (config.partition.hyperedge_size_threshold == -1) {
        config.partition.hyperedge_size_threshold = std::numeric_limits<HyperedgeID>::max();
      }
    }
    if (vm.count("ctype")) {
      if (vm["ctype"].as<std::string>() == "heavy_full") {
        config.partition.coarsening_algorithm = CoarseningAlgorithm::heavy_full;
      } else if (vm["ctype"].as<std::string>() == "heavy_partial") {
        config.partition.coarsening_algorithm = CoarseningAlgorithm::heavy_partial;
      } else if (vm["ctype"].as<std::string>() == "heavy_lazy") {
        config.partition.coarsening_algorithm = CoarseningAlgorithm::heavy_lazy;
      } else if (vm["ctype"].as<std::string>() == "hyperedge") {
        config.partition.coarsening_algorithm = CoarseningAlgorithm::hyperedge;
      } else {
        std::cout << "Illegal ctype option! Exiting..." << std::endl;
        exit(0);
      }
    }
    if (vm.count("s")) {
      config.coarsening.max_allowed_weight_multiplier = vm["s"].as<double>();
    }
    if (vm.count("t")) {
      config.coarsening.contraction_limit_multiplier = vm["t"].as<HypernodeID>();
    }
    if (vm.count("stopFM")) {
      if (vm["stopFM"].as<std::string>() == "simple") {
        config.fm_local_search.stopping_rule = RefinementStoppingRule::simple;
        config.her_fm.stopping_rule = RefinementStoppingRule::simple;
      } else if (vm["stopFM"].as<std::string>() == "adaptive1") {
        config.fm_local_search.stopping_rule = RefinementStoppingRule::adaptive1;
        config.her_fm.stopping_rule = RefinementStoppingRule::adaptive1;
      } else if (vm["stopFM"].as<std::string>() == "adaptive2") {
        config.fm_local_search.stopping_rule = RefinementStoppingRule::adaptive2;
        config.her_fm.stopping_rule = RefinementStoppingRule::adaptive2;
      } else {
        std::cout << "Illegal stopFM option! Exiting..." << std::endl;
        exit(0);
      }
    }
    if (vm.count("FMreps")) {
      config.fm_local_search.num_repetitions = vm["FMreps"].as<int>();
      config.her_fm.num_repetitions = vm["FMreps"].as<int>();
      if (config.fm_local_search.num_repetitions == -1) {
        config.fm_local_search.num_repetitions = std::numeric_limits<int>::max();
        config.her_fm.num_repetitions = std::numeric_limits<int>::max();
      }
    }
    if (vm.count("i")) {
      config.fm_local_search.max_number_of_fruitless_moves = vm["i"].as<int>();
      config.her_fm.max_number_of_fruitless_moves = vm["i"].as<int>();
    }
    if (vm.count("alpha")) {
      config.fm_local_search.alpha = vm["alpha"].as<double>();
      if (config.fm_local_search.alpha == -1) {
        config.fm_local_search.alpha = std::numeric_limits<double>::max();
      }
    }
    if (vm.count("verbose")) {
      config.partition.verbose_output = vm["verbose"].as<bool>();
    }
    if (vm.count("init-remove-hes")) {
      config.partition.initial_parallel_he_removal = vm["init-remove-hes"].as<bool>();
    }
    if (vm.count("lp_refiner_max_iterations")) {
      config.lp_refiner.max_number_iterations = vm["lp_refiner_max_iterations"].as<int>();
    }
    if (vm.count("rtype")) {
      if (vm["rtype"].as<std::string>() == "twoway_fm") {
        config.partition.refinement_algorithm = RefinementAlgorithm::twoway_fm;
      } else if (vm["rtype"].as<std::string>() == "kway_fm") {
        config.partition.refinement_algorithm = RefinementAlgorithm::kway_fm;
      } else if (vm["rtype"].as<std::string>() == "kway_fm_maxgain") {
        config.partition.refinement_algorithm = RefinementAlgorithm::kway_fm_maxgain;
      } else if (vm["rtype"].as<std::string>() == "hyperedge") {
        config.partition.refinement_algorithm = RefinementAlgorithm::hyperedge;
      } else if (vm["rtype"].as<std::string>() == "sclap") {
        config.partition.refinement_algorithm = RefinementAlgorithm::label_propagation;
      } else {
        std::cout << "Illegal stopFM option! Exiting..." << std::endl;
        exit(0);
      }
    }
  } else {
    std::cout << "Parameter error! Exiting..." << std::endl;
    exit(0);
  }
}

void setDefaults(Configuration& config) {
  config.partition.k = 2;
  config.partition.epsilon = 0.05;
  config.partition.seed = -1;
  config.partition.initial_partitioning_attempts = 1;
  config.partition.global_search_iterations = 1;
  config.partition.hyperedge_size_threshold = -1;
  config.partition.coarsening_algorithm = CoarseningAlgorithm::heavy_lazy;
  config.partition.refinement_algorithm = RefinementAlgorithm::kway_fm;
  config.coarsening.contraction_limit_multiplier = 160;
  config.coarsening.max_allowed_weight_multiplier = 2.5;
  config.coarsening.contraction_limit =
    config.coarsening.contraction_limit_multiplier * config.partition.k;
  config.coarsening.hypernode_weight_fraction = config.coarsening.max_allowed_weight_multiplier
                                                / config.coarsening.contraction_limit;
  config.fm_local_search.stopping_rule = RefinementStoppingRule::simple;
  config.fm_local_search.num_repetitions = -1;
  config.fm_local_search.max_number_of_fruitless_moves = 200;
  config.fm_local_search.alpha = 8;
  config.her_fm.stopping_rule = RefinementStoppingRule::simple;
  config.her_fm.num_repetitions = 1;
  config.her_fm.max_number_of_fruitless_moves = 10;
  config.lp_refiner.max_number_iterations = 3;
}

static Registrar<CoarsenerFactory> reg_heavy_lazy_coarsener(
  CoarseningAlgorithm::heavy_lazy,
  [](Hypergraph& hypergraph, const Configuration& config,
     const HypernodeWeight weight_of_heaviest_node)  -> ICoarsener* {
  return new RandomWinsLazyUpdateCoarsener(hypergraph, config, weight_of_heaviest_node);
});

static Registrar<CoarsenerFactory> reg_heavy_partial_coarsener(
  CoarseningAlgorithm::heavy_partial,
  [](Hypergraph& hypergraph, const Configuration& config,
     const HypernodeWeight weight_of_heaviest_node) -> ICoarsener* {
  return new RandomWinsHeuristicCoarsener(hypergraph, config,
                                          weight_of_heaviest_node);
});

static Registrar<CoarsenerFactory> reg_heavy_full_coarsener(
  CoarseningAlgorithm::heavy_full,
  [](Hypergraph& hypergraph, const Configuration& config,
     const HypernodeWeight weight_of_heaviest_node) -> ICoarsener* {
  return new RandomWinsFullCoarsener(hypergraph, config, weight_of_heaviest_node);
});

static Registrar<CoarsenerFactory> reg_hyperedge_coarsener(
  CoarseningAlgorithm::hyperedge,
  [](Hypergraph& hypergraph, const Configuration& config,
     const HypernodeWeight weight_of_heaviest_node) -> ICoarsener* {
  return new HyperedgeCoarsener2(hypergraph, config, weight_of_heaviest_node);
});

static Registrar<RefinerFactory> reg_twoway_fm_local_search(
  RefinementAlgorithm::twoway_fm,
  [](Hypergraph& hypergraph, const Configuration& config) {
  return TwoWayFMFactoryDispatcher::go(
    PolicyRegistry<RefinementStoppingRule>::getInstance().getPolicy(
      config.fm_local_search.stopping_rule),
    NullPolicy(),
    TwoWayFMFactoryExecutor(), hypergraph, config);
});

static Registrar<RefinerFactory> reg_kway_fm_maxgain_local_search(
  RefinementAlgorithm::kway_fm_maxgain,
  [](Hypergraph& hypergraph, const Configuration& config) {
  return MaxGainNodeKWayFMFactoryDispatcher::go(
    PolicyRegistry<RefinementStoppingRule>::getInstance().getPolicy(
      config.fm_local_search.stopping_rule),
    NullPolicy(),
    MaxGainNodeKWayFMFactoryExecutor(), hypergraph, config);
});

static Registrar<RefinerFactory> reg_kway_fm_local_search(
  RefinementAlgorithm::kway_fm,
  [](Hypergraph& hypergraph, const Configuration& config) {
  return KWayFMFactoryDispatcher::go(
    PolicyRegistry<RefinementStoppingRule>::getInstance().getPolicy(
      config.fm_local_search.stopping_rule),
    NullPolicy(),
    KWayFMFactoryExecutor(), hypergraph, config);
});

static Registrar<RefinerFactory> reg_lp_local_search(
  RefinementAlgorithm::label_propagation,
  [](Hypergraph& hypergraph, const Configuration& config) -> IRefiner* {
  return new LPRefiner(hypergraph, config);
});

static Registrar<RefinerFactory> reg_hyperedge_local_search(
  RefinementAlgorithm::hyperedge,
  [](Hypergraph& hypergraph, const Configuration& config) -> IRefiner* {
  return HyperedgeFMFactoryDispatcher::go(
    PolicyRegistry<RefinementStoppingRule>::getInstance().getPolicy(
      config.her_fm.stopping_rule),
    NullPolicy(),
    HyperedgeFMFactoryExecutor(), hypergraph, config);
});

static Registrar<PolicyRegistry<RefinementStoppingRule> > reg_simple_stopping(
  RefinementStoppingRule::simple, new NumberOfFruitlessMovesStopsSearch());

static Registrar<PolicyRegistry<RefinementStoppingRule> > reg_adaptive1_stopping(
  RefinementStoppingRule::adaptive1, new RandomWalkModelStopsSearch());

static Registrar<PolicyRegistry<RefinementStoppingRule> > reg_adaptive2_stopping(
  RefinementStoppingRule::adaptive2, new nGPRandomWalkStopsSearch());

static Registrar<InitialPartitioningFactory> reg_random(
		InitialPartitionerAlgorithm::random,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RandomInitialPartitioner(hypergraph,config);
		});

static Registrar<InitialPartitioningFactory> reg_bfs(
		InitialPartitionerAlgorithm::bfs,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new BFSInitialPartitioner<BFSStartNodeSelectionPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_lp(
		InitialPartitionerAlgorithm::lp,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new LabelPropagationInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy(
		InitialPartitionerAlgorithm::greedy,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingSequentialInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_global(
		InitialPartitionerAlgorithm::greedy_global,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingGlobalInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_round(
		InitialPartitionerAlgorithm::greedy_round,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingRoundRobinInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_maxpin(
		InitialPartitionerAlgorithm::greedy_maxpin,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingSequentialInitialPartitioner<BFSStartNodeSelectionPolicy,MaxPinGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_global_maxpin(
		InitialPartitionerAlgorithm::greedy_global_maxpin,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingGlobalInitialPartitioner<BFSStartNodeSelectionPolicy,MaxPinGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_round_maxpin(
		InitialPartitionerAlgorithm::greedy_round_maxpin,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingRoundRobinInitialPartitioner<BFSStartNodeSelectionPolicy,MaxPinGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_maxnet(
		InitialPartitionerAlgorithm::greedy_maxnet,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingSequentialInitialPartitioner<BFSStartNodeSelectionPolicy,MaxNetGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_global_maxnet(
		InitialPartitionerAlgorithm::greedy_global_maxnet,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingGlobalInitialPartitioner<BFSStartNodeSelectionPolicy,MaxNetGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_greedy_round_maxnet(
		InitialPartitionerAlgorithm::greedy_round_maxnet,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new GreedyHypergraphGrowingRoundRobinInitialPartitioner<BFSStartNodeSelectionPolicy,MaxNetGainComputationPolicy>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_random(
		InitialPartitionerAlgorithm::rb_random,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<RandomInitialPartitioner>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_bfs(
		InitialPartitionerAlgorithm::rb_bfs,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<BFSInitialPartitioner<BFSStartNodeSelectionPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy(
		InitialPartitionerAlgorithm::rb_greedy,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingSequentialInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_global(
		InitialPartitionerAlgorithm::rb_greedy_global,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingGlobalInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_round(
		InitialPartitionerAlgorithm::rb_greedy_round,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingRoundRobinInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_lp(
		InitialPartitionerAlgorithm::rb_lp,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<LabelPropagationInitialPartitioner<BFSStartNodeSelectionPolicy,FMGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_maxpin(
		InitialPartitionerAlgorithm::rb_greedy_maxpin,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingSequentialInitialPartitioner<BFSStartNodeSelectionPolicy,MaxPinGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_global_maxpin(
		InitialPartitionerAlgorithm::rb_greedy_global_maxpin,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingGlobalInitialPartitioner<BFSStartNodeSelectionPolicy,MaxPinGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_round_maxpin(
		InitialPartitionerAlgorithm::rb_greedy_round_maxpin,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingRoundRobinInitialPartitioner<BFSStartNodeSelectionPolicy,MaxPinGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_maxnet(
		InitialPartitionerAlgorithm::rb_greedy_maxnet,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingSequentialInitialPartitioner<BFSStartNodeSelectionPolicy,MaxNetGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_global_maxnet(
		InitialPartitionerAlgorithm::rb_greedy_global_maxnet,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingGlobalInitialPartitioner<BFSStartNodeSelectionPolicy,MaxNetGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_rb_greedy_round_maxnet(
		InitialPartitionerAlgorithm::rb_greedy_round_maxnet,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<GreedyHypergraphGrowingRoundRobinInitialPartitioner<BFSStartNodeSelectionPolicy,MaxNetGainComputationPolicy>>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_nlevel(
		InitialPartitionerAlgorithm::nLevel,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new RecursiveBisection<nLevelInitialPartitioner>(hypergraph,config);
		});
static Registrar<InitialPartitioningFactory> reg_pool(
		InitialPartitionerAlgorithm::pool,
		[](Hypergraph& hypergraph, Configuration& config) -> IInitialPartitioner* {
			return new PoolInitialPartitioner(hypergraph,config);
		});

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show help message")
    ("verbose", po::value<bool>(), "Verbose partitioner output")
    ("hgr", po::value<std::string>(), "Filename of the hypergraph to be partitioned")
    ("k", po::value<PartitionID>(), "Number of partitions")
    ("e", po::value<double>(), "Imbalance parameter epsilon")
    ("seed", po::value<int>(), "Seed for random number generator")
    ("mode", po::value<std::string>(), "(rb) recursive bisection, (direct) direct k-way")
    ("init-remove-hes", po::value<bool>(), "Initially remove parallel hyperedges before partitioning")
    ("nruns", po::value<int>(),
    "# initial partition trials, the final bisection corresponds to the one with the smallest cut")
    ("part", po::value<std::string>(), "Initial Partitioner: hMetis (default), PaToH or KaHyPar")
    ("init-mode", po::value<std::string>(), "If part=KaHyPar: direct or nLevel recursive-bisection Initial Partitioning")
    ("init-algo", po::value<std::string>(), "If part=KaHyPar, used initial partitioning algorithm")
    ("init-alpha", po::value<double>(), "Restrict initial partitioning epsilon to init-alpha*epsilon")
    ("part-path", po::value<std::string>(), "Path to Initial Partitioner Binary")
    ("vcycles", po::value<int>(), "# v-cycle iterations")
    ("cmaxnet", po::value<HyperedgeID>(), "Any hyperedges larger than cmaxnet are removed from the hypergraph before partition (disable:-1 (default))")
    ("ctype", po::value<std::string>(), "Coarsening: Scheme to be used: heavy_full (default), heavy_heuristic, heavy_lazy, hyperedge")
    ("s", po::value<double>(),
    "Coarsening: The maximum weight of a hypernode in the coarsest is:(s * w(Graph)) / (t * k)")
    ("t", po::value<HypernodeID>(), "Coarsening: Coarsening stops when there are no more than t * k hypernodes left")
    ("rtype", po::value<std::string>(), "Refinement: 2way_fm (default for k=2), her_fm, max_gain_kfm, kfm, lp_refiner")
    ("lp_refiner_max_iterations", po::value<int>(), "Refinement: maximum number of iterations for label propagation based refinement")
    ("stopFM", po::value<std::string>(), "2-Way-FM | HER-FM: Stopping rule \n adaptive1: new implementation based on nGP \n adaptive2: original nGP implementation \n simple: threshold based")
    ("FMreps", po::value<int>(), "2-Way-FM | HER-FM: max. # of local search repetitions on each level (default:1, no limit:-1)")
    ("i", po::value<int>(), "2-Way-FM | HER-FM: max. # fruitless moves before stopping local search (simple)")
    ("alpha", po::value<double>(), "2-Way-FM: Random Walk stop alpha (adaptive) (infinity: -1)")
    ("file", po::value<std::string>(), "filename of result file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  std::string result_file;
  if (vm.count("file")) {
    result_file = vm["file"].as<std::string>();
  }

  Configuration config;
  setDefaults(config);
  configurePartitionerFromCommandLineInput(config, vm);

  Randomize::setSeed(config.partition.seed);

  Hypergraph hypergraph(io::createHypergraphFromFile(config.partition.graph_filename,
                                                     config.partition.k));

  config.partition.total_graph_weight = hypergraph.totalWeight();

  config.coarsening.contraction_limit = config.coarsening.contraction_limit_multiplier
                                        * config.partition.k;

  config.coarsening.hypernode_weight_fraction = config.coarsening.max_allowed_weight_multiplier
                                                / config.coarsening.contraction_limit;

  config.partition.max_part_weight = (1 + config.partition.epsilon)
                                     * ceil(config.partition.total_graph_weight /
                                            static_cast<double>(config.partition.k));


  config.coarsening.max_allowed_node_weight = config.coarsening.hypernode_weight_fraction *
                                              config.partition.total_graph_weight;
  config.fm_local_search.beta = log(hypergraph.numNodes());

  // We use hMetis-RB as initial partitioner. If called to partition a graph into k parts
  // with an UBfactor of b, the maximal allowed partition size will be 0.5+(b/100)^(log2(k)) n.
  // In order to provide a balanced initial partitioning, we determine the UBfactor such that
  // the maximal allowed partiton size corresponds to our upper bound i.e.
  // (1+epsilon) * ceil(total_weight / k).
  double exp = 1.0 / log2(config.partition.k);
  config.partition.hmetis_ub_factor =
    50.0 * (2 * pow((1 + config.partition.epsilon), exp)
            * pow(ceil(static_cast<double>(config.partition.total_graph_weight)
                       / config.partition.k) / config.partition.total_graph_weight, exp) - 1);

  io::printPartitionerConfiguration(config);
  io::printHypergraphInfo(hypergraph, config.partition.graph_filename.substr(
                            config.partition.graph_filename.find_last_of("/") + 1));

#ifdef GATHER_STATS
  LOG("*******************************");
  LOG("***** GATHER_STATS ACTIVE *****");
  LOG("*******************************");
#endif

  Partitioner partitioner;

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  partitioner.partition(hypergraph, config);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

#ifndef NDEBUG
  for (const HyperedgeID he : hypergraph.edges()) {
    ASSERT([&]() -> bool {
      HypernodeID num_pins = 0;
      for (PartitionID i = 0; i < config.partition.k; ++i) {
        num_pins += hypergraph.pinCountInPart(he, i);
      }
      return num_pins == hypergraph.edgeSize(he);
    } (),
           "Incorrect calculation of pin counts");
  }
#endif

  std::chrono::duration<double> elapsed_seconds = end - start;

  io::printPartitioningStatistics();
  io::printPartitioningResults(hypergraph, elapsed_seconds, partitioner.timings());
  io::writePartitionFile(hypergraph, config.partition.graph_partition_filename);

  std::remove(config.partition.coarse_graph_filename.c_str());
  std::remove(config.partition.coarse_graph_partition_filename.c_str());

  SQLPlotToolsSerializer::serialize(config, hypergraph, partitioner,
                                    elapsed_seconds, partitioner.timings(), result_file);
  return 0;
}
