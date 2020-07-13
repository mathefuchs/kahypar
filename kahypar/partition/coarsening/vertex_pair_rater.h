/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <limits>
#include <stack>
#include <vector>

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/partition/coarsening/policies/fixed_vertex_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_community_policy.h"
#include "kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "kahypar/partition/coarsening/policies/rating_partition_policy.h"
#include "kahypar/partition/coarsening/policies/rating_score_policy.h"
#include "kahypar/partition/coarsening/policies/rating_tie_breaking_policy.h"
#include "kahypar/partition/context.h"

namespace kahypar {
    template<class ScorePolicy = HeavyEdgeScore,
            class HeavyNodePenaltyPolicy = MultiplicativePenalty,
            class CommunityPolicy = UseCommunityStructure,
            class RatingPartitionPolicy = NormalPartitionPolicy,
            class AcceptancePolicy = BestRatingWithTieBreaking<>,
            class FixedVertexPolicy = AllowFreeOnFixedFreeOnFreeFixedOnFixed,
            typename RatingType = RatingType>
    class VertexPairRater {
    private:
        static constexpr bool debug = false;

        class VertexPairRating {
        public:
            VertexPairRating(HypernodeID trgt, RatingType val, bool is_valid) :
                    target(trgt), value(val), valid(is_valid) {
            }

            VertexPairRating() :
                    target(std::numeric_limits<HypernodeID>::max()),
                    value(std::numeric_limits<RatingType>::min()),
                    valid(false) {}

            VertexPairRating(const VertexPairRating &) = delete;

            VertexPairRating &operator=(const VertexPairRating &) = delete;

            VertexPairRating(VertexPairRating &&) = default;

            VertexPairRating &operator=(VertexPairRating &&) = delete;

            ~VertexPairRating() = default;

            HypernodeID target;
            RatingType value;
            bool valid;
        };

    public:
        using Rating = VertexPairRating;

        VertexPairRater(Hypergraph &hypergraph, const Context &context) :
                _hg(hypergraph),
                _context(context),
                _tmp_ratings(_hg.initialNumNodes()),
                _tmp_num_common_incident_nets(_hg.initialNumNodes()),
                _already_matched(_hg.initialNumNodes()),
                _pin_already_visited(ScorePolicy::requires_flag_array ? hypergraph.initialNumNodes() : 1),
                _pin_already_visited_for_calc(hypergraph.initialNumNodes()),
                _avg_deg(0), _avg_node_weight(0) {

            if (ScorePolicy::requires_additional_fields) {
                HyperedgeID sumDeg = 0;
                HypernodeWeight sumWeight = 0;
                for (const auto node : hypergraph.nodes()) {
                    sumDeg += hypergraph.nodeDegree(node);
                    sumWeight += hypergraph.nodeWeight(node);
                }
                _avg_deg = static_cast<RatingType>(sumDeg) / static_cast<RatingType>(hypergraph.currentNumNodes());
                _avg_node_weight =
                        static_cast<RatingType>(sumWeight) / static_cast<RatingType>(hypergraph.currentNumNodes());
            }
        }

        VertexPairRater(const VertexPairRater &) = delete;

        VertexPairRater &operator=(const VertexPairRater &) = delete;

        VertexPairRater(VertexPairRater &&) = delete;

        VertexPairRater &operator=(VertexPairRater &&) = delete;

        ~VertexPairRater() = default;

        VertexPairRating rate(const HypernodeID u) {
            DBG << "Calculating rating for HN" << u;
            _pin_already_visited_for_calc.reset();
            const HypernodeWeight weight_u = _hg.nodeWeight(u);

            // Preload heavy edge score if necessary
            if (ScorePolicy::preload_edge_metrics) {
                for (const HyperedgeID &he : _hg.incidentEdges(u)) {
                    ASSERT(_hg.edgeSize(he) > 1, V(he));
                    if (_hg.edgeSize(he) <= _context.partition.hyperedge_size_threshold) {
                        const RatingType score = HeavyEdgeScore::calculate_rating(_hg, he);
                        for (const HypernodeID &v : _hg.pins(he)) {
                            if (v != u && belowThresholdNodeWeight(weight_u, _hg.nodeWeight(v)) &&
                                RatingPartitionPolicy::accept(_hg, _context, u, v)) {
                                _tmp_ratings[v] += score;
                                ++_tmp_num_common_incident_nets[v];
                            }
                        }
                    }
                }
            }

            // Process edge metrics and score node pairs
            for (const HyperedgeID &he : _hg.incidentEdges(u)) {
                ASSERT(_hg.edgeSize(he) > 1, V(he));
                if (_hg.edgeSize(he) <= _context.partition.hyperedge_size_threshold) {

                    // Load per edge score if required
                    RatingType score = 0.0;
                    if (!ScorePolicy::preload_edge_metrics) {
                        score = ScorePolicy::score(_hg, _context, he, u, -1, _pin_already_visited, _avg_deg,
                                                   _avg_node_weight, -1, 0);
                    }

                    for (const HypernodeID &v : _hg.pins(he)) {
                        if (v != u && !_pin_already_visited_for_calc[v] &&
                            belowThresholdNodeWeight(weight_u, _hg.nodeWeight(v)) &&
                            RatingPartitionPolicy::accept(_hg, _context, u, v)) {

                            // Whether to use per edge or per node scoring
                            if (ScorePolicy::preload_edge_metrics) {
                                _pin_already_visited_for_calc.set(v, true);
                                _tmp_ratings[v] = ScorePolicy::score(
                                        _hg, _context, he, u, v, _pin_already_visited, _avg_deg, _avg_node_weight,
                                        _tmp_ratings[v], _tmp_num_common_incident_nets[v]);
                            } else {
                                _tmp_ratings[v] += score;
                            }
                        }
                    }
                }
            }

            // Find max rating
            RatingType max_rating = std::numeric_limits<RatingType>::min();
            HypernodeID target = std::numeric_limits<HypernodeID>::max();
            for (auto it = _tmp_ratings.end() - 1; it >= _tmp_ratings.begin(); --it) {
                const HypernodeID tmp_target = it->key;
                const HypernodeWeight target_weight = _hg.nodeWeight(tmp_target);
                HypernodeWeight penalty = HeavyNodePenaltyPolicy::penalty(weight_u,
                                                                          target_weight);
                penalty = penalty == 0 ? std::max(std::max(weight_u, target_weight), 1) : penalty;
                const RatingType tmp_rating = it->value / static_cast<double>(penalty);
                DBG << "r(" << u << "," << tmp_target << ")=" << tmp_rating;
                if (CommunityPolicy::sameCommunity(_hg.communities(), u, tmp_target) &&
                    AcceptancePolicy::acceptRating(tmp_rating, max_rating,
                                                   target, tmp_target, _already_matched) &&
                    FixedVertexPolicy::acceptContraction(_hg, _context, u, tmp_target)) {
                    max_rating = tmp_rating;
                    target = tmp_target;
                }
            }

            // Return max rating if available
            VertexPairRating ret;
            if (max_rating != std::numeric_limits<RatingType>::min()) {
                ASSERT(target != std::numeric_limits<HypernodeID>::max(), "invalid contraction target");
                ret.value = max_rating;
                ret.target = target;
                ret.valid = true;
                ASSERT(_hg.communities()[u] == _hg.communities()[ret.target]);
            }
            ASSERT(!ret.valid || (_hg.partID(u) == _hg.partID(ret.target)));
            _tmp_ratings.clear();
            _tmp_num_common_incident_nets.clear();
            DBG << "rating=(" << ret.value << "," << ret.target << "," << ret.valid << ")";
            return ret;
        }

        void markAsMatched(const HypernodeID hn) {
            _already_matched.set(hn, true);
        }

        void resetMatches() {
            _already_matched.reset();
        }

        HypernodeWeight thresholdNodeWeight() const {
            return _context.coarsening.max_allowed_node_weight;
        }

    private:
        bool belowThresholdNodeWeight(const HypernodeWeight weight_u,
                                      const HypernodeWeight weight_v) const {
            return weight_v + weight_u <= _context.coarsening.max_allowed_node_weight;
        }

        Hypergraph &_hg;
        const Context &_context;
        ds::SparseMap<HypernodeID, RatingType> _tmp_ratings;
        ds::SparseMap<HypernodeID, size_t> _tmp_num_common_incident_nets;
        ds::FastResetFlagArray<> _already_matched;
        ds::FastResetFlagArray<> _pin_already_visited;
        ds::FastResetFlagArray<> _pin_already_visited_for_calc;
        ds::FastResetFlagArray<> _edge_already_visited;
        RatingType _avg_deg;
        RatingType _avg_node_weight;
    };
}  // namespace kahypar
