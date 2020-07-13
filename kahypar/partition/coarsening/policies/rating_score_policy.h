/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/definitions.h"
#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"
#include "kahypar/partition/context.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"

#include <cmath>

namespace kahypar {
    class HeavyEdgeScore final : public meta::PolicyBase {
    public:
        static constexpr bool preload_edge_metrics = true;
        static constexpr bool requires_flag_array = false;
        static constexpr bool requires_additional_fields = false;

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType score(
                const Hypergraph &, const Context &, const HyperedgeID, const HypernodeID, const HypernodeID,
                ds::FastResetFlagArray<> &, const RatingType, const RatingType, const RatingType heavy_edge,
                const size_t) {
            return heavy_edge;
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType calculate_rating(
                const Hypergraph &hypergraph, const HyperedgeID he) {
            return static_cast<RatingType>(hypergraph.edgeWeight(he)) / (hypergraph.edgeSize(he) - 1);
        }
    };

    class EdgeFrequencyScore final : public meta::PolicyBase {
    public:
        static constexpr bool preload_edge_metrics = false;
        static constexpr bool requires_flag_array = false;
        static constexpr bool requires_additional_fields = false;

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType score(
                const Hypergraph &hypergraph, const Context &context, const HyperedgeID he, const HypernodeID,
                const HypernodeID, ds::FastResetFlagArray<> &, const RatingType, const RatingType, const RatingType,
                const size_t) {
            return static_cast<RatingType>(exp(-context.evolutionary.gamma * context.evolutionary.edge_frequency[he]) /
                                           hypergraph.edgeSize(he));
        }
    };

    class CombinedMetricScore final : public meta::PolicyBase {
    public:
        static constexpr bool preload_edge_metrics = true;
        static constexpr bool requires_flag_array = true;
        static constexpr bool requires_additional_fields = true;

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType score(
                const Hypergraph &hypergraph, const Context &context, const HyperedgeID, const HypernodeID u,
                const HypernodeID v, ds::FastResetFlagArray<> &pinVisitedBefore, const RatingType avg_deg,
                const RatingType avg_node_weight, const RatingType heavy_edge, const size_t common_incident_nets) {

            // Feature 17
            RatingType final_rating = weight_rating(scale_rating(average_deg_feat_17(hypergraph, u, v), feat17_index),
                                                    context, feat17_index);

            // Feature 20, 21, 22
            pinVisitedBefore.reset();
            size_t countNeighbors = 0;
            HyperedgeID sumNodeDegrees = 0;
            RatingType sumChiSquared = 0.0;

            for (const auto &incidentEdge : hypergraph.incidentEdges(u)) {
                for (const auto &neighbor : hypergraph.pins(incidentEdge)) {
                    if (!pinVisitedBefore[neighbor]) {
                        pinVisitedBefore.set(neighbor, true);

                        // Feature 20
                        ++countNeighbors;
                        sumNodeDegrees += hypergraph.nodeDegree(neighbor);

                        // Feature 21
                        const auto sub_deg = static_cast<RatingType>(hypergraph.nodeDegree(neighbor) - avg_deg);
                        sumChiSquared += (sub_deg * sub_deg) / static_cast<RatingType>(avg_deg);
                    }
                }
            }

            for (const auto &incidentEdge : hypergraph.incidentEdges(v)) {
                for (const auto &neighbor : hypergraph.pins(incidentEdge)) {
                    if (!pinVisitedBefore[neighbor]) {
                        pinVisitedBefore.set(neighbor, true);

                        // Feature 20
                        ++countNeighbors;
                        sumNodeDegrees += hypergraph.nodeDegree(neighbor);

                        // Feature 21
                        const auto sub_deg = static_cast<RatingType>(hypergraph.nodeDegree(neighbor) - avg_deg);
                        sumChiSquared += (sub_deg * sub_deg) / static_cast<RatingType>(avg_deg);
                    }
                }
            }

            // Feature 20
            const auto feat20 = static_cast<RatingType>(sumNodeDegrees) / static_cast<RatingType>(countNeighbors);
            final_rating += weight_rating(scale_rating(feat20, feat20_index), context, feat20_index);

            // Feature 21
            final_rating += weight_rating(scale_rating(sumChiSquared, feat21_index), context, feat21_index);

            // Feature 22
            const auto min_deg = static_cast<RatingType>(std::min(hypergraph.nodeDegree(u), hypergraph.nodeDegree(v)));
            const auto feat22 = static_cast<RatingType>(common_incident_nets) / min_deg
                                - static_cast<RatingType>(hypergraph.nodeWeight(u) + hypergraph.nodeWeight(v)) /
                                  static_cast<RatingType>(2 * avg_node_weight) + 1.0;
            final_rating += weight_rating(scale_rating(feat22, feat22_index), context, feat22_index);

            // Feature 24
            final_rating += weight_rating(scale_rating(strawman_feat_24(hypergraph, u, v, heavy_edge), feat24_index),
                                          context, feat24_index);

            return final_rating;
        }

    private:
        static constexpr size_t weight_table_size_num_metrics = 5; // Number of metrics used
        static constexpr size_t feat17_index = 0;
        static constexpr size_t feat20_index = 1;
        static constexpr size_t feat21_index = 2;
        static constexpr size_t feat22_index = 3;
        static constexpr size_t feat24_index = 4;
        // F17, F20, F21, F22, F24
        static constexpr RatingType weights[] = {0.134351, 0.018833, 0.169911, 0.365923, 0.310983, // k=2
                                                 0.167567, 0.079032, 0.164949, 0.374497, 0.213956, // k=4
                                                 0.203674, 0.215209, 0.165988, 0.362918, 0.052212, // k=8
                                                 0.237586, 0.288530, 0.176861, 0.357446, -0.060423}; // k>=16
        static constexpr RatingType boxcox_transform[] = {-0.19401318752826535, -0.022742996530458,
                                                          -0.08076547562124864, -0.32272110830116896,
                                                          -1155.864789069469};
        static constexpr RatingType mean_feat[] = {3.083458868248361, 2.983550498729181, 4.299104879799698,
                                                   0.38882348209716333, 0.018380497977313514};
        static constexpr RatingType std_feat[] = {1.5801378888179056, 1.2745601157855129, 4.03428153242448,
                                                  0.24413999690814012, 0.05575034968541556};

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType average_deg_feat_17(
                const Hypergraph &hypergraph, const HypernodeID u, const HypernodeID v) {
            return static_cast<RatingType>(hypergraph.nodeDegree(u) + hypergraph.nodeDegree(v)) / 2.0;
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType strawman_feat_24(
                const Hypergraph &hypergraph, const HypernodeID u, const HypernodeID v, const RatingType heavy_edge) {
            return static_cast<RatingType>(heavy_edge) /
                   static_cast<RatingType>((hypergraph.nodeDegree(u) - heavy_edge) *
                                           (hypergraph.nodeDegree(v) - heavy_edge));
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType scale_rating(
                const RatingType rating, const size_t feat_index) {
            const auto bc_transformed =
                    (std::pow(rating + 1.0, boxcox_transform[feat_index]) - 1.0) / boxcox_transform[feat_index];
            return static_cast<RatingType>((bc_transformed - mean_feat[feat_index]) / std_feat[feat_index]);
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType weight_rating(
                const RatingType scaled_rating, const Context &context, const size_t feat_index) {
            const size_t num_k_offset = getNumKOffset(context);
            return scaled_rating * weights[num_k_offset * weight_table_size_num_metrics + feat_index];
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline size_t getNumKOffset(const Context &context) {
            switch (context.partition.k) {
                case 2:
                    return 0;
                case 4:
                    return 1;
                case 8:
                    return 2;
                default:
                    return 3;
            }
        }
    };

    class SimpleCombinedMetricScore final : public meta::PolicyBase {
    public:
        static constexpr bool preload_edge_metrics = true;
        static constexpr bool requires_flag_array = false;
        static constexpr bool requires_additional_fields = true;

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType score(
                const Hypergraph &hypergraph, const Context &context, const HyperedgeID, const HypernodeID u,
                const HypernodeID v, ds::FastResetFlagArray<> &, const RatingType, const RatingType avg_node_weight,
                const RatingType heavy_edge, const size_t common_incident_nets) {

            // Feature 17
            RatingType final_rating = weight_rating(scale_rating(average_deg_feat_17(hypergraph, u, v), feat17_index),
                                                    context, feat17_index);

            // Feature 22
            const auto min_deg = static_cast<RatingType>(std::min(hypergraph.nodeDegree(u), hypergraph.nodeDegree(v)));
            const auto feat22 = static_cast<RatingType>(common_incident_nets) / min_deg
                                - static_cast<RatingType>(hypergraph.nodeWeight(u) + hypergraph.nodeWeight(v)) /
                                  static_cast<RatingType>(2 * avg_node_weight) + 1.0;
            final_rating += weight_rating(scale_rating(feat22, feat22_index), context, feat22_index);

            // Feature 24
            final_rating += weight_rating(scale_rating(strawman_feat_24(hypergraph, u, v, heavy_edge), feat24_index),
                                          context, feat24_index);

            return final_rating;
        }

    private:
        static constexpr size_t weight_table_size_num_metrics = 3; // Number of metrics used
        static constexpr size_t feat17_index = 0;
        static constexpr size_t feat22_index = 1;
        static constexpr size_t feat24_index = 2;
        // F17, F22, F24
        static constexpr RatingType weights[] = {0.165608171583601, 0.451056954702276, 0.383334873714123, // k=2
                                                 0.221643663573812, 0.495353535643722, 0.283002800782467, // k=4
                                                 0.329141939520741, 0.586483091031745, 0.084374969447513, // k=8
                                                 0.444411006261509, 0.66861101340892, -0.113022019670429}; // k>=16
        static constexpr RatingType boxcox_transform[] = {-0.19401318752826535, -0.32272110830116896,
                                                          -1155.864789069469};
        static constexpr RatingType mean_feat[] = {3.083458868248361, 0.38882348209716333, 0.018380497977313514};
        static constexpr RatingType std_feat[] = {1.5801378888179056, 0.24413999690814012, 0.05575034968541556};

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType average_deg_feat_17(
                const Hypergraph &hypergraph, const HypernodeID u, const HypernodeID v) {
            return static_cast<RatingType>(hypergraph.nodeDegree(u) + hypergraph.nodeDegree(v)) / 2.0;
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType strawman_feat_24(
                const Hypergraph &hypergraph, const HypernodeID u, const HypernodeID v, const RatingType heavy_edge) {
            return static_cast<RatingType>(heavy_edge) /
                   static_cast<RatingType>((hypergraph.nodeDegree(u) - heavy_edge) *
                                           (hypergraph.nodeDegree(v) - heavy_edge));
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType scale_rating(
                const RatingType rating, const size_t feat_index) {
            const auto bc_transformed =
                    (std::pow(rating + 1.0, boxcox_transform[feat_index]) - 1.0) / boxcox_transform[feat_index];
            return static_cast<RatingType>((bc_transformed - mean_feat[feat_index]) / std_feat[feat_index]);
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType weight_rating(
                const RatingType scaled_rating, const Context &context, const size_t feat_index) {
            const size_t num_k_offset = getNumKOffset(context);
            return scaled_rating * weights[num_k_offset * weight_table_size_num_metrics + feat_index];
        }

        KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline size_t getNumKOffset(const Context &context) {
            switch (context.partition.k) {
                case 2:
                    return 0;
                case 4:
                    return 1;
                case 8:
                    return 2;
                default:
                    return 3;
            }
        }
    };

    using RatingScorePolicies = meta::Typelist<HeavyEdgeScore, EdgeFrequencyScore, CombinedMetricScore, SimpleCombinedMetricScore>;
}  // namespace kahypar
