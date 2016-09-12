/***************************************************************************
 *  Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "datastructure/FastResetBitVector.h"
#include "definitions.h"
#include "meta/Mandatory.h"
#include "meta/TemplateParameterToString.h"
#include "utils/Stats.h"
#include "partition/coarsening/ICoarsener.h"
#include "partition/coarsening/VertexPairCoarsenerBase.h"

using defs::Hypergraph;
using defs::HypernodeID;
using utils::Stats;
using datastructure::FastResetBitVector;

namespace partition {
template <class Rater = Mandatory>
class LazyVertexPairCoarsener final : public ICoarsener,
                                           private VertexPairCoarsenerBase<>{
 private:
  using Base = VertexPairCoarsenerBase;
  using Base::rateAllHypernodes;
  using Base::performContraction;
  using Rating = typename Rater::Rating;

  class NullMap {
 public:
    void insert(std::pair<HypernodeID, HypernodeID>) { }
  };

 public:
  LazyVertexPairCoarsener(Hypergraph& hypergraph, const Configuration& config,
                               const HypernodeWeight weight_of_heaviest_node) noexcept :
    Base(hypergraph, config, weight_of_heaviest_node),
    _rater(_hg, _config),
    _outdated_rating(hypergraph.initialNumNodes()),
    _target(_hg.initialNumNodes()) {
    LOG("Coarsener does not choose highest degree node as representative!");
    LOG("This could be slow on instances with skewed incidence structure.");
    LOG("Press any key to continue.");
  }

  virtual ~LazyVertexPairCoarsener() { }

  LazyVertexPairCoarsener(const LazyVertexPairCoarsener&) = delete;
  LazyVertexPairCoarsener& operator= (const LazyVertexPairCoarsener&) = delete;

  LazyVertexPairCoarsener(LazyVertexPairCoarsener&&) = delete;
  LazyVertexPairCoarsener& operator= (LazyVertexPairCoarsener&&) = delete;

 private:
  FRIEND_TEST(ALazyUpdateCoarsener, InvalidatesAdjacentHypernodesInsteadOfReratingThem);

  void coarsenImpl(const HypernodeID limit) noexcept override final {
    _pq.clear();

    NullMap null_map;
    rateAllHypernodes(_rater, _target, null_map);

    while (!_pq.empty() && _hg.currentNumNodes() > limit) {
      const HypernodeID rep_node = _pq.top();

      if (_outdated_rating[rep_node]) {
        DBG(dbg_coarsening_coarsen,
            "Rating for HN" << rep_node << " is invalid: " << _pq.topKey() << "--->"
            << _rater.rate(rep_node).value << " target=" << _rater.rate(rep_node).target
            << " valid= " << _rater.rate(rep_node).valid);
        updatePQandContractionTarget(rep_node, _rater.rate(rep_node));
      } else {
        const HypernodeID contracted_node = _target[rep_node];

        DBG(dbg_coarsening_coarsen, "Contracting: (" << rep_node << ","
            << _target[rep_node] << ") prio: " << _pq.topKey() <<
            " deg(" << rep_node << ")=" << _hg.nodeDegree(rep_node) <<
            " deg(" << contracted_node << ")=" << _hg.nodeDegree(contracted_node));

        ASSERT(_hg.nodeWeight(rep_node) + _hg.nodeWeight(_target[rep_node])
               <= _rater.thresholdNodeWeight(),
               "Trying to contract nodes violating maximum node weight");
        ASSERT(_pq.topKey() == _rater.rate(rep_node).value,
               "Key in PQ != rating calculated by rater:" << _pq.topKey() << "!="
               << _rater.rate(rep_node).value);

        performContraction(rep_node, contracted_node);
        ASSERT(_pq.contains(contracted_node), V(contracted_node));
        _pq.remove(contracted_node);

        // this also invalidates rep_node, however rep_node
        // will be re-rated and updated afterwards
        invalidateAffectedHypernodes(rep_node);

        updatePQandContractionTarget(rep_node, _rater.rate(rep_node));
      }
    }
  }

  bool uncoarsenImpl(IRefiner& refiner) noexcept override final {
    return Base::doUncoarsen(refiner);
  }

  std::string policyStringImpl() const noexcept override final {
    return std::string(" ratingFunction=" + templateToString<Rater>());
  }

  void invalidateAffectedHypernodes(const HypernodeID rep_node) noexcept {
    for (const HyperedgeID he : _hg.incidentEdges(rep_node)) {
      for (const HypernodeID pin : _hg.pins(he)) {
        _outdated_rating.set(pin, true);
      }
    }
  }

  void updatePQandContractionTarget(const HypernodeID hn, const Rating& rating) noexcept {
    _outdated_rating.set(hn, false);
    if (rating.valid) {
      ASSERT(_pq.contains(hn),
             "Trying to update rating of HN " << hn << " which is not in PQ");
      _pq.updateKey(hn, rating.value);
      _target[hn] = rating.target;
    } else {
      // In this case, no explicit contaiment check is necessary because the
      // method is only called on rep_node, which is definetly in the PQ.
      ASSERT(_pq.contains(hn),
             "Trying to remove rating of HN " << hn << " which is not in PQ");
      _pq.remove(hn);
      DBG(dbg_coarsening_no_valid_contraction, "Progress [" << _hg.currentNumNodes() << "/"
          << _hg.initialNumNodes() << "]:HN " << hn
          << " \t(w=" << _hg.nodeWeight(hn) << "," << " deg=" << _hg.nodeDegree(hn)
          << ") did not find valid contraction partner.");
#ifdef GATHER_STATS
      Stats::instance().add(_config, "numHNsWithoutValidContractionPartner", 1);
#endif
    }
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_config;
  using Base::_history;
  Rater _rater;
  FastResetBitVector<> _outdated_rating;
  std::vector<HypernodeID> _target;
};
}              // namespace partition