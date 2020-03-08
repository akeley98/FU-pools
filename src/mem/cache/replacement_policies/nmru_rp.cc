/**
 * Copyright (c) 2018 Inria
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Authors: Daniel Carvalho
 */

#include "mem/cache/replacement_policies/nmru_rp.hh"

#include <cstddef>
#include <cstdint>
#include <cassert>
#include <memory>

#include <iostream> /* XXX */
extern "C" { int akeley_verbose; }

#include "base/random.hh"
#include "params/NMRURP.hh"

NMRURP::NMRURP(const Params *p)
    : BaseReplacementPolicy(p)
{
}

void
NMRURP::invalidate(const std::shared_ptr<ReplacementData>& replacement_data)
const
{
    // Reset last touch timestamp
    std::static_pointer_cast<NMRUReplData>(
        replacement_data)->lastTouchTick = Tick(0);
}

void
NMRURP::touch(const std::shared_ptr<ReplacementData>& replacement_data) const
{
    // Update last touch timestamp
    std::static_pointer_cast<NMRUReplData>(
        replacement_data)->lastTouchTick = curTick();
}

void
NMRURP::reset(const std::shared_ptr<ReplacementData>& replacement_data) const
{
    // Set last touch timestamp
    std::static_pointer_cast<NMRUReplData>(
        replacement_data)->lastTouchTick = curTick();
}

ReplaceableEntry*
NMRURP::getVictim(const ReplacementCandidates& candidates) const
{
    // There must be at least one replacement candidate
    assert(candidates.size() > 0);

    // Degenerate case: only 1 candidate. Then we have no choice but
    // to evict the most-recently-used (and least-recently-used!) entry.
    if (candidates.size() == 1) return candidates[0];

    // Visit all candidates to find most recently used.
    // Prioritize replacing invalid candidates.
    ReplaceableEntry* mru = candidates[0];
    std::size_t mru_index = 0;
    std::size_t loop_index = 0;
    for (const auto& candidate : candidates) {
        auto candidate_tick = std::static_pointer_cast<NMRUReplData>(
            candidate->replacementData)->lastTouchTick;
        auto current_mru_tick = std::static_pointer_cast<NMRUReplData>(
            mru->replacementData)->lastTouchTick;

        // Return this candidate if invalid.
        if (!candidate_tick) return candidate;
        
        if (akeley_verbose) { // XXX
            std::cout << loop_index << ' ' << candidate_tick << '\n';
        }

        // Update mru entry if necessary
        if (candidate_tick > current_mru_tick) {
            mru = candidate;
            mru_index = loop_index;
        }
        ++loop_index;
    }

    // Of the candidates that are not the most-recently used, pick one
    // randomly to evict.
    auto max_index = unsigned(candidates.size() - 1);
    auto evict_index = random_mt.random<unsigned>(0, max_index - 1);
    if (evict_index >= mru_index) ++evict_index;
    
    if (akeley_verbose) { // XXX
        std::cout << "Evicting " << evict_index << '\n';
    }
    
    assert(evict_index < candidates.size() && evict_index != mru_index);
    return candidates[evict_index];
}

std::shared_ptr<ReplacementData>
NMRURP::instantiateEntry()
{
    return std::shared_ptr<ReplacementData>(new NMRUReplData());
}

NMRURP*
NMRURPParams::create()
{
    return new NMRURP(this);
}
