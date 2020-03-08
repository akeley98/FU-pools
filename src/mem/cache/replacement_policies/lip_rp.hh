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

/**
 * @file
 * Declaration of a LRU insertion policy replacement policy.
 */

#ifndef __MEM_CACHE_REPLACEMENT_POLICIES_LIP_RP_HH__
#define __MEM_CACHE_REPLACEMENT_POLICIES_LIP_RP_HH__

#include "mem/cache/replacement_policies/base.hh"
#include <cstdint>
#include <limits>

struct LIPRPParams;

class LIPRP : public BaseReplacementPolicy
{
  protected:
    static constexpr int64_t invalid_priority =
        std::numeric_limits<int64_t>::min();

    /** LIP-specific implementation of replacement data. */
    struct LIPReplData : ReplacementData
    {
        /** Eviction priority. Lowest is the first to be evicted. If
            == invalid_priority, this cache entry is not occupied. */
        std::int64_t priority = invalid_priority;
    };
    
    // Monotonically-increasing counter, used to fill in LIPReplData::priority.
    mutable int64_t positive_counter = 1;

  public:
    /** Convenience typedef. */
    typedef LIPRPParams Params;

    /**
     * Construct and initiliaze this replacement policy.
     */
    LIPRP(const Params *p);

    /**
     * Destructor.
     */
    ~LIPRP() {}

    /**
     * Invalidate replacement data to set it as the next probable victim.
     * Sets its last touch tick as the starting tick.
     *
     * @param replacement_data Replacement data to be invalidated.
     */
    void invalidate(const std::shared_ptr<ReplacementData>& replacement_data)
                                                              const override;

    /**
     * Touch an entry to update its replacement data.
     * Sets its last touch tick as the current tick.
     *
     * @param replacement_data Replacement data to be touched.
     */
    void touch(const std::shared_ptr<ReplacementData>& replacement_data) const
                                                                     override;

    /**
     * Reset replacement data. Used when an entry is inserted.
     * Sets its last touch tick as the current tick.
     *
     * @param replacement_data Replacement data to be reset.
     */
    void reset(const std::shared_ptr<ReplacementData>& replacement_data) const
                                                                     override;

    /**
     * Find replacement victim using LIP timestamps.
     *
     * @param candidates Replacement candidates, selected by indexing policy.
     * @return Replacement entry to be replaced.
     */
    ReplaceableEntry* getVictim(const ReplacementCandidates& candidates) const
                                                                     override;

    /**
     * Instantiate a replacement data entry.
     *
     * @return A shared pointer to the new replacement data.
     */
    std::shared_ptr<ReplacementData> instantiateEntry() override;
};

#endif // __MEM_CACHE_REPLACEMENT_POLICIES_LIP_RP_HH__
