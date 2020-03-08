#include "mem/cache/replacement_policies/lip_rp.hh"

#include <cassert>
#include <memory>

#include "params/LIPRP.hh"

LIPRP::LIPRP(const Params *p)
    : BaseReplacementPolicy(p)
{
}

void
LIPRP::invalidate(const std::shared_ptr<ReplacementData>& replacement_data)
const
{
    std::static_pointer_cast<LIPReplData>(
        replacement_data)->priority = invalid_priority;
}

void
LIPRP::touch(const std::shared_ptr<ReplacementData>& replacement_data) const
{
    // Update cache entry to MRU position (last to be evicted).
    std::static_pointer_cast<LIPReplData>(
        replacement_data)->priority = ++positive_counter;
}

void
LIPRP::reset(const std::shared_ptr<ReplacementData>& replacement_data) const
{
    // Update cache entry to LRU position (lowest priority -- first to
    // be evicted).
    std::static_pointer_cast<LIPReplData>(
        replacement_data)->priority = -(++positive_counter);
}

ReplaceableEntry*
LIPRP::getVictim(const ReplacementCandidates& candidates) const
{
    // There must be at least one replacement candidate
    assert(candidates.size() > 0);

    // Visit all candidates to find victim
    ReplaceableEntry* victim = candidates[0];
    for (const auto& candidate : candidates) {
        // Update victim entry if necessary
        if (std::static_pointer_cast<LIPReplData>(
                    candidate->replacementData)->priority <
                std::static_pointer_cast<LIPReplData>(
                    victim->replacementData)->priority) {
            victim = candidate;
        }
    }

    // Prioritize invalid entries -- this is guaranteed by
    // invalid_priority being the lowest possible value.
    assert(invalid_priority == std::numeric_limits<decltype(invalid_priority)>::min());

    return victim;
}

std::shared_ptr<ReplacementData>
LIPRP::instantiateEntry()
{
    return std::shared_ptr<ReplacementData>(new LIPReplData());
}

LIPRP*
LIPRPParams::create()
{
    return new LIPRP(this);
}
