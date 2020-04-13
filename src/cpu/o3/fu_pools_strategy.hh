/*
 * Copyright (c) 2020 David Zhao Akeley
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
 * Authors: David Zhao Akeley
 */

#ifndef __CPU_O3_FU_POOLS_STRATEGY_HH__
#define __CPU_O3_FU_POOLS_STRATEGY_HH__

#include "base/random.hh"
#include "base/types.hh"
#include "cpu/o3/fu_pool.hh"
#include "debug/IQ.hh"
#include "debug/IQFU.hh"
#include "debug/SIMDFU.hh"
#include "enums/OpClass.hh"
#include "params/DerivO3CPU.hh"

/** Status code used by the FU Pool selection function, to
 * indicate the reason for not scheduling an instruction on an FU
 * Pool where bypassing the needed operands is possible.
 */
enum class BypassStatus
{
    // Successfully scheduled inst on an FU Pool with bypassed values.
    Bypassed = 0,
    // Could not bypass because the FU Pool producing the needed
    // operand contains the needed FU type to execute this
    // instruction, but all of those FUs are busy.
    Congestion = 1,
    // Could not bypass because the FU Pool producing the needed
    // operand does not have an FU capable of executing this
    // instruction.
    Capability = 2,
    // Could not bypass because multiple different FU Pools
    // produced the needed operands for this instruction.
    Confluence = 3,
    // No bypassing was needed -- all operands are old enough.
    NotNeeded = 4
};

/** Return type of FU pool selection function.  */
struct ReserveResult
{
    // Reason for failure to bypass, if any.
    BypassStatus bypassStatus;
    // Index of the FU within the FU Pool that's chosen to execute
    // this instruction. NoFreeFU if we failed to choose an FU and
    // should try again later; NoCapableFU if we will execute the
    // instruction without an FU.
    int fuIdx;
    // Index/number of the FU Pool that will be chosen to execute
    // this instruction. Meaningless if fuIdx < 0.
    int poolIdx;
    // Latency of this instruction's execution on the chosen FU.
    Cycles opLatency;
};

inline bool isSimdOpClass(OpClass op_class)
{
    return (Enums::SimdAdd <= op_class)
        && (op_class <= Enums::SimdPredAlu);
}

/** Functional units are split into different pools. The FU Pools
 *  Strategy determines
 *
 *  1. How instructions requiring functional units to execute are
 *  assigned to different pools.
 *
 *  2. Whether and how values are bypassed between functional units.
 *
 *  3. What additional data is stored in each dynamic instruction in
 *     order to correctly implement the bypassing model.
 *
 *  This default strategy models multiple functional unit pools, with
 *  complete bypassing between FUs within a single pool, and no
 *  bypassing between different pools. If a value is ready on a
 *  certain FU pool on clock cycle C, then, for that clock cycle,
 *  ready instructions dependent on that value may only begin
 *  execution on said FU pool. This restriction is lifted on cycle C +
 *  1, at which point the produced value is considered visible to all
 *  FUs. This matches physical implementations that spend one cycle
 *  writing a result back to a renamed physical register file visibled
 *  to all functional units.
 *
 *  The InstructionRecord type, and reserveFU and
 *  InstructionRecord::resultOnFuPool functions should be replaced to
 *  customize the FUPoolStrategy.
 */
template <class Impl>
class DefaultFUPoolsStrategy
{
    typedef typename Impl::O3CPU O3CPU;
    typedef typename Impl::DynInstPtr DynInstPtr;

    /** Functional unit pools available to choose from. */
    std::vector<FUPool*> fuPools;

    /** Back pointer to CPU. */
    O3CPU *cpu;

    /** Integer that decides which fu pool selection algorithm is used. */
    unsigned fuPoolStrategy;

  public:
    DefaultFUPoolsStrategy(O3CPU *cpu_ptr, DerivO3CPUParams *params)
    {
        cpu = cpu_ptr;
        fuPools = cpu_ptr->fuPools;
        fuPoolStrategy = params->fuPoolStrategy;
        assert(fuPools.size() != 0);
    }

    /** Information that is kept in each dynamic instruction record,
     * which is used to keep track of which FU pools produce the
     * instruction's input values, and whether the instruction can be
     * executed with bypassing.
     */
    struct InstructionRecord
    {
        /** The cpu cycle number of the last call to resultOnFuPool */
        Cycles latestSrcRegCycle = Cycles(0);

        /** The FU pool index of the last call to resultOnFuPool,
         * except if there were multiple such calls within one CPU
         * cycle, and the FU pools aren't the same, this value is
         * -1. (This indicates that bypassing is not possible this
         * cycle as no one FU Pool has all the data needed for the
         * inst to execute.)
         */
        int srcRegFuPoolIndex = -1;

        /** Called whenever an FU pool produces a result needed by
         * this instruction. Requires the index of the FU Pool (in the
         * fuPools list) and the clock cycle that the result was
         * generated on.
         */
        void resultOnFuPool(int fu_pool_idx, Cycles current_cycle)
        {
            assert(fu_pool_idx >= 0);
            // New cycle? Older registers are now visible on all FUs, so
            // just record the FU producing the new result for this cycle.
            if (latestSrcRegCycle != current_cycle) {
                latestSrcRegCycle = current_cycle;
                srcRegFuPoolIndex = fu_pool_idx;
            }
            // Otherwise, if multiple results on different FU Pools are
            // available this cycle, record via srcRegFuPoolIndex that
            // bypassing is not possible this cycle.
            else if (srcRegFuPoolIndex != fu_pool_idx) {
                srcRegFuPoolIndex = -1;
            }
        }

        /** Given that the instruction is ready to issue, return
         * whether bypassing is needed to execute the instruction this
         * cycle (due to src reg results being too new). If so, write
         * through (*outFuPoolIndex) the number of the FU Pool that
         * this instruction may execute on (or -1, if no FU Pool can
         * execute the instruction).  If not, the instruction's needed
         * inputs are available at all FU Pools.
        */
        bool
        getBypassFuPoolIndex(int* out_fu_pool_idx, Cycles current_cycle) const
        {
            if (current_cycle != latestSrcRegCycle) return false;
            *out_fu_pool_idx = srcRegFuPoolIndex;
            return true;
        }
    };

    /** Given that this instruction has the given OpClass and is ready
     * to execute, select a suitable FU from an FU pool and mark it as
     * busy (unless no FU is needed).
     *
     * This function could use some improvement (in particular the
     * reliance on the "magic number" fuPoolStrategy parameter).
     */
    ReserveResult reserveFU(OpClass op_class, DynInstPtr issuing_inst)
    {
        ReserveResult result;
        result.poolIdx = -2;
        result.fuIdx = FUPool::NoCapableFU;
        result.opLatency = Cycles(1);

        bool debug_print = DTRACE(IQFU) ||
            (DTRACE(SIMDFU) && isSimdOpClass(op_class));

        // Change this if it turns out bypassing would be needed to
        // execute this instruction this cycle.
        result.bypassStatus = BypassStatus::NotNeeded;

        // Lambda that returns true iff the FU Pool with the given index
        // (in_pool_idx) is available for executing this instruction. If
        // so, reserve the FU and fill in the ReserveResult variables with
        // the correct values. Otherwise, set fuIdx to NoFreeFU only if
        // that's the reason (this has the effect of making fuIdx ==
        // NoCapableFU only when ALL FU pools are not capable of executing
        // this op class).
        const auto& fuPools_ = fuPools;
        auto try_fu =
            [&result, &fuPools_, op_class]
            (int in_pool_idx) -> bool
        {
            FUPool* fuPool = fuPools_[in_pool_idx];
            const int idx = fuPool->getUnit(op_class);
            if (idx >= 0) {
                result.fuIdx = idx;
                result.poolIdx = in_pool_idx;
                result.opLatency = fuPool->getOpLatency(op_class);
                return true;
            } else if (idx == FUPool::NoFreeFU) {
                result.fuIdx = idx;
            }
            return false;
        };

        bool need_bypassing = issuing_inst->fuPoolsRecord.getBypassFuPoolIndex(
            &result.poolIdx, cpu->cycleCounter);

        // Did we find a FU pool able to execute this instruction?
        bool success = false;
        assert(fuPools.size() > 0);

        // No FU needed case indicated by NoCapableFU (for some reason).
        if (op_class == No_OpClass) {
            // Nothing to do in this case.
        }

        // If bypassing would be needed, try to schedule this instruction
        // on the FU pool producing the bypassed value.
        else if (need_bypassing) {
            // Needed bypassed values on different FU pools -- fail.
            if (result.poolIdx < 0) {
                result.bypassStatus = BypassStatus::Confluence;
                result.fuIdx = FUPool::NoFreeFU;
            }
            // See if the FU pool producing the bypassed value is
            // ready to execute this instruction. If not for any
            // reason, fail, because maybe next cycle we'll find
            // another FU pool that is capable.
            else {
                try_fu(result.poolIdx);
                if (result.fuIdx == FUPool::NoFreeFU) {
                    result.bypassStatus = BypassStatus::Congestion;
                }
                else if (result.fuIdx == FUPool::NoCapableFU) {
                    result.bypassStatus = BypassStatus::Capability;
                    result.fuIdx = FUPool::NoFreeFU;
                }
                else {
                    success = true;
                    result.bypassStatus = BypassStatus::Bypassed;
                }
            }
            assert(result.bypassStatus != BypassStatus::NotNeeded);
        }
        // Below this point, bypassing is not needed, so all FU pools see
        // all the needed operands; we are free to choose any FU pool to
        // execute this instruction, using various strategies.

        // Greedy scheme.  Try FU Pools in reverse order until we find
        // one that can execute this instruction. (I don't remember
        // why I did this in reverse order but I don't want to change
        // it now and possibly mess up my previous sim results).
        else if (fuPoolStrategy <= 0) {
            for (int i = int(fuPools.size()) - 1; !success && i >= 0; --i) {
                success = try_fu(i);
            }
        }

        // Random scheme. Like the greedy scheme, except that we start probing
        // FU pools starting from a random index.
        else if (fuPoolStrategy == 1) {
            auto fu_pool_count = int(fuPools.size());
            int random_idx = random_mt.random<int>(0, fu_pool_count-1);
            for (int i = random_idx; !success && i >= 0; --i) {
                success = try_fu(i);
            }
            for (int i = fu_pool_count-1; !success && i > random_idx; --i) {
                success = try_fu(i);
            }
        }

        // Load balancing scheme. Select the FU pool that has the
        // highest number of FUs ready to execute this class of
        // instruction. (This is definitely a lazy 5 AM design; peek
        // behind the function call boundaries and you'll see this is
        // unreasonably expensive in simulation time).
        else {
            auto fu_pool_count = int(fuPools.size());
            int maxFreeUnitCount = FUPool::NoCapableFU;
            for (int i = 0; i < fu_pool_count; ++i) {
                int free_count = fuPools[i]->getFreeUnitCount(op_class);
                success |= (free_count >= 1);
                if (free_count > maxFreeUnitCount) {
                    maxFreeUnitCount = free_count;
                    result.poolIdx = i;
                }
                if (debug_print) {
                    DPRINTF_UNCONDITIONAL(
                        IQ, "op_class %d  pool %d  free %d\n",
                        int(op_class), i, free_count);
                }
            }
            if (success) {
                FUPool* fuPool = fuPools[result.poolIdx];
                result.fuIdx = fuPool->getUnit(op_class); // assert later
                result.opLatency = fuPool->getOpLatency(op_class);
                if (debug_print) {
                    DPRINTF_UNCONDITIONAL(
                        IQ, "Chose pool %d.\n", result.poolIdx);
                }
            }
            else {
                if (debug_print) {
                    DPRINTF_UNCONDITIONAL(IQ, "No pool chosen.\n");
                }
                result.fuIdx = maxFreeUnitCount == FUPool::NoCapableFU
                    ? FUPool::NoCapableFU
                    : FUPool::NoFreeFU;
            }
        }

        // Finally time to return the result (after doing a bunch of
        // debug stuff).
        if (success) {
            assert(result.poolIdx >= 0);
            assert(result.fuIdx >= 0);
            assert(result.poolIdx < int(fuPools.size()));
        }
        else {
            assert(
               result.fuIdx == FUPool::NoFreeFU ||
               result.fuIdx == FUPool::NoCapableFU);
        }

        if (debug_print) {
            std::stringstream ss;
            ss << "cycleCounter " << long(cpu->cycleCounter) << ":";
            const int src_reg_count = issuing_inst->numSrcRegs();
            for (int i = 0; i < src_reg_count; ++i) {
                ss << " r" << int(issuing_inst->renamedSrcRegIdx(i)->index());
            }
            ss << " ->";
            const int dest_reg_count = issuing_inst->numDestRegs();
            for (int i = 0; i < dest_reg_count; ++i) {
                ss << " r" << int(issuing_inst->renamedDestRegIdx(i)->index());
            }
            if (result.fuIdx == FUPool::NoCapableFU) {
                ss << (op_class == No_OpClass
                    ? "\n\t no FU needed" : "\n\t no capable FU");
            }
            else if (result.fuIdx == FUPool::NoFreeFU) {
                switch (result.bypassStatus) {
                    case BypassStatus::Congestion:
                        ss << "\n\t congestion bypass fail on FU pool ";
                        ss << result.poolIdx;
                        break;
                    case BypassStatus::Capability:
                        ss << "\n\t capability bypass fail on FU pool ";
                        ss << result.poolIdx;
                    break;
                case BypassStatus::Confluence:
                    ss << "\n\t confluence bypass fail";
                    break;
                case BypassStatus::NotNeeded:
                    ss << "\n\t no free FU";
                    break;
                case BypassStatus::Bypassed:
                    assert(0);
                default:
                    assert(0);
                }
            }
            else {
                assert(result.fuIdx >= 0);
                ss << "\n\t executing on FU pool ";
                ss << result.poolIdx;
                ss << (result.bypassStatus == BypassStatus::Bypassed
                    ? " (bypassed)" : " (no bypassing)");
            }
            DPRINTF_UNCONDITIONAL(IQ, "%s.\n", ss.str().c_str());
        }

        return result;
    }
};

#endif /* !__CPU_O3_FU_POOLS_STRATEGY_HH__ */
