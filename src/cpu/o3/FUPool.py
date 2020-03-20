# Copyright (c) 2017 ARM Limited
# All rights reserved
#
# The license below extends only to copyright in the software and shall
# not be construed as granting a license to any other intellectual
# property including but not limited to intellectual property relating
# to a hardware implementation of the functionality of the software
# licensed hereunder.  You may use the software subject to the license
# terms below provided that you ensure that this notice is replicated
# unmodified and in its entirety in all distributions of the software,
# modified or unmodified, in source code or in binary form.
#
# Copyright (c) 2006-2007 The Regents of The University of Michigan
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer;
# redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution;
# neither the name of the copyright holders nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Authors: Kevin Lim

from m5.SimObject import SimObject
from m5.params import *
from m5.objects.FuncUnit import *
from m5.objects.FuncUnitConfig import *

class FUPool(SimObject):
    type = 'FUPool'
    cxx_header = "cpu/o3/fu_pool.hh"
    FUList = VectorParam.FUDesc("list of FU's for this pool")

class DefaultFUPool(FUPool):
    FUList = [ IntALU(), IntMultDiv(), FP_ALU(), FP_MultDiv(), ReadPort(),
               SIMD_Unit(), PredALU(), WritePort(), RdWrPort(), IprPort() ]

# Make a subclass of FUPool that contains the given number of each
# type of functional unit. For example,
#
# make_fu_pool_class(int_alu=4, int_mult_div=1)()
#
# Creates an FU Pool type consisting of 4 int ALUs and 1 int mult/div
# unit, and then instantiates one instance of said FU Pool.
def make_fu_pool_class(int_alu=0, int_mult_div=0, fp_alu=0,
                       fp_mult_div=0, read_port=0, simd_unit=0,
                       pred_alu=0, write_port=0, rdwr_port=0, ipr_port=0):
    unit_types = [ IntALU, IntMultDiv, FP_ALU, FP_MultDiv, ReadPort,
                   SIMD_Unit, PredALU, WritePort, RdWrPort, IprPort ]
    args = (int_alu, int_mult_div, fp_alu, fp_mult_div, read_port,
            simd_unit, pred_alu, write_port, rdwr_port, ipr_port)
    outer_fu_list = []
    for i, count in enumerate(args):
        if count <= 0: continue
        unit = unit_types[i]()
        unit.count = count
        outer_fu_list.append(unit)

    class GeneratedFUPool(FUPool):
        FUList = outer_fu_list

    return GeneratedFUPool
