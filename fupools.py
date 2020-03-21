import math
import os
import os.path
import pathlib

max_processes = 16

get_default_args = lambda: [
    "configs/akeley/cache-se.py",
    "--cpu-type=DerivO3CPU",
    "--mem-type=DDR3_1600_8x8",
    "--caches",
    "--cpu-clock=4.0GHz",
    "--l1i_size=16kB",
    "--l1i_assoc=2",
    "--l1d_size=64kB",
    "--l1d_assoc=8",
    "--l2cache",
    "--l2_size=2048kB",
    "--l2_assoc=8",
    "--fu_pools=1" ]

benchmark_dir = "./akeley-workloads"
top_outdir = "./fu-pool-bench"
benchmark_executables = ["cholesky", "cnn", "particles", "3mm"]

cached_stats = {}

# Given benchmark name, look up its stats.txt file and return a
# dictionary of stat_name:stat_value pairs. Caches already loaded
# benchmarks.
def get_stats(benchmark_name):
    try:
        return cached_stats[benchmark_name]
    except KeyError:
        pass

    result_dict = {}
    path = os.path.join("%s/%s" % (top_outdir, benchmark_name), "stats.txt")
    for line in open(path):
        tokens = line.split()
        if len(tokens) >= 2:
            try:
                value = float(tokens[1])
            except ValueError:
                continue
            result_dict[tokens[0]] = value
    cached_stats[benchmark_name] = result_dict
    return result_dict


def fork_exec_gem5(args):
    """Fork and exec gem5 in child proc with the given arguments,
excluding argv[0] (automatically filled in with gem5 executable
path). Returns pid to wait for.
    """
    args = ["build/X86/gem5.opt"] + args
    outdir = None
    for a in args:
        if a.startswith("--outdir="): outdir = a[9:]

    print('\x1b[34m' + ' '.join(args) + '\x1b[0m')
    pid = os.fork()
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    if pid == 0:
        if outdir:
            flags = os.O_TRUNC | os.O_WRONLY | os.O_CREAT
            mode = 0o600
            stdout = os.open(os.path.join(outdir, "stdout.txt"), flags, mode)
            stderr = os.open(os.path.join(outdir, "stderr.txt"), flags, mode)
            os.dup2(stdout, 1)
            os.dup2(stderr, 2)
            os.close(stdout)
            os.close(stderr)
        os.execvp("build/X86/gem5.opt", args)
    else:
        return pid

# Dictionary mapping running command PIDs to their command line args.
pid_arg_dict = {}

def start_benchmarks(base_args, benchmark_name, benchmark_executables):
    """Given a list of base args for gem5, a name for the benchmark, and
    a list of benchmark executables (in the benchmark_dir directory),
    launch gem5 in parallel for all the named benchmark executables.
    Respects the max_process limit.

    Each gem5 run has two extra arguments in addition to base_args:
    one --outdir arg, based on the benchmark and executable name,
    and one --cmd arg, to specify one of the benchmark executables."""
    print("\x1b[36mRunning benchmark %r\x1b[0m" % (benchmark_name,))
    for f in benchmark_executables:
        wait_below_process_count(max_processes)
        benchmark_path = os.path.join(benchmark_dir, f)
        dirarg = "--outdir=%s/%s.%s" % (top_outdir, benchmark_name, f)
        args = [dirarg] + base_args + ["--cmd=%s" % benchmark_path]
        pid = fork_exec_gem5(args)
        pid_arg_dict[pid] = args

def wait_below_process_count(process_count):
    """Wait for child processes (recorded in pid_arg_dict). If there are
already (>= process_count) processes running, block until one finishes.

Throw exception if a child process terminates with nonzero status.
    """
    while len(pid_arg_dict) > 0:
        blocking = len(pid_arg_dict) >= process_count
        options = 0 if blocking else os.WNOHANG
        pid, status = os.waitpid(-1, options)
        if pid == 0: return
        if status != 0:
            raise Exception("Status %d from gem5 with args %r"
                % (status, pid_arg_dict[pid]))
        else:
            print("\x1b[32mCompleted gem5 with args %r\x1b[0m"
                % (pid_arg_dict[pid],))
        del pid_arg_dict[pid]

def replace_arg(args, new_arg):
    """Given a list of args, and a new_arg in the form --key=value, replace
    an arg in args that starts with --key= with the new_arg."""
    prefix = new_arg.split('=')[0] + "="
    for i, s in enumerate(args):
        if s.startswith(prefix):
            args[i] = new_arg
            return
    raise ValueError("No arg starts with %r" % prefix)

def baseline_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    start_benchmarks(args, "baseline", exe)

def more_mem_1_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    start_benchmarks(args, "mm-1", exe)

def greedy_2_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    replace_arg(args, "--fu_pools=2")
    args.append("--fu_pool_strategy=greedy")
    start_benchmarks(args, "greedy2", exe)

def more_mem_greedy_2_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    replace_arg(args, "--fu_pools=2")
    args.append("--fu_pool_strategy=greedy")
    start_benchmarks(args, "mm-greedy2", exe)

def greedy_4_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    replace_arg(args, "--fu_pools=4")
    args.append("--fu_pool_strategy=greedy")
    start_benchmarks(args, "greedy4", exe)

def more_mem_greedy_4_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    replace_arg(args, "--fu_pools=4")
    args.append("--fu_pool_strategy=greedy")
    start_benchmarks(args, "mm-greedy4", exe)

def random_2_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    replace_arg(args, "--fu_pools=2")
    args.append("--fu_pool_strategy=random")
    start_benchmarks(args, "random2", exe)

def more_mem_random_2_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    replace_arg(args, "--fu_pools=2")
    args.append("--fu_pool_strategy=random")
    start_benchmarks(args, "mm-random2", exe)

def random_4_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    replace_arg(args, "--fu_pools=4")
    args.append("--fu_pool_strategy=random")
    start_benchmarks(args, "random4", exe)

def more_mem_random_4_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    replace_arg(args, "--fu_pools=4")
    args.append("--fu_pool_strategy=random")
    start_benchmarks(args, "mm-random4", exe)

def balance_2_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    replace_arg(args, "--fu_pools=2")
    args.append("--fu_pool_strategy=balance")
    start_benchmarks(args, "balance2", exe)

def more_mem_balance_2_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    replace_arg(args, "--fu_pools=2")
    args.append("--fu_pool_strategy=balance")
    start_benchmarks(args, "mm-balance2", exe)

def balance_4_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    replace_arg(args, "--fu_pools=4")
    args.append("--fu_pool_strategy=balance")
    start_benchmarks(args, "balance4", exe)

def more_mem_balance_4_pool_benchmarks(exe=None):
    if exe is None: exe = benchmark_executables
    args = get_default_args()
    args.append("--fu_pool_more_mem=1")
    replace_arg(args, "--fu_pools=4")
    args.append("--fu_pool_strategy=balance")
    start_benchmarks(args, "mm-balance4", exe)

# Run this file with python3 -i and call main() to start the simulations.
def main(j=None):
    if j is not None:
        global max_processes
        max_processes = j

    baseline_benchmarks()

    greedy_2_pool_benchmarks()
    greedy_4_pool_benchmarks()
    random_2_pool_benchmarks()
    random_4_pool_benchmarks()
    balance_2_pool_benchmarks()
    balance_4_pool_benchmarks()

    # more_mem_1_pool_benchmarks()
    # more_mem_greedy_2_pool_benchmarks()
    # more_mem_greedy_4_pool_benchmarks()
    # more_mem_random_2_pool_benchmarks()
    # more_mem_random_4_pool_benchmarks()
    # more_mem_balance_2_pool_benchmarks()
    # more_mem_balance_4_pool_benchmarks()

    # Wait for all gem5 child processes to finish.
    wait_below_process_count(1)

# Given a benchmark_category (identifies the kind of CPU simulated;
# the part before . in a benchmark name), return the average speedup
# the benchmark had over thebaseline for the given benchmark
# executable, based on the final_tick stat (total picoseconds?).
def speed_over_baseline(benchmark_category, executable_name, baseline_config):
    if executable_name is not None:
        benchmark_name = benchmark_category + '.' + executable_name
        baseline_name = baseline_config + '.' + executable_name
        return (get_stats(baseline_name)["final_tick"] /
                get_stats(benchmark_name)["final_tick"])
    else:
        total = 0.0
        for executable_name in benchmark_executables:
            speedup = speed_over_baseline(benchmark_category, executable_name)
            total += speedup
        return total / len(benchmark_executables)

# For a given configuration and workload (e.g. 'greedy2' and 'cnn'),
# return latex for 5 columns of a table:
#
# * Percent slowdown from baseline (perfect bypassing)
# * Percentage of instructions executed with bypassing
# * Percentage points of said stat below baseline
# * Number of congestion bypass failures, in millions
# * Number of confluence bypass failures, in millions.
def stats_row(config, workload, baseline_config='baseline'):
    # Chandler Bing: Could there BE more escape characters???
    format = "%.2f\\%% & %.2f\\%% & %.2f\\%% & %.2f & %.2f \\\\"
    stats = get_stats(config + '.' + workload)
    baseline_stats = get_stats(baseline_config + '.' + workload)
    numer = stats["system.cpu.iq.insts_with_bypassing"]
    denom = stats["system.cpu.iq.insts_without_bypassing"] + numer
    percent_bypassed = 100 * numer / denom
    numer = baseline_stats["system.cpu.iq.insts_with_bypassing"]
    denom = baseline_stats["system.cpu.iq.insts_without_bypassing"] + numer
    delta_percent = 100 * numer / denom - percent_bypassed

    return format % (
        100 * (1-speed_over_baseline(config, workload, baseline_config)),
        percent_bypassed,
        delta_percent,
        stats["system.cpu.iq.congestion_bypass_fails"] * 1e-6,
        stats["system.cpu.iq.confluence_bypass_fails"] * 1e-6,
    )

# Generate a table of stats rows for every system configuration.
def generate_table_for_workload(workload):
    lines = []
    lines.append(r"\begin{tabular}{|r c|r r r r r|}")
    lines.append(r"\hline")
    lines.append(r"& & \% slowdown & \% insts & \(-\Delta\)\% points & congestion & confluence \\")
    lines.append(r"strategy & FU pools & from complete & bypassed & bypassed & (million insts) & (million insts) \\")
    lines.append(r"\hline")

    lines.append(" complete & 1 & " + stats_row("baseline", workload))
    for pools in (2,4):
        lines.append(r"\hline")
        for strat in ("greedy", "random", "balance"):
            row = "%s & %i & " % (strat, pools)
            row += stats_row(strat+str(pools), workload)
            lines.append(row)

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append("")
    return '\n'.join(lines)

def write_tex_tables():
    for workload in benchmark_executables:
        file_name = "./%s-fu-pool-table.tex" % workload
        open(file_name, 'w').write(generate_table_for_workload(workload))

# Generate a table of stats rows for every system configuration
# (variation with higher number of memory units).
def mm_generate_table_for_workload(workload):
    lines = []
    lines.append(r"\begin{tabular}{|r c|r r r r r|}")
    lines.append(r"\hline")
    lines.append(r"& & \% slowdown & \% insts & \(-\Delta\)\% points & congestion & confluence \\")
    lines.append(r"strategy & FU pools & from complete & bypassed & bypassed & (million insts) & (million insts) \\")
    lines.append(r"\hline")

    lines.append(" complete & 1 & " + stats_row("baseline", workload))
    lines.append(" mm-complete & 1 & " + stats_row("mm-1", workload, "mm-1"))
    for pools in (2,4):
        lines.append(r"\hline")
        for strat in ("greedy", "random", "balance"):
            row = "%s & %i & " % (strat, pools)
            row += stats_row(strat+str(pools), workload)
            lines.append(row)
        lines.append(r"\hline")
        for strat in ("greedy", "random", "balance"):
            row = "%s & %i & " % ("mm-" + strat, pools)
            row += stats_row("mm-"+strat+str(pools), workload, "mm-1")
            lines.append(row)

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append("")
    return '\n'.join(lines)

def mm_write_tex_tables():
    for workload in benchmark_executables:
        file_name = "./mm-%s-fu-pool-table.tex" % workload
        open(file_name, 'w').write(mm_generate_table_for_workload(workload))
