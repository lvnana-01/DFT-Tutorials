# 第18章 自动化、并行性能与报错排查（Automation, Parallel Performance & Troubleshooting）

本章是 Quantum ESPRESSO 部分的收尾，讲三件事：第一，并行层级（`-ni` images、`-nk` pools、`-nb` band groups、`-nt` task groups、`-nd` 线性代数组）、MPI + OpenMP 混合并行、GPU 运行、内存估计（`Estimated max dynamical RAM`）与 `disk_io`，以及如何根据时间报告（`c_bands`、`cdiaghg`、`fft_scatter` 等）调整参数；第二，自动化：Slurm/PBS 作业模板，`max_seconds` + `restart_mode='restart'` + `prefix.EXIT` 的正常停止与续算，ASE 的 `EspressoProfile`、pymatgen 的 `PWInput`/`PWOutput`、qe-tools、qeschema 读 `data-file-schema.xml`，AiiDA 的 `PwBaseWorkChain`、aiidalab-qe 和 Materials Cloud 输入生成器；第三，一张"真实报错 → 原因 → 修复"的速查表，汇总全系列的常见问题。

> [!NOTE]
> 本章示例在 [examples/18-automation-parallel](../examples/18-automation-parallel/README.md)。仓库**没有运行**这些计算，输出片段均为"（示意输出，数值仅供参考）"。报错信息原文是从 QE 7.6 源码（`errore` 调用）和官方 PWscf 用户指南中摘出来的；不同版本的措辞可能略有差别。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 并行层级 | 理解 images / pools / 平面波（FFT） / task groups / 线性代数组各自分摊什么、能加速什么 |
| 选参数 | 根据 k 点数 N_k、FFT 第三维 N₃ 和能带数 M 快速给出 `-nk`、`-nt`、`-nd`；理解 7.1 起的自动猜测 |
| OpenMP 与 GPU | 设置 `OMP_NUM_THREADS` 和绑定；GPU 上每块卡一个 MPI 进程、优先 pool 并行 |
| 内存与 I/O | 读 `Estimated max dynamical RAM`，用 `disk_io`、`-nd`、减少 pools 来省内存；scratch 目录放在哪里 |
| 作业与续算 | Slurm/PBS 模板；`max_seconds`、`prefix.EXIT`、信号；`restart_mode='restart'` 的限制 |
| 自动化 | ASE、pymatgen、qe-tools、qeschema、AiiDA、Materials Cloud 工具各适合做什么 |
| 排错 | 速查表：namelist 读取、赝势、SCF 不收敛、cdiaghg、对称性、并行、ph.x、DFT+U 旧语法、杂化泛函等 |

与其它代码的对应关系（帮助从 VASP/CP2K 迁移）：

| 概念 | QE | VASP | CP2K |
|---|---|---|---|
| k 点并行 | `-nk`（pools） | `KPAR` | 无（一般只用 Γ 点） |
| 能带并行 | `-nb`（band groups，主要用于杂化泛函） | `NCORE`/`NPAR` | — |
| FFT 并行 | 默认（平面波/实空间分割）+ `-nt` | `NCORE` | 默认 |
| 线性代数 | `-nd`（ScaLAPACK/ELPA） | 自动（ScaLAPACK） | `&GLOBAL PREFERRED_DIAG_LIBRARY` 等 |
| 图像并行（NEB/声子） | `-ni` | `IMAGES` | `&REPLICA` / `NPROC_REP` |
| 正常停止 | `prefix.EXIT`、`max_seconds` | `STOPCAR` | `EXIT` 文件、`WALLTIME` |
| 续算 | `restart_mode='restart'` | `ISTART`/`ICHARG` + `CONTCAR` | `&EXT_RESTART` |

VASP 的后处理与排错见 [VASP 第18章](../../VASP/tutorials/18-postprocess-troubleshooting.md)，CP2K 见 [CP2K 第18章](../../CP2K/tutorials/18-troubleshooting.md)。作业调度系统的基础见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)，ASE/pymatgen 的基础见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)，通用工作流工具见 [Common 第6章](../../Common/tutorials/06-workflow-automation.md)。

---

## 2. 原理速览：QE 的并行层级

### 2.1 并行层级

QE 总用户手册把 MPI 进程按层级划分（命令行选项及缩写）：

| 层级 | 选项 | 缩写 | 分摊的对象 | 通信量 | 是否分摊内存 |
|---|---|---|---|---|---|
| images | `-nimage` | `-ni` | NEB 的各个像、ph.x 的不可约 q 点/表示 | 很少 | 否 |
| k 点 pools | `-npools` | `-nk` | k 点（及自旋） | 少 | 只分摊波函数，其余不分摊 |
| band groups | `-nband` | `-nb` | 能带（主要用于杂化泛函的 EXX） | 中等 | 否（用户指南：只有 FFT 与线性代数并行分摊内存） |
| 平面波 / FFT | （默认，pool 内部） | — | G 矢量与实空间网格（沿 z 分平面） | 很大 | 是（几乎线性下降） |
| task groups | `-ntg` | `-nt` | 同时处理多个波函数的 FFT | 中等 | 否 |
| 线性代数组 | `-ndiag` / `-northo` | `-nd` | 子空间对角化/正交化矩阵（方形网格 n²） | 中等 | 分摊 M×M 矩阵 |

官方示例：

```bash
mpirun -np 4096 ./neb.x -ni 8 -nk 2 -nt 4 -nd 144 -i my.input
```

这条命令的含义是：8 个 image，每个 512 进程；每个 image 分成 2 个 pool，每个 pool 256 进程；FFT 用 4 个 task group（每组 64 进程）；子空间对角化分布在 12×12 = 144 个进程上。

官方手册给出的默认值：`-ni 1 -nk 1 -nt 1`；`-nd` 在没有编译 ScaLAPACK 时为 1，否则取不超过每个 pool 进程数的最大平方数。注意 7.1 起，命令行没有给出这些选项时 pw.x 会自动猜测（见 2.3 节），所以"不写"不一定等于"取 1"，以输出开头的并行信息为准。

### 2.2 各层级加速了哪些例程

PWscf 用户指南的总结：

| 并行方式 | 加速 | 不加速 |
|---|---|---|
| FFT（平面波） | 几乎所有例程 | `cdiaghg`/`rdiaghg` |
| k 点 pools | `c_bands` 及其子例程（近线性）、部分 `sum_band` | `v_of_rho`、`newd`、`mix_rho` |
| 线性代数 | `cdiaghg`/`rdiaghg`（不总是） | — |
| task groups | `fftw` | — |
| OpenMP | `fftw` 及部分例程、部分线性代数 | — |

内存方面：只有 FFT（平面波）并行和线性代数并行会分摊内存，**其它层级都不分摊**。所以 `cannot allocate` 时的第一反应应该是"减少 pool 数、增加每个 pool 的进程数"，而不是增加 pool。

### 2.3 快速估计参数

需要知道三个数（都在 pw.x 输出的开头）：k 点数 N_k（`number of k points=`）、平滑 FFT 网格第三维 N₃（`FFT dimensions`）、Kohn-Sham 态数 M（`number of Kohn-Sham states=`）。PWscf 用户指南给出的界限：

- `-nk` ≤ N_k，且最好是 N_k 的约数；进程总数必须能被 `-nk` 整除（否则报 `invalid number of pools`）；
- 每个 pool 的进程数不超过 N₃，最好不超过 N₃/2–N₃/4，并且最好是 N₃ 的约数；
- M 不到几百时不必用线性代数并行（`-nd 1`）；
- 进程数超过 FFT 平面数时，用 `-nt 2`/`4`/`8` 等 task groups。

QE 7.1 起，如果命令行没有给出 `-nk`、`-nd`、`-nt`，pw.x 会自动猜测：进程数 N_p > N₃ 时，取最小的 N_k 使 N_p/N_k ≤ N₃/2；若仍超过 N₃，再用 task groups 使 N_p/N_k/N_t ≤ N₃/4；线性代数组取 n_d² ≤ N_p/N_k/N_t，且 M/n_d ≈ 100。自动猜测是一个合理的起点，但仍然建议对生产计算做一次小规模扫描（第 4 节）。

> [!TIP]
> 以下是常见体系的经验起点（只作参考，最终以扫描结果为准）：
> - **小原胞、很多 k 点**（体相金属、声子）：`-nk` 尽量大，例如每个 pool 2–8 个进程；
> - **大超胞、Γ 点或少量 k 点**（表面、MD、缺陷）：`-nk 1`，靠 FFT 并行；进程数超过 N₃ 时加 `-nt`；能带数上千时加 `-nd`；
> - **杂化泛函**：`-nb`（band groups）对 EXX 部分很有效（[第14章](14-hybrid-metagga.md)）；
> - **NEB、ph.x**：`-ni` 并行各个 image / q 点（[第11章](11-neb.md)、[第12章](12-phonon-dfpt.md)）；注意 ph.x 的电声计算不支持 `-ni`（[第16章](16-elph-epw.md)）。

### 2.4 MPI + OpenMP

编译时需要 `--enable-openmp`（configure）或相应的 CMake 选项。运行时：

```bash
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
srun -n 16 --cpus-per-task=2 --cpu-bind=cores pw.x -nk 4 -in pw.in > pw.out
```

QE 无法控制进程和线程绑定在哪个核上，这由操作系统和 MPI 启动器决定。官方用户指南的建议是：如果混合并行性能很差，很可能是多个线程挤在同一个核上，此时先设 `OMP_NUM_THREADS=1`。

### 2.5 GPU

- 编译：NVIDIA GPU 需要 NVHPC SDK。QE 7.6 的 CMake 使用新的 `QE_GPU`（如 `"openacc;cuda"`）和 `QE_GPU_ARCHS`（如 `sm_90`），旧的 `QE_ENABLE_CUDA` 等选项已删除，残留在 `CMakeCache.txt` 中也会报错；configure 用 `--with-cuda`、`--with-cuda-cc`、`--with-cuda-mpi`（[第1章](01-install.md)）；
- 运行：一般**每块 GPU 一个 MPI 进程**，优先用 pool 并行（`-nk` = GPU 数，在 k 点足够时）；
- 多 GPU 之间的通信依赖 CUDA-aware MPI；
- 输出开头会注明是否启用了 GPU 加速。7.5 起已经移除 PPCG 对角化，GPU 上常用默认的 Davidson。

### 2.6 内存与 I/O

`pw.x` 在初始化后打印内存估计（格式来自 `PW/src/memory_report.f90`；示意输出，数值仅供参考）：

```text
     Estimated max dynamical RAM per process >      xxx.xx MB

     Estimated total dynamical RAM >       xx.xx GB
```

设 `nstep=0`（对 relax/md）或创建 `prefix.EXIT` 后启动，程序会在初始化后停下，可以先看内存估计再正式提交。内存不足时（PWscf 用户指南的建议）：

- 使用更多进程，或在相同进程数下减少 pools（pools 不分摊内存）；
- 把 `nbnd` 降到必要的最小值；
- 设 `diago_david_ndim=2` 缩小 Davidson 工作空间，或改用 `diagonalization='cg'`（慢但几乎不需要工作空间）；
- 电荷密度占内存较多时，把 `mixing_ndim` 从 8 降到 4；
- 增加 `-nd`（分摊 M×M 子空间矩阵）。

`disk_io` 控制二进制文件和 XML 的写盘量：

| `disk_io` | 行为 | 用途 |
|---|---|---|
| `'high'` | 每步 SCF 写电荷，波函数以"distributed"格式留在磁盘，同时保存 mixing 数据 | 基本不需要；续算并不依赖它 |
| `'medium'` | 每步 SCF 写电荷；每个进程上 k 点多于一个时把波函数留在磁盘，结束时写盘（非 scf 计算的默认值） | k 点很多、内存紧张时 |
| `'low'` | 每步 SCF 写电荷；波函数留在内存，结束时写盘（scf 的默认值） | 一般情况 |
| `'nowf'` | 只保存 XML 和电荷密度，不写波函数 | 不需要后续 nscf/后处理时；**不能续算** |
| `'minimal'` / `'none'` | 只写 XML / 什么都不写 | 大批量单点、MLFF 标注 |

`outdir` 最好放在并行文件系统上，避免 NFS。也可以用环境变量 `ESPRESSO_TMPDIR` 设置默认的 `outdir`，用 `ESPRESSO_PSEUDO` 设置默认的 `pseudo_dir`（输入中没有写时生效）。

---

## 3. 读时间报告，诊断并行瓶颈

`pw.x` 结束时会打印各例程的耗时，最后一行类似（格式来自 `UtilXlib/clocks_handler.f90`；示意输出，数值仅供参考）：

```text
     PWSCF        :   xxm xx.xxs CPU   xxm xx.xxs WALL
```

正常的 NC 赝势计算中，大部分时间花在 `electrons` → `c_bands` → `cegterg`（或 Γ 点的 `regterg`）→ `h_psi`，FFT（`fftw`）是主要开销；`cdiaghg`/`rdiaghg` 只应占一小部分。US/PAW 赝势下 `sum_band`、`newd` 的比例会大一些，`tqr=.true.` 可以显著降低这部分开销。

PWscf 用户指南列出的"并行不当"症状和对策：

| 症状 | 对策 |
|---|---|
| `v_of_rho`、`newd`、`mix_rho` 占比很大；增加 pools 不再加速 | 同时使用 FFT 并行（每个 pool 多几个进程） |
| `cdiaghg`/`rdiaghg` 占比异常大 | 使用 k 点并行；使用线性代数并行（`-nd`，最好编译 ScaLAPACK/ELPA） |
| `fft_scatter` 占比大，或其 CPU 与 WALL 时间相差很大 | 网络不够快时，FFT 并行不要超过 4–8 个进程；改用 k 点并行 |
| 增加 FFT 并行的进程数不再加速 | 使用 task groups：`-nt 4` 或 `-nt 8` |
| CPU 时间远小于 WALL 时间（未开 OpenMP 时） | I/O 或通信瓶颈；检查 scratch 目录和 `disk_io` |

```bash
grep -A1 -e "     electrons    :" -e "     c_bands      :" -e "     cdiaghg      :" -e "     fft_scatter  :" pw.out
grep "PWSCF *:.*WALL" pw.out
```

---

## 4. 实战一：对 Si 做并行参数扫描

示例输入 [si.scf.in](../examples/18-automation-parallel/si.scf.in) 是一个 2 原子的 Si SCF（8×8×8 k 网格）。脚本 [scan_parallel.sh](../examples/18-automation-parallel/scan_parallel.sh) 在固定总核数（`NTOT`）下，遍历 `-nk` ∈ {1, 2, 4, 8} 和 `OMP_NUM_THREADS` ∈ {1, 2}，为每个组合生成一个目录和作业脚本：

```bash
cd QE/examples/18-automation-parallel
bash scan_parallel.sh si.scf.in                     # prepare only
SUBMIT=1 NTOT=16 bash scan_parallel.sh si.scf.in    # submit to Slurm
bash scan_parallel.sh collect                       # summarise wall time and RAM
grep '^!' scan_*/pw.out                             # energies must agree
```

汇总格式（示意输出，数值仅供参考）：

```text
run                      PWSCF_WALL  max_RAM_estimate
scan_nk1_omp1                 xx.xxs        xx.xx MB
scan_nk2_omp1                 xx.xxs        xx.xx MB
scan_nk4_omp1                 xx.xxs        xx.xx MB
scan_nk8_omp1                 xx.xxs        xx.xx MB
scan_nk1_omp2                 xx.xxs        xx.xx MB
...
```

怎么读：

- 不同并行设置的总能量（`!` 行）应该在 `conv_thr` 的量级内一致；差别大说明编译或库有问题；
- 对这种小体系，`-nk` 越大通常越快，直到每个 pool 只剩 1–2 个进程；
- 每进程的内存估计随 `-nk` 增大而增大，这符合"pools 不分摊内存"；
- 把扫描得到的最优设置写进生产作业模板。对大体系，同样的思路要换成扫描 `-nt`、`-nd` 和每个 pool 的进程数。

> [!TIP]
> 扫描用的输入应与生产计算相同的体系和截断能，但可以把 `electron_maxstep` 限制为几步，或只比较前几个 SCF 迭代的平均耗时，以节省机时。

---

## 5. 实战二：作业模板与正常停止、续算

### 5.1 Slurm 模板

[job_qe_slurm.sh](../examples/18-automation-parallel/job_qe_slurm.sh) 中值得注意的几行：

```bash
#SBATCH --ntasks-per-node=16       # MPI ranks per node
#SBATCH --cpus-per-task=2          # OpenMP threads per rank (1 = pure MPI)
#SBATCH --signal=B:USR1@300        # send USR1 to this script 5 min before the limit

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
rm -f "tmp/${PREFIX}.EXIT" "${PREFIX}.EXIT"   # a stale EXIT file stops the run immediately
trap 'echo "USR1: requesting clean stop"; touch "tmp/${PREFIX}.EXIT"' USR1

srun -n "${NP}" --cpu-bind=cores pw.x -nk "${NK}" -in "${INPUT}" > "${OUTPUT}" &
```

- `pw.x` 用 `-in`（或 `-i`、`-inp`、`-input`）读输入比 `< input` 重定向更可靠，官方用户指南提到某些 MPI 实现在并行时处理标准输入有问题；
- 后台运行 + `wait` 是为了让 bash 能在程序运行期间响应 `USR1` 信号；
- 仓库通用模板 `QE/scripts/job_qe_slurm.sh`、`QE/scripts/job_qe_pbs.sh` 提供了相同的功能，用 `bash ../scripts/job_qe_slurm.sh -h` 或阅读文件头注释查看用法。

### 5.2 正常停止与续算

PWscf 用户指南（"Restarting"）规定，QE 5.1 起**只能从正常停止的状态续算**。正常停止的三种方式：

1. 到达 `max_seconds`（`&CONTROL`，默认 1e7 秒，相当于不限时）；
2. 在工作目录或 `outdir` 中创建 `prefix.EXIT`；
3. （实验性）编译时加 `-D__TERMINATE_GRACEFULLY`，程序会拦截 SIGINT、SIGTERM、SIGUSR1、SIGUSR2、SIGXCPU 等信号。

续算：

```fortran
&CONTROL
  restart_mode = 'restart'     ! continue an interrupted run; same number of processes and parallelization
  max_seconds  = 82800
  ...
/
```

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `restart_mode` | `'from_scratch'` | `'restart'` 只能用来接着跑被正常停止的同一个计算，**必须使用相同的进程数和并行方式**；不能用来开始新计算或做 nscf；此时 `startingwfc`、`startingpot` 默认为 `'file'` |
| `max_seconds` | 1e7 | 设为作业时限减去写盘所需时间（大体系预留 10–30 分钟） |
| `disk_io` | scf 为 `'low'` | `'nowf'` 时无法续算 |

> [!WARNING]
> 不要把 `restart_mode='restart'` 当作"从上一个计算的波函数/电荷开始"的开关。要用已有电荷密度开始一个新计算，应设置 `startingpot='file'`（和/或 `startingwfc='file'`），并保持 `prefix`、`outdir` 一致。

relax/vc-relax/md 的续算也遵循同样的规则。另一种更稳妥的做法是：从输出中取出最后一个结构，作为新的 `from_scratch` 计算的起点（ASE 可以直接读出，见下一节）。

---

## 6. 实战三：Python 自动化

### 6.1 ASE：`EspressoProfile`（ASE ≥ 3.23）

ASE 3.23 起，Espresso 计算器改用 profile 写法；旧的 `Espresso(command=...)` 会直接报错（错误信息提示改用 `Espresso(profile=EspressoProfile(...))`）。示例 [ase_si_scf.py](../examples/18-automation-parallel/ase_si_scf.py) 的核心部分：

```python
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile

atoms = bulk("Si", "diamond", a=5.43)
profile = EspressoProfile(command="mpirun -np 4 pw.x -nk 2", pseudo_dir="/path/to/pseudo")
input_data = {
    "control": {"calculation": "scf", "prefix": "si_ase", "outdir": "./tmp",
                "tprnfor": True, "tstress": True},
    "system": {"ecutwfc": 30.0, "ecutrho": 240.0, "occupations": "fixed"},
    "electrons": {"conv_thr": 1.0e-10, "mixing_beta": 0.7},
}
atoms.calc = Espresso(profile=profile, directory="ase_si", input_data=input_data,
                      pseudopotentials={"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"},
                      kpts=(8, 8, 8))
print(atoms.get_potential_energy())   # eV
```

要点：

- 输入写到 `directory/espresso.pwi`，输出读自 `espresso.pwo`；`EspressoProfile` 的 `pseudo_dir` 会作为 `&CONTROL` 中 `pseudo_dir` 的默认值；
- `input_data` 中的数值**不做单位换算**（Ry、Bohr 等 QE 原生单位）；
- ASE 写出的结构总是 `ibrav=0` + `CELL_PARAMETERS`，k 点支持 Monkhorst-Pack 网格（`kpts`、`koffset`）或 `kspacing`；
- 读输出：`ase.io.read("pw.out", index=":", format="espresso-out")` 可以读出 relax/MD 的每一步（能量 eV、力 eV/Å）；
- 示例脚本默认只写输入（不加 `--run`），方便先检查。

### 6.2 pymatgen、qe-tools、qeschema、postqe

| 工具 | 功能 | 适合 |
|---|---|---|
| [pymatgen](https://github.com/materialsproject/pymatgen) `pymatgen.io.pwscf` | `PWInput`（由 `Structure` 写输入）、`PWOutput`（正则解析输出） | 已经用 pymatgen 做结构处理和高通量 |
| [pymatgen-io-espresso](https://github.com/Griffin-Group/pymatgen-io-espresso) | 把 QE 的 XML 转成类似 `Vasprun` 的对象 | 习惯 VASP 后处理流程的用户 |
| [qe-tools](https://github.com/aiidateam/qe-tools) | 纯 Python 解析 pw.x/cp.x 输入输出 | 轻量脚本 |
| [qeschema](https://github.com/QEF/qeschema) | 官方 XML Schema 与读取 `data-file-schema.xml` 的库 | 需要严格按 schema 读取 |
| [postqe](https://github.com/QEF/postqe) | 官方 Python 后处理（EOS、能带、DOS、电荷） | 不想写 Fortran 后处理输入 |

不想引入依赖时，标准库就够了。`pw.x` 的 XML 结构（`qes.xsd`）中常用的元素：

| XPath（相对根元素） | 内容 | 单位 |
|---|---|---|
| `output/convergence_info/scf_conv/convergence_achieved` | SCF 是否收敛 | — |
| `output/total_energy/etot` | 总能 | Ha |
| `output/band_structure/fermi_energy` | 费米能（金属） | Ha |
| `output/band_structure/highestOccupiedLevel`、`lowestUnoccupiedLevel` | HOMO / LUMO（绝缘体） | Ha |
| `output/atomic_structure/atomic_positions/atom` | 原子（`name` 属性 + 坐标） | Bohr |
| `output/atomic_structure/cell/a1`…`a3` | 晶格矢量 | Bohr |
| `output/forces` | 力（3×nat） | Ha/Bohr |
| `output/stress` | 应力（3×3） | Ha/Bohr³ |

示例 [parse_xml.py](../examples/18-automation-parallel/parse_xml.py) 用 `xml.etree.ElementTree` 读出这些量并换算成 eV、Å、kbar：

```bash
python parse_xml.py tmp/si.save/data-file-schema.xml --json si_result.json
```

注意 XML 中的能量单位是 **Hartree**，而 `pw.x` 文本输出用 Ry。本仓库的 `qe_parse_out.py` 解析文本输出（总能、费米能、SCF 历史、力、应力、磁矩、最终结构），在 `QE/tutorials/` 目录下用 `python ../scripts/qe_parse_out.py -h` 查看用法；写输入前可以用 `qe_check_input.py` 做静态检查（检查 namelist 结尾的 `/`、CRLF、单位缺失、废弃关键字等），用 `python ../scripts/qe_check_input.py -h` 查看用法。

### 6.3 工作流框架

| 工具 | 说明 |
|---|---|
| [aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso) | AiiDA 插件；`PwBaseWorkChain` 能自动处理常见错误（SCF 不收敛、超时、对角化失败等）并重启；另有弛豫、能带、声子工作流；数据库自动记录完整溯源 |
| [aiidalab-qe](https://github.com/aiidalab/aiidalab-qe) | 浏览器图形界面：结构 → 弛豫 → 能带/PDOS，适合教学和非专业用户 |
| [Materials Cloud QE 输入生成器](https://qeinputgenerator.materialscloud.io/) | 上传结构，自动给出基于 SSSP 的截断能、k 点和赝势的 pw.x 输入 |
| [aiida-pseudo](https://github.com/aiidateam/aiida-pseudo) | 一键安装 SSSP/PseudoDojo 赝势族并给出推荐截断 |
| 本仓库 `qe_inputgen.py` | 从 CIF/POSCAR 生成 pw.x 输入（显式单位、SSSP 截断、k 间距），在 `QE/tutorials/` 目录下用 `python ../scripts/qe_inputgen.py -h` 查看用法 |

经验：几十个计算以内，bash 循环 + 模板替换（[第4章](04-convergence.md) 的 `qe_conv_test.sh`）或 ASE 脚本就足够；上百个计算、需要自动纠错和数据溯源时，再考虑 AiiDA。

---

## 7. 报错速查表

下面汇总全系列中最常见的问题。QE 的错误统一以如下格式输出：

```text
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine <routine> (<code>):
     <message>
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

先看 routine 名和 message，再到下表中查找。并行运行时错误信息有时只出现在某个进程的 `CRASH` 文件里，作业日志中看不到：

```bash
grep -A3 "Error in routine" pw.out CRASH 2>/dev/null
tail -30 pw.out
```

### 7.1 输入与文件

| 报错原文 | 原因 | 修复 |
|---|---|---|
| `Error in routine read_namelists: bad line in namelist &system: "..." (error could be in the previous line)` | 关键字拼错、放错 namelist、字符串没加引号、上一行少了逗号/引号 | 对照 INPUT_PW 检查报出的行**和它的上一行**；用 `qe_check_input.py` 做静态检查 |
| `could not find namelist &...` | namelist 顺序错、名字拼错，或上一个 namelist 漏了 `/` | 按 `&CONTROL`、`&SYSTEM`、`&ELECTRONS`、`&IONS`、`&CELL` 的顺序写，每个以 `/` 结束 |
| 读输入时出现莫名其妙的错误 | Windows 的 CRLF 换行、制表符、非 ASCII 引号（如中文引号） | `file pw.in` 检查；`dos2unix pw.in`；用纯文本编辑器 |
| `Error in routine readpp: file ./pseudo/XX.UPF not found`（文件名部分为实际路径） | `pseudo_dir` 路径错、文件名大小写不符 | 检查 `pseudo_dir`（默认读 `$ESPRESSO_PSEUDO`，再退回 `$HOME/espresso/pseudo/`）和 `ATOMIC_SPECIES` 中的文件名 |
| `Error in routine readpp: inconsistent DFT read from PP files` | 不同元素的赝势泛函不一致（如 LDA 与 PBE 混用） | 统一使用同一泛函的赝势族；`input_dft` 强制覆盖要慎用 |
| `Error in routine card_atomic_positions: wrong number of columns ...`、`nat out of range` 等 `card_atomic_species`/`card_atomic_positions` 报错 | `ntyp`/`nat` 与卡片行数不一致，行中列数不对，或物种名不匹配 | 核对 `ntyp`、`nat` 与卡片内容 |
| `Error in routine cell_base_init: ibrav=0: must read cell parameters` | `ibrav=0` 但没有 `CELL_PARAMETERS` | 补上 `CELL_PARAMETERS` 并写明单位（angstrom、bohr 或 alat）；`celldm(1)` 与 `A` 不能同时给 |
| 坐标或晶胞单位缺失的 DEPRECATED 警告 | `ATOMIC_POSITIONS`/`CELL_PARAMETERS` 没写单位 | 显式写 `angstrom`、`crystal`、`bohr` |
| 旧 DFT+U 输入报错（`lda_plus_u`、`Hubbard_U` 等） | 7.1 起必须使用 `HUBBARD` 卡，旧关键字已删除 | 改写成 `HUBBARD (ortho-atomic)` + `U Fe-3d 4.0`（[第8章](08-dftu-hp.md)） |
| `wf_collect`、`lkpoint_dir` 被忽略 | 已为 OBSOLETE | 删除 |
| nscf/bands/后处理找不到数据 | `prefix` 或 `outdir` 与 scf 不一致；scf 用了 `disk_io='nowf'` | 统一 `prefix`/`outdir`；需要后处理时不要用 `'nowf'` |

### 7.2 电子结构与收敛

| 报错原文 | 原因 | 修复 |
|---|---|---|
| `convergence NOT achieved after 100 iterations: stopping` | 金属没加 smearing；`mixing_beta` 太大；磁性初值不当；表面/长胞的电荷振荡 | 降低 `mixing_beta`（0.7 → 0.3/0.1）；增大 `electron_maxstep`；金属用 `occupations='smearing'`；检查 `starting_magnetization`；表面用 `mixing_mode='local-TF'`；适当增加 `mixing_ndim`；US 赝势可提高 `ecutrho` |
| `Message from routine setup: the system is metallic, specify occupations`（提示信息，不终止运行） | 电子数为奇数而 `occupations='fixed'`（默认） | 用 `occupations='smearing'`（DOS 用 `'tetrahedra_opt'`）；自旋极化且未给 smearing 时则直接报错 `spin-polarized system, specify occupations` |
| `Error in routine cdiaghg: S matrix not positive definite`（或 `rdiaghg`；也可能是 `problems computing cholesky`、`eigenvectors failed to converge`） | 原子距离过近或单位写错；赝势有 ghost 态；US 赝势电荷出现负值；数学库问题 | 先检查结构（单位、坐标）；改 `diagonalization='cg'`（或 `'paro'`）；提高 `ecutrho`；换赝势；换用非优化的 BLAS/LAPACK 排查。注意 7.5 起没有 `'ppcg'` |
| `Error in routine c_bands: too many bands are not converged` | 对角化失败，常与结构问题或 `nbnd` 过大有关 | `diagonalization='cg'`；检查结构；适当减小 `nbnd` |
| `c_bands:  N eigenvalues not converged` | SCF 早期常见 | 只出现在前几步时可忽略；最后几步仍出现则说明有问题 |
| `Error in routine efermig: internal error, cannot bracket Ef` | 电子数/能带数不对；smearing 值不合理；少量 k 点下用 m-p smearing | 增加 `nbnd`；改用 `'gaussian'` 或 `'mv'` smearing；检查 `tot_charge` |
| `Warning: negative or imaginary charge...` 等 | NLCC 或 US 赝势在有限截断下的正常现象（PWscf 用户指南 Troubleshooting） | 量级接近 0.1 时说明有严重问题：检查结构与赝势，提高 `ecutrho` |
| 对称性少于预期 / `fractional translation` 相关警告 | 坐标精度不够；分数平移与 FFT 网格不兼容；结构是超胞 | 用 spglib 对称化结构；给足有效数字；必要时设 `nosym`（k 点会增多） |
| 结果在不同机器上略有不同 | FFT 网格、运算顺序不同 | 差别很小时正常；需要更高一致性时收紧 `conv_thr` |

### 7.3 结构优化、MD

| 报错/现象 | 原因 | 修复 |
|---|---|---|
| relax/vc-relax 中途出现 cdiaghg 错误 | 离子步太大，原子过近 | 减小 `trust_radius_ini`；先固定晶胞做 `relax`（[第5章](05-relax-eos.md)） |
| BFGS 反复重置、迟迟不收敛 | 软模、平坦势能面；力有噪声 | 收紧 `conv_thr`；放宽 `forc_conv_thr`；换 `ion_dynamics='fire'` |
| vc-relax 结束后最后一次 scf 的能量/压强差别大 | Pulay 应力（截断能不够） | 提高 `ecutwfc`/`ecutrho`，再做一次 vc-relax |
| 续算没有接上：`Message from routine iosys: restart disabled: needed files not found`，或读重启文件时报错 | 上次没有正常停止、`prefix`/`outdir` 变了，或换了并行设置 | 按 5.2 节正常停止；保持 `prefix`/`outdir`/并行设置不变；否则从最后结构重新开始 |
| `Error in routine davcio: error writing file "..."`（或 `error reading file`） | `outdir` 不可写、磁盘满、多个作业共用同一 `prefix`+`outdir`、NFS 不稳定 | 检查权限和配额；每个作业用独立的 `outdir` |
| MD 能量漂移 | `conv_thr` 太松、`dt` 太大 | 见 [第17章](17-md-cp.md) |

### 7.4 并行与资源

| 报错原文 | 原因 | 修复 |
|---|---|---|
| `Error in routine mp_start_pools: invalid number of pools, parent_nproc /= nproc_pool * npool` | 进程数不能被 `-nk` 整除 | 让 `-nk` 整除总进程数 |
| `invalid number of pools, out of range` | `-nk` 大于进程数或小于 1 | 调整 `-nk` |
| `Error in routine setup: Gamma-only calculations not allowed with pools` | `K_POINTS gamma` + `-nk > 1` | 去掉 `-nk` |
| `Berry Phase/electric fields not implemented with pools` | `lberry`/`lelfield` + pools | 去掉 `-nk` |
| `stress for hybrid functionals not available with pools` | 杂化泛函 + `tstress` + pools | 去掉 `-nk` 或不算应力；另注意 7.6.1（开发版预告）起 USPP/PAW 杂化泛函不再计算应力 |
| `pools with no k-points not allowed for hybrid functionals` | pool 数多于 k 点数 | 减少 `-nk` |
| `cannot allocate` / 被系统 OOM 杀掉 / `segmentation fault` | 内存不足或栈太小 | 增加进程数、减少 pools、加 `-nd`、减小 `nbnd`/`mixing_ndim`/`diago_david_ndim`；`ulimit -s unlimited` |
| 程序没有报错就消失，或只在某个进程出错后挂起 | 并行环境问题；输出丢失 | 交互式小规模运行重现；查看 `CRASH` 文件；检查是否混用了不同的 MPI |
| 并行运行比串行还慢 | 混合了不同 MPI 库，或在串行编译的程序上用 `mpirun` | 输出开头确认 `Parallel version (MPI)` 与进程数；统一编译器和 MPI 的 module |
| 7.5 + Open MPI 5.0.x 的部分测试挂起（社区报告） | MPI 版本兼容问题 | 改用 Open MPI 4.1.x（[第1章](01-install.md)） |
| CMake 报 `QE_ENABLE_CUDA` 等已移除的选项（7.6） | 旧的 GPU 编译选项 | 删除 `CMakeCache.txt`，改用 `QE_GPU`、`QE_GPU_ARCHS` |

### 7.5 后处理与专题程序

| 报错/现象 | 原因 | 修复 | 章节 |
|---|---|---|---|
| bands/DOS 缺少导带 | `nbnd` 太小 | 显式设置 `nbnd` | [第6章](06-band-dos-pdos.md) |
| bands 计算给出的费米能不对 | `calculation='bands'` 不重新确定费米能 | 以 scf/nscf 的费米能为准 | [第6章](06-band-dos-pdos.md) |
| DOS 呈锯齿状 | nscf 的 k 点太少 | 加密 k 网格，用 `occupations='tetrahedra_opt'` | [第6章](06-band-dos-pdos.md) |
| 磁性体系收敛到非磁解 | 没设 `nspin=2` 或 `starting_magnetization` | 设置初始磁矩；AFM 要把同种元素拆成两个 species | [第7章](07-magnetism-soc.md) |
| hp.x 与 pw.x 的 U 不一致 | Hubbard 投影类型不同 | `atomic` 与 `ortho-atomic` 得到的 U 不可互换 | [第8章](08-dftu-hp.md) |
| 声子虚频 | 结构未充分弛豫；`conv_thr`、`tr2_ph` 太松；没加 ASR | scf `conv_thr` ≤ 1e-12，`tr2_ph` 1e-14–1e-16，`asr='crystal'` | [第12章](12-phonon-dfpt.md) |
| `The phonon code with US-PP and raman or elop not yet available`（PAW 时为 `...with paw and raman or elop...`） | `lraman`/`elop` 只支持模守恒赝势；GPU 版本也不支持 Raman | 换模守恒赝势（泛函方面的限制见第13章）；用 CPU 版本 | [第13章](13-dielectric-ir-raman.md) |
| 金属 ph.x 报错或结果异常 | scf 没用 smearing | 金属 scf 用 `occupations='smearing'` | [第12章](12-phonon-dfpt.md) |
| 版本之间结果略有差别 | 7.5 修改了 Simpson 径向积分；7.6 修改了 `epsilon.x` | 比较时注明版本 | [第1章](01-install.md) |

> [!TIP]
> 向邮件列表（users@lists.quantum-espresso.org）或论坛提问前，准备好：QE 版本和编译方式、完整的输入、报错前后的输出、并行命令行，以及你已经试过的办法。官方用户指南也强调，"it crashes with ...(obscure MPI error)" 这类信息几乎无法得到回答。

---

## 参数选择与收敛建议

| 场景 | 建议 |
|---|---|
| 首次在新集群运行 | 先跑官方 `test-suite` 或 `PW/examples/example01`，再跑本章的并行扫描 |
| k 点多的小体系 | `-nk` = k 点数的约数，每个 pool 2–8 个进程；`-nd 1` |
| 大超胞（数百原子以上） | `-nk` 小或为 1；每个 pool 的进程数 ≤ N₃/2；能带上千时 `-nd` 取平方数使 M/n_d ≈ 100；进程数 > N₃ 时加 `-nt` |
| 内存紧张 | 减少 pools、增加每个 pool 的进程数、加 `-nd`；`disk_io='medium'` |
| OpenMP | 每个 MPI 进程 2–4 个线程起步；绑定到核；表现异常时退回 `OMP_NUM_THREADS=1` |
| GPU | 每块 GPU 一个 MPI 进程；k 点够多时 `-nk` = GPU 数；使用 CUDA-aware MPI |
| 长作业 | `max_seconds` = 时限减去写盘时间；配合 `prefix.EXIT` 与信号；续算保持并行设置不变 |
| scratch | `outdir` 放在并行文件系统或节点本地盘（注意续算时能否访问）；避免 NFS |
| 高通量 | 统一的模板 + 脚本化检查（`JOB DONE`、`convergence has been achieved`）；需要纠错和溯源时用 AiiDA |

---

## 常见问题与报错

这一部分的完整内容见第 7 节的速查表。下面补充几条容易被忽略的检查习惯：

| 现象 | 原因 | 处理 |
|---|---|---|
| 输出里没有 `JOB DONE` | 程序异常退出或被调度系统杀掉 | 查看作业日志、`CRASH` 文件和输出末尾；用 `grep -c "JOB DONE" *.out` 批量检查 |
| 作业"成功"但结果不对 | SCF 没收敛也继续往下走（例如 relax 中途） | 检查每一步都有 `convergence has been achieved`；relax 检查 `bfgs converged` |
| 同时提交的多个作业互相覆盖文件 | 共用了同一 `prefix` + `outdir` | 每个作业使用独立的 `outdir` |
| 换了机器后速度慢很多 | 数学库/FFT 库不同，或绑定错误 | 输出开头核对库和并行信息；检查 `OMP_NUM_THREADS` 与 `--cpus-per-task` 一致 |
| 程序启动后立刻停止，没有报错 | 残留的 `prefix.EXIT`，或 `nstep=0` | 删除 EXIT 文件；检查 `nstep` |

---

## 练习

1. 用 `scan_parallel.sh` 在你的集群上扫描 Si 的 `-nk` 与 `OMP_NUM_THREADS`，画出墙钟时间对 `-nk` 的曲线，并说明内存估计为什么随 `-nk` 增大。
2. 用一个 100 原子以上的超胞（例如第9章的 slab，Γ 点或少量 k 点），比较 `-nd 1`、`-nd 4`、`-nd 16` 时 `cdiaghg` 的耗时。
3. 在一个 vc-relax 作业中设 `max_seconds` 为 5 分钟，确认程序正常停止后用 `restart_mode='restart'` 续算；再故意改变进程数续算，观察会发生什么。
4. 用 ASE 的 `EspressoProfile` 写一个脚本，对 Si 扫描 5 个晶格常数并拟合 EOS（可以配合 `qe_eos_fit.py`，在 `QE/tutorials/` 目录下用 `python ../scripts/qe_eos_fit.py -h` 查看用法）。
5. 用 `parse_xml.py` 和 `qe_parse_out.py` 分别读取同一个计算，核对总能（注意 Ha 与 Ry 的换算）。
6. 故意制造三种错误（漏掉 `/`、赝势文件名大小写错误、`-nk` 不整除进程数），记录真实报错信息，与第 7 节的表格对照。

---

## 参考资料

**官方文档**

- [QE 总用户手册（并行化、CMake、GPU）](https://www.quantum-espresso.org/Doc/user_guide/)：并行层级、`-ni/-nk/-nb/-nt/-nd`、并行 I/O、MPI-OpenMP 问题
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)（[PDF](https://www.quantum-espresso.org/Doc/user_guide_PDF/pw_user_guide.pdf)）：并行参数的快速估计、时间报告的解读、重启、Troubleshooting 一节
- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)：`restart_mode`、`max_seconds`、`disk_io`、`pseudo_dir`、`outdir`
- [QE GitLab Wiki（CMake/GPU 编译）](https://gitlab.com/QEF/q-e/-/wikis/home) / [release notes](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes)
- [users 邮件列表](https://lists.quantum-espresso.org/mailman/listinfo/users)
- [ASE Espresso 计算器](https://docs.ase-lib.org/ase/calculators/espresso.html)（ASE 文档 2026 年迁到 docs.ase-lib.org，旧的 wiki.fysik.dtu.dk 链接已失效）
- [aiida-quantumespresso 文档](https://aiida-quantumespresso.readthedocs.io)

**课程与教程（只链接）**

- [ENCCS / MaX: Efficient materials modelling on HPC](https://enccs.github.io/efficient-materials-modelling-on-hpc/)：`-nk`、`-nd`、`-nt`、OpenMP 和 GPU 的上机练习，强烈推荐配合本章阅读
- [MaX e-School 2021](https://indico.ictp.it/event/9616/)：QE 官方暑校，含并行与 GPU 的讲座
- [Advanced QE school 2023（Materials Cloud）](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023)：含 AiiDA 上机
- [Imperial MSE404-MM](https://imperialcollegelondon.github.io/MSE404-MM/)：Lab 1–2（Linux 与 QE 输入输出）适合补齐运行环境和读输出的基础
- Cornell PARADIM [Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials)：在集群（MARCC）上配置和运行 QE 的流程；编译和版本部分已过时，只作参考

**相关仓库（只链接）**

- [QEF/q-e](https://github.com/QEF/q-e)（官方 GitHub 镜像，开发在 GitLab）/ [QEF/qeschema](https://github.com/QEF/qeschema) / [QEF/postqe](https://github.com/QEF/postqe)
- [aiidateam/aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso) / [aiidalab/aiidalab-qe](https://github.com/aiidalab/aiidalab-qe) / [aiidateam/qe-tools](https://github.com/aiidateam/qe-tools)
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) / [Griffin-Group/pymatgen-io-espresso](https://github.com/Griffin-Group/pymatgen-io-espresso)
- [hfp/xconfigure](https://github.com/hfp/xconfigure)：Intel 平台编译 QE 的配方
- [sebranchett/MoS2-QE-Yambo-pnas](https://github.com/sebranchett/MoS2-QE-Yambo-pnas)：集群上 QE + Yambo 的作业脚本示例

**本仓库脚本**：`qe_inputgen.py`、`qe_check_input.py`、`qe_parse_out.py`、`qe_conv_test.sh`、`job_qe_slurm.sh`、`job_qe_pbs.sh` 等，用 `python ../scripts/<name> -h`（`.sh` 用 `bash ../scripts/<name> -h`）查看用法。

**关键文献**

- P. Giannozzi et al., *J. Phys.: Condens. Matter* **21**, 395502 (2009). DOI: [10.1088/0953-8984/21/39/395502](https://doi.org/10.1088/0953-8984/21/39/395502)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)
- P. Giannozzi et al., *J. Chem. Phys.* **152**, 154105 (2020)，面向百亿亿次计算的 QE（并行与 GPU）. DOI: [10.1063/5.0005082](https://doi.org/10.1063/5.0005082)
- I. Carnimeo et al., *J. Chem. Theory Comput.* **19**, 6992 (2023)，QE 的 GPU 实现. DOI: [10.1021/acs.jctc.3c00249](https://doi.org/10.1021/acs.jctc.3c00249)
- S. P. Huber et al., *Sci. Data* **7**, 300 (2020)，AiiDA 1.0. DOI: [10.1038/s41597-020-00638-4](https://doi.org/10.1038/s41597-020-00638-4)
- G. Prandini et al., *npj Comput. Mater.* **4**, 72 (2018)，SSSP. DOI: [10.1038/s41524-018-0127-2](https://doi.org/10.1038/s41524-018-0127-2)

---

[← 上一章](17-md-cp.md) | [返回目录](../README.md)
