# 第1章 QE 概述、版本与安装（Installing Quantum ESPRESSO: conda-forge, Source Build, GPU & WSL）

本章介绍 Quantum ESPRESSO（QE）的定位：平面波 + 赝势（plane waves & pseudopotentials）的开源 DFT 套件，核心程序 pw.x、ph.x、pp.x、neb.x、cp.x、hp.x、EPW。我们梳理 7.0 → 7.6 的版本线（7.6 引入 `QE_GPU`/`QE_GPU_ARCHS`，旧的 `QE_ENABLE_CUDA` 等选项被删除），给出四种安装路径的完整命令：conda-forge、发行版包与容器、`configure` + `make` 源码编译、CMake 编译（含 NVIDIA GPU），还有 Windows 上的 WSL2。最后用 test-suite、官方 `PW/examples` 和一个 Si 冒烟测试确认安装可用。

> [!NOTE]
> 本章所有脚本与输入在示例目录 [01-install](../examples/01-install/README.md)。本仓库**只整理教程，不提供任何计算结果**：文中的终端输出都标注为“示意输出”，数值仅供参考。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 QE 与 VASP、CP2K 在方法和许可上的区别，知道 QE 的主要可执行程序分别做什么；
- 读懂版本号与发布节奏，知道 7.1、7.5、7.6 有哪些**会让旧输入或旧编译脚本报错**的变化；
- 按自己的机器（笔记本、工作站、无 root 权限的集群、GPU 节点、Windows）选对安装方式并装好；
- 用 `pw.x -i x.in > x.out` 和 `mpirun -np N pw.x -nk K -i x.in` 运行，并用 test-suite 和冒烟测试判断安装是否真的可用。

| 你的情况 | 推荐安装方式 | 本章小节 |
|---|---|---|
| 只想尽快跑通第一个算例（Linux / macOS / WSL2） | conda-forge 的 `qe` 包（当前为 7.6） | §4 |
| 学校集群已经有 `module load quantum-espresso` | 直接用 module，先用 §9 的方法检查版本和并行库 | §9 |
| 集群生产计算（CPU），需要 ScaLAPACK/ELPA、Intel MKL 等优化库 | `configure` + `make` 或 CMake 源码编译 | §6、§7 |
| NVIDIA GPU 节点 | NVHPC SDK + CMake（`QE_GPU`）或 `configure --with-cuda` | §8 |
| Windows 笔记本 | WSL2 + Ubuntu，再用 conda 或源码编译 | §10 |
| 教学演示、不想装任何东西 | 浏览器版 nanoHUB、虚拟机 Quantum Mobile | §5 |

> [!TIP]
> Linux 基础、conda、SSH 和文件传输见 [公共篇第1章](../../Common/tutorials/01-linux-hpc.md)，Slurm/PBS 作业系统见 [公共篇第2章](../../Common/tutorials/02-job-scheduler.md)。

---

## 2. 原理速览：QE 是什么

### 2.1 一句话定位

Quantum ESPRESSO（opEn-Source Package for Research in Electronic Structure, Simulation, and Optimization）是一个以 GPL-2.0 许可发布的程序套件。Kohn–Sham 轨道用**平面波**展开，离子实用**赝势**（模守恒 NC、超软 US、PAW）表示。波函数的展开受截断能 `ecutwfc` 控制：

$$
\psi_{n\mathbf k}(\mathbf r)=\sum_{|\mathbf k+\mathbf G|^2/2\,<\,E_\mathrm{cut}} c_{n,\mathbf k+\mathbf G}\; e^{i(\mathbf k+\mathbf G)\cdot\mathbf r}
$$

QE 内部使用 **Rydberg 原子单位**：能量单位是 Ry（1 Ry = 13.6057 eV），长度单位是 Bohr（1 Bohr = 0.529177 Å）。这和 VASP（eV、Å）不同，是新用户最常见的误区之一，[第2章](02-input-anatomy.md) 会详细讲。

### 2.2 主要程序一览

| 程序 | 所在包 | 作用 | 本系列章节 |
|---|---|---|---|
| `pw.x` | PWscf | SCF、nscf、能带、结构优化、变胞优化、BOMD | 第2–9、14、17章 |
| `ph.x`、`q2r.x`、`matdyn.x`、`dynmat.x` | PHonon | DFPT 声子、介电张量、Born 有效电荷、拉曼 | [第12章](12-phonon-dfpt.md)、[第13章](13-dielectric-ir-raman.md) |
| `pp.x`、`dos.x`、`projwfc.x`、`bands.x`、`average.x` | PP | 电荷密度、势、DOS/PDOS、能带整理、平面平均 | [第6章](06-band-dos-pdos.md)、[第10章](10-charge-elf-bader.md) |
| `neb.x` | NEB | 微动弹性带过渡态搜索 | [第11章](11-neb.md) |
| `hp.x` | HP | 线性响应计算 Hubbard U/V | [第8章](08-dftu-hp.md) |
| `cp.x` | CPV | Car–Parrinello 分子动力学 | [第17章](17-md-cp.md) |
| `epw.x` | EPW | 基于 Wannier 插值的电声耦合、迁移率、超导 | [第16章](16-elph-epw.md) |
| `pw2wannier90.x` | PP | 与 Wannier90 的接口 | [第15章](15-wannier90.md) |
| `ld1.x` | atomic | 原子计算与赝势生成（pslibrary 用它） | [第3章](03-pseudopotentials.md) |
| `ev.x` | PW/tools | 状态方程拟合 | [第5章](05-relax-eos.md) |

### 2.3 与 VASP、CP2K 的对比

| 方面 | QE | VASP | CP2K |
|---|---|---|---|
| 许可 | GPL-2.0，免费、开源 | 商业许可 | GPL，开源 |
| 基组 | 平面波 | 平面波 | 高斯基 + 平面波（GPW） |
| 赝势 | NC / US / PAW，UPF 格式，可自由下载（SSSP、PseudoDojo、pslibrary） | PAW（POTCAR，随许可证分发） | GTH 赝势 |
| 强项 | DFPT 声子与电声（ph.x、EPW）、线性响应 U（hp.x）、Wannier 接口、源码可改 | 功能全面、文档统一、工业界常用 | 大体系、AIMD、QM/MM |
| 输入 | 单个文本文件（namelist + card） | INCAR/POSCAR/KPOINTS/POTCAR | 嵌套 section 文本 |
| 单位 | Ry、Bohr（坐标单位可选 Å） | eV、Å | Hartree、Å |

从 VASP 转过来的读者可以对照 [VASP 第1章](../../VASP/tutorials/01-setup-compile.md)，从 CP2K 转过来的读者可以对照 [CP2K 第1章](../../CP2K/tutorials/01-install.md)。

---

## 3. 版本线与关键变化

### 3.1 发布历史

QE 的主仓库在 **GitLab**（<https://gitlab.com/QEF/q-e>），GitHub 上的 `QEF/q-e` 只是镜像。issue 和 merge request 都要在 GitLab 上提。

| 版本 | 发布时间 | 与本教程相关的关键变化 |
|---|---|---|
| 7.0 | 2021-12 | GPU 支持大幅扩展；新增 RMM-DIIS 对角化（`diagonalization='rmm-davidson'`/`'rmm-paro'`） |
| 7.1 | 2022-06 | **DFT+U 改用 `HUBBARD` 卡**，删除 `lda_plus_u`、`Hubbard_U(i)` 等旧关键词（旧输入直接报错）；不给 `-nk`/`-nd` 时自动选择部分并行参数 |
| 7.2 | 2023-03 | ph.x 支持 GPU；开始通过 OpenMP 支持 AMD/Intel GPU |
| 7.3 / 7.3.1 | 2023-12 / 2024-02 | 非共线 DFT+U/U+V（含力和应力） |
| 7.4 / 7.4.1 | 2024-10 / 2025-02 | 双化学势声子；`fixed_magnetization` 的 lambda 含义修改 |
| 7.5 | 2025-08 | 轨道分辨 DFT+U；**删除 PPCG 对角化**；CUDA Fortran 大部分换成 OpenACC；pw.x 新增 Nose–Hoover NVT/NPT；偶数网格 Simpson 径向积分修改，数值可能有微小差异 |
| **7.6** | **2026-07-09** | **当前最新正式版**。EPW v6.1；NC meta-GGA 支持 NLCC；vdW-DF3-mc；epsilon.x 支持 USPP/PAW（结果与旧版不同）；autoconf 支持源码目录外构建；**CMake 删除 `QE_ENABLE_CUDA/ROCM/ONEAPI/OFFLOAD/OPENACC`，改用 `QE_GPU` + `QE_GPU_ARCHS`**；删除 `__SPIN_BALANCE` 预处理开关 |

> [!NOTE]
> 截至 2026-11-20，develop 分支的 release notes 中已经写入了 7.6.1 和 8.0 的条目，但两者**都还没有打 tag**。本系列以 7.6 为准，开发版的内容只作为“预告”提及（例如 7.6.1 起 USPP/PAW 下不再计算杂化泛函的应力，见 [第14章](14-hybrid-metagga.md)）。

### 3.2 已废弃、不要再写的写法

| 写法 | 状态 | 处理 |
|---|---|---|
| `wf_collect`、`lkpoint_dir` | INPUT_PW 中标注为 OBSOLETE（不再实现） | 删掉；新版总是以可移植格式写波函数 |
| `ATOMIC_POSITIONS`、`CELL_PARAMETERS` 后面不写单位 | DEPRECATED，运行时会打印 `DEPRECATED: no units specified in ATOMIC_POSITIONS card` | 显式写 `angstrom`、`crystal`、`bohr` 或 `alat` |
| `lda_plus_u`、`Hubbard_U(i)`、`U_projection_type` | 7.1 起删除 | 改用 `HUBBARD` 卡（[第8章](08-dftu-hp.md)） |
| `diagonalization='ppcg'` | 7.5 起不再支持 | 改用 `'david'`（默认）、`'cg'`、`'paro'` |
| 旧 XML/iotk 输出 | 6.x 起只写新的 `data-file-schema.xml` | 用 qeschema 等工具读新 XML（[第18章](18-automation-parallel-troubleshooting.md)） |

> [!IMPORTANT]
> 论文方法部分请写明**确切版本号**（如 “Quantum ESPRESSO 7.6”）和每个元素的赝势文件名。7.5（Simpson 积分）、7.6（epsilon.x）都有“结果与旧版略有不同”的说明，跨版本比较数据时必须注明版本。

---

## 4. 方式一：conda-forge（最快上手）

conda-forge 上的包名是 `qe`。截至 2026-11，最新版本是 7.6，提供 `linux-64`、`linux-aarch64`、`osx-64`、`osx-arm64` 四个平台，并分别有基于 Open MPI 和 MPICH 的构建（build 字符串分别以 `mpi_openmpi_` 和 `mpi_mpich_` 开头），另外还有 build 字符串里不带 `mpi_` 的非 MPI 构建。想要并行版本，就用下面注释里的写法显式指定 MPI 构建，安装后再用 `conda list qe` 看 build 字符串确认。

### 4.1 安装

```bash
# 1) create an isolated environment (mamba/micromamba accept the same syntax)
conda create -n qe -c conda-forge qe=7.6
conda activate qe

# 2) optional: force a specific MPI flavour
# conda create -n qe -c conda-forge "qe=7.6=mpi_openmpi_*"
# conda create -n qe -c conda-forge "qe=7.6=mpi_mpich_*"

# 3) check which executables are on PATH
which pw.x ph.x pp.x projwfc.x dos.x bands.x neb.x hp.x
conda list | grep -E "^(qe|openmpi|mpich|fftw|libblas|scalapack)"
```

### 4.2 运行

```bash
conda activate qe
export OMP_NUM_THREADS=1
mpirun -np 4 pw.x -i si_test.in > si_test.out
```

> [!WARNING]
> conda 环境里的 `mpirun` 必须是**同一个环境**里的那个 MPI。如果系统里还装了另一套 MPI（例如集群的 Intel MPI module），`which mpirun` 可能指向错误的那一个，表现是每个进程都各自独立运行一遍：输出头部 `Program PWSCF` 出现 N 次，进程数显示为 1。处理方法：`conda activate qe` 之后不要再 `module load` 其它 MPI，并用 `which mpirun` 确认路径在 conda 环境内。

### 4.3 适用范围

| 优点 | 局限 |
|---|---|
| 几分钟可用，无需 root，Linux/macOS 都能用 | 使用通用 BLAS/FFTW，性能不如针对本机优化的编译版本 |
| 版本跟进快（7.6 已上架） | 不含 GPU 加速 |
| 适合学习、调试输入、小体系计算 | 跨节点多机并行要和集群的网络/MPI 配合，通常不如集群自己编译的版本 |

---

## 5. 方式二：发行版包、容器与虚拟机

### 5.1 Ubuntu/Debian 的 apt 包（不推荐用于新项目）

```bash
sudo apt install quantum-espresso
```

截至 2026-11，Ubuntu 22.04、24.04、26.04 仓库里的 `quantum-espresso` 都是 **6.7**。它不支持 7.1 起的 `HUBBARD` 卡、7.5 起的 Nose–Hoover MD 等新功能，与本系列的输入不完全兼容。只建议在“系统里临时需要一个 pw.x”时使用。

### 5.2 NVIDIA NGC 容器（GPU）

NVIDIA NGC 上有官方 GPU 容器 `nvcr.io/hpc/quantum_espresso`。截至 2026-11，最新标签仍是 `qe-7.3.1`，比 7.6 落后好几个版本。它适合快速测试 GPU 性能，但输入语法以 7.3.1 为准。

```bash
# Apptainer / Singularity on a GPU cluster
apptainer pull qe-7.3.1.sif docker://nvcr.io/hpc/quantum_espresso:qe-7.3.1
apptainer exec --nv qe-7.3.1.sif pw.x -h

# Docker on a workstation with nvidia-container-toolkit
docker run --rm --gpus all -v "$PWD":/work -w /work \
    nvcr.io/hpc/quantum_espresso:qe-7.3.1 pw.x -i si_test.in
```

> [!NOTE]
> 命令中的 `pw.x -h` 只用于检查可执行文件能否启动。容器内的 MPI 与宿主机 MPI 如何配合，请以 NGC 页面上的说明为准。

### 5.3 不安装也能体验

- **Quantum Mobile**：Materials Cloud 提供的虚拟机，预装 QE 及常用工具，用 VirtualBox 打开即可。QE 用户手册在“Download”一节推荐了它。
- **nanoHUB dftqe**：<https://nanohub.org/resources/dftqe>，在浏览器里运行 QE，适合课堂演示。

---

## 6. 方式三：源码编译之 configure + make

### 6.1 前置条件

根据 QE 用户手册（Installation → Prerequisites）：

| 需求 | 最低要求 |
|---|---|
| Shell 与工具 | bash、`make`、`awk`、`sed`；脚本假设 `LC_ALL=C` |
| 构建系统 | `configure`（autoconf ≥ 2.64）或 CMake ≥ 3.20 |
| 编译器 | 支持 Fortran 2008 的 Fortran 编译器 + C 编译器 |
| MPI | 实现 MPI-3 标准（并行版必需） |
| 从 git 获取源码时 | git ≥ 2.13（用来拉取外部子模块） |
| NVIDIA GPU | NVIDIA HPC SDK ≥ 21.7 |
| 数学库（推荐） | BLAS/LAPACK（OpenBLAS 或 MKL）、FFTW3（或 MKL DFTI）、ScaLAPACK；可选 ELPA、Libxc、HDF5 |

Ubuntu 24.04 上准备一套 GNU 工具链（示例）：

```bash
sudo apt update
sudo apt install -y build-essential gfortran git cmake \
     libopenmpi-dev openmpi-bin \
     libopenblas-dev liblapack-dev libfftw3-dev libscalapack-openmpi-dev
```

集群上通常改为 `module load` 编译器、MPI 和 MKL，具体模块名问管理员。

### 6.2 获取源码

```bash
mkdir -p $HOME/src && cd $HOME/src
git clone --branch qe-7.6 --depth 1 https://gitlab.com/QEF/q-e.git qe-7.6
cd qe-7.6
```

也可以在官网的 Download 页面下载发布 tarball，然后用 `tar zxvf qe-X.Y.Z.tar.gz` 解压。编译时部分外部包会**按需联网下载**。如果计算节点不能上网，手册的做法是：把需要的包手动下载到源码的 `archive/` 子目录（**不要解压**），再运行 `make`。

### 6.3 configure（推荐在源码目录外构建）

7.6 起 autoconf 也支持源码目录外构建（out-of-source build），这样同一份源码可以同时编出 CPU 版和 GPU 版：

```bash
cd $HOME/src/qe-7.6
mkdir build-cpu && cd build-cpu
export LC_ALL=C

# MPI + OpenMP, with ScaLAPACK (default: --with-scalapack=yes)
../configure --enable-openmp --prefix=$HOME/opt/qe-7.6 \
             MPIF90=mpif90 CC=gcc
```

configure 常用选项（括号内为默认值，取自用户手册）：

| 选项 | 默认 | 说明 |
|---|---|---|
| `--enable-parallel` | yes | 检测到 MPI 时编译并行版 |
| `--enable-openmp` | no | 启用 OpenMP 线程 |
| `--with-scalapack=yes\|no\|intel` | yes | 使用 ScaLAPACK（`intel` 强制 Intel MPI + BLACS，已过时） |
| `--with-elpa-include`、`--with-elpa-lib` | 无 | ELPA 的头文件/库路径 |
| `--with-libxc`、`--with-libxc-prefix` | no | 启用 Libxc（部分 meta-GGA 需要，见 [第14章](14-hybrid-metagga.md)） |
| `--with-hdf5` | no | HDF5 输出支持 |
| `--enable-static` | no | 静态链接 |
| `--prefix=PATH` | 无 | `make install` 的安装位置 |

> [!WARNING]
> 手册特别强调：第一次运行 configure **不要加多余的环境变量**，只设置编译器和库能正常工作所需的那些。不要给编译器写完整路径（例如 `F90=/path/to/f90xyz`），否则可能出现难以理解的编译错误。用 Open MPI 的 `mpif90` 时，可以用 `OMPI_FC` 指定它调用的串行编译器。

### 6.4 检查 make.inc

configure 会生成 `make.inc`（编译规则）、`configure.msg`（配置报告）和 `config.log`（详细日志）。编译前请先看一眼：

```bash
cat configure.msg
grep -E "^(DFLAGS|MPIF90|FFLAGS|BLAS_LIBS|LAPACK_LIBS|SCALAPACK_LIBS|FFT_LIBS)" make.inc
```

要确认的内容：

| 变量 | 期望 | 说明（预处理宏取自 `include/defs.h.README`） |
|---|---|---|
| `DFLAGS` | 含 `-D__MPI`；用外部 FFTW3 时含 `-D__FFTW3`，用 MKL FFT 时含 `-D__DFTI`；用 ScaLAPACK 时含 `-D__SCALAPACK` | 只有 `-D__FFTW` 表示用的是 QE 内置的 FFTW（`defs.h.README` 中标为默认回退），生产计算建议改用外部 FFTW3 或 MKL |
| `MPIF90` | `mpif90` / `mpiifort` / `mpiifx` 等 | 串行版会是 `gfortran` 等 |
| `BLAS_LIBS`、`LAPACK_LIBS` | 指向 OpenBLAS 或 MKL | 为空或指向 QE 内置的参考实现时，性能会差很多 |
| `SCALAPACK_LIBS` | 非空 | 大体系子空间对角化需要（`-nd` 并行） |
| `FFT_LIBS` | 指向 `-lfftw3` 或 MKL | |

如果改了 `make.inc` 中的预处理或编译选项，手册要求先 `make clean` 再重新编译。

### 6.5 编译目标

```bash
make -j 8 pw          # pw.x only
make -j 8 ph pp       # phonons + post-processing
make -j 8 pwall       # pw neb ph pp pwcond (+ acfdt in the 7.6 Makefile)
make -j 8 all         # pwall cp ld1 tddfpt hp xspectra (+ gwl kcw pioud in 7.6)
make -j 8 epw         # EPW (also builds wannier90)
make install          # copy executables to $PREFIX/bin (needs --prefix)
```

直接输入 `make`（不带目标）会打印全部可用目标。可执行文件以符号链接的形式出现在构建目录的 `bin/` 下。如果 `make` 在解压外部库时被中断，再次编译可能出错，这时先在顶层目录运行 `make veryclean`（必要时 `make distclean`），然后重新 configure。

完整的示例脚本见 [`build_qe_configure.sh`](../examples/01-install/build_qe_configure.sh)（**仅供阅读和修改，本仓库从未执行过它**）。

---

## 7. 方式四：源码编译之 CMake

### 7.1 CPU 版

```bash
cd $HOME/src/qe-7.6
cmake -B build-cmake \
      -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_INSTALL_PREFIX=$HOME/opt/qe-7.6 \
      -DQE_ENABLE_OPENMP=ON \
      -DQE_ENABLE_SCALAPACK=ON \
      -DQE_FFTW_VENDOR=FFTW3
cmake --build build-cmake -j 8 --target pw ph pp
# or everything defined in 'pwall': pw ph pp pwcond neb
# cmake --build build-cmake -j 8 --target pwall
cmake --install build-cmake
ls build-cmake/bin/
```

### 7.2 常用 CMake 选项（取自 7.6 的 `CMakeLists.txt`）

| 选项 | 默认 | 说明 |
|---|---|---|
| `QE_ENABLE_MPI` | ON | MPI 并行 |
| `QE_ENABLE_OPENMP` | OFF（启用 GPU 时默认 ON） | OpenMP 线程 |
| `QE_ENABLE_SCALAPACK` | OFF | ScaLAPACK。注意 CMake 默认关闭，而 configure 默认开启 |
| `QE_ENABLE_ELPA` | OFF | ELPA 本征值求解器 |
| `QE_ENABLE_LIBXC` | OFF | Libxc 泛函库 |
| `QE_ENABLE_HDF5` | OFF | HDF5 |
| `QE_FFTW_VENDOR` | AUTO | 可选 `Intel_DFTI`、`Intel_FFTW3`、`ArmPL`、`IBMESSL`、`FFTW3`、`NVPL`、`Internal` |
| `QE_ENABLE_TEST` | ON | 生成 ctest 测试 |
| `QE_GPU` | 空 | 7.6 新增。分号分隔的 GPU 特性列表，可选 `openacc`、`openmp`、`cuda`、`rocm`、`oneapi`；`cuda`、`rocm`、`oneapi` 只能三选一 |
| `QE_GPU_ARCHS` | 空（尝试自动检测） | 7.6 新增。GPU 架构，如 `sm_80`、`sm_90`（NVIDIA）或 `gfx90a`（AMD） |
| `QE_ENABLE_MPI_GPU_AWARE` | OFF | GPU-aware（CUDA-aware）MPI |

> [!IMPORTANT]
> **7.6 的不兼容变化**：如果命令行或旧的 `CMakeCache.txt` 里还有 `QE_ENABLE_CUDA`、`QE_ENABLE_ROCM`、`QE_ENABLE_ONEAPI`、`QE_ENABLE_OFFLOAD`、`QE_ENABLE_OPENACC` 中的任何一个（**即使设成 OFF**），CMake 会直接报 FATAL_ERROR：
> `QE_ENABLE_CUDA, QE_ENABLE_ROCM, QE_ENABLE_ONEAPI, QE_ENABLE_OFFLOAD and QE_ENABLE_OPENACC CMake options are removed. Use QE_GPU and/or QE_GPU_ARCHS instead.`
> 解决方法：删掉这些选项，并**删除构建目录里的 `CMakeCache.txt`**（或干脆换一个新的构建目录）后重新配置。

### 7.3 ctest

CMake 构建默认生成测试。当前 CMake 测试只覆盖 `pw` 和 `cp` 两类，标签形如 `system--pw`：

```bash
cd build-cmake
ctest -L system--pw -j 4 --output-on-failure
```

示例脚本 [`build_qe_cmake.sh`](../examples/01-install/build_qe_cmake.sh) 把 CPU 和 GPU 两种配置写在一起，用 `MODE=cpu|gpu` 切换（**同样只供阅读，没有被执行过**）。

---

## 8. GPU 版本（NVIDIA）

### 8.1 CMake（7.6 写法）

```bash
module load nvhpc            # NVIDIA HPC SDK >= 21.7; module name is site-specific
cd $HOME/src/qe-7.6
rm -rf build-gpu
cmake -B build-gpu \
      -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_C_COMPILER=mpicc \
      -DQE_GPU="openacc;cuda" \
      -DQE_GPU_ARCHS=sm_90 \
      -DQE_ENABLE_MPI_GPU_AWARE=ON \
      -DQE_ENABLE_SCALAPACK=OFF
cmake --build build-gpu -j 16 --target pw ph pp
```

`QE_GPU_ARCHS` 要和硬件匹配：A100 对应 `sm_80`，H100 对应 `sm_90`。不确定时可以运行 `nvidia-smi --query-gpu=name,compute_cap --format=csv` 查看 compute capability。7.6 的 CMake 逻辑还规定：只给 `QE_GPU_ARCHS=sm_xx` 而不给 `QE_GPU` 时，`QE_GPU` 会自动设为 `"openacc;cuda"`。

7.5 及更早版本用 `-DQE_ENABLE_CUDA=ON`（OpenACC 随之开启）。不要把两种写法混在一起。

### 8.2 configure（GPU）

```bash
module load nvhpc
cd $HOME/src/qe-7.6
mkdir build-gpu-make && cd build-gpu-make
../configure --with-cuda=$NVHPC_CUDA_HOME --with-cuda-cc=90 \
             --with-cuda-mpi=yes --enable-openmp
make -j 16 pw
```

| 选项 | 说明（用户手册 2.4.5） |
|---|---|
| `--with-cuda=PATH` | CUDA toolkit 路径，如 `$NVHPC_CUDA_HOME` |
| `--with-cuda-cc=NN` | compute capability，要与硬件和驱动一致（如 80、90；不写时默认 35，对新 GPU 不合适） |
| `--with-cuda-runtime=X.Y` | 可选，CUDA 运行时版本 |
| `--with-cuda-mpi=yes` | 使用 CUDA-aware MPI。手册提醒：没有 NVLink 或 InfiniBand RDMA 等快速 GPU 间通信时可能在运行时崩溃 |

### 8.3 GPU 运行要点

- 一般**每个 GPU 配 1 个 MPI 进程**，优先用 k 点池（`-nk`）并行；
- `OMP_NUM_THREADS` 可设为每个进程分到的 CPU 核数；
- 详细的并行调优见 [第18章](18-automation-parallel-troubleshooting.md)。

```bash
# 1 node, 4 GPUs, 4 MPI ranks, k-point pools
export OMP_NUM_THREADS=8
mpirun -np 4 pw.x -nk 4 -i si_test.in > si_test.gpu.out
grep -i "GPU" si_test.gpu.out | head
```

---

## 9. 集群 module 与安装后检查

拿到一个现成的 QE（module、conda、别人编译的），先确认它的版本和并行能力：

```bash
module avail 2>&1 | grep -i espresso       # find available modules
module load quantum-espresso/7.6           # name is site-specific
which pw.x
ldd $(which pw.x) | grep -Ei "mkl|openblas|fftw|scalapack|elpa|mpi|cuda" 
```

真正可靠的检查是跑一个小算例，再看输出开头：

```bash
mpirun -np 4 pw.x -nk 2 -i si_test.in > si_test.out
grep -E "Program PWSCF|Parallel version|Serial version|Number of MPI processes|Threads/MPI process|K-points division|R & G space division" si_test.out
```

（示意输出，数值仅供参考）

```text
     Program PWSCF v.7.6 starts on ...
     Parallel version (MPI & OpenMP), running on       4 processor cores
     Number of MPI processes:                 4
     Threads/MPI process:                     1
     K-points division:     npool     =       2
     R & G space division:  proc/nbgrp/npool/nimage =       2
```

| 检查项 | 正常 | 异常与原因 |
|---|---|---|
| 版本号 | `v.7.6`（或你期望的版本） | 版本过旧：module/PATH 指错 |
| `Program PWSCF` 行数 | 只出现 1 次 | 出现 N 次：`mpirun` 和 pw.x 链接的 MPI 不一致，每个进程都在独立运行 |
| 进程数 | MPI + OpenMP 编译版打印 `Number of MPI processes:`，应等于 `-np`；纯 MPI 编译版只打印 `Parallel version (MPI), running on N processors` | 为 1：同上；出现 `Serial version`：pw.x 是串行编译的 |
| `K-points division: npool` | 等于 `-nk`（只有 `npool > 1` 时才打印这一行） | 不一致或缺失：命令行参数没传进去 |
| `Threads/MPI process`（仅 OpenMP 编译版打印） | 等于 `OMP_NUM_THREADS` | 线程超额订阅会严重拖慢速度 |

---

## 10. Windows：WSL2

QE 的 configure 列出了 `cygwin` 和 MinGW 交叉编译（`mingw32`/`mingw64`）等 Windows 相关架构，但这些路径的测试很少。对大多数用户来说，最省事的做法是使用 WSL2（Windows Subsystem for Linux）：

```powershell
# PowerShell (administrator)
wsl --list --online          # list installable distributions
wsl --install -d <DistroName> # e.g. one of the Ubuntu LTS names shown above
```

进入 Ubuntu 后，按 §4（conda）或 §6/§7（源码）安装即可。WSL2 下的注意事项：

| 问题 | 建议 |
|---|---|
| 计算文件放在 `/mnt/c/...` | I/O 很慢。计算目录放在 WSL 自己的文件系统（如 `~/calc`），需要时再复制回 Windows |
| Windows 编辑器保存的输入带 CRLF 换行 | 用 `file x.in` 检查；用 `dos2unix x.in` 或 `sed -i 's/\r$//' x.in` 转换。PW 用户手册 Troubleshooting 5.0.0.4 把 CRLF、Tab、非 ASCII 字符（如“智能”引号）和缺少末尾换行列为输入读取错误的常见原因 |
| 中文全角引号、全角空格 | namelist 读取失败，见 [第2章](02-input-anatomy.md) 的报错表 |
| 内存与核数 | WSL2 默认只用部分内存，可以在 Windows 用户目录的 `.wslconfig` 中调整 |
| Open MPI 版本 | IMS（日本分子科学研究所）报告过 QE 7.5 与 Open MPI 5.0.x 组合时部分测试挂起，建议优先用 Open MPI 4.1.x |

---

## 11. 验证安装

### 11.1 test-suite（make 构建）

```bash
cd $HOME/src/qe-7.6/test-suite
make pseudo                       # download PPs needed by tests into ESPRESSO_PSEUDO
make run-tests-pw NPROCS=4        # pw.x tests on 4 MPI ranks
make run-tests NPROCS=4           # all main categories (takes longer)
```

`test-suite/Makefile` 中列出的目标包括 `run-tests-XX`（XX = pw、cp、ph、pp、hp、tddfpt、kcw、epw 等）、`run-tests`、`run-custom-test testdir=DIR`、`compare`、`pseudo`、`clean`。手册说明：测试只为检查正确性，每个都很快，**不要用它们测并行性能**。串行编译的版本不要设置 `NPROCS` 大于 1。

### 11.2 官方 examples

```bash
cd $HOME/src/qe-7.6
# edit environment_variables: BIN_DIR, PSEUDO_DIR, TMP_DIR, PARA_PREFIX, PARA_POSTFIX
#   PARA_PREFIX="mpirun -np 4"
cd PW/examples/example01
./run_example
ls results/ reference/
```

`example01` 计算 Si、Al、Cu、Ni 的总能和能带。`reference/` 中是官方的参考输出。手册提醒：不同机器选择的 FFT 网格可能不同，结果会有微小差别，所以直接 `diff` 不可行，需要人工对比关键数值。`TMP_DIR` 必须是专用目录，example 脚本会清空它。

### 11.3 本仓库的冒烟测试

示例目录中的 [`si_test.in`](../examples/01-install/si_test.in) 是金刚石 Si 的 SCF（PBE，SSSP 1.3 Efficiency 的 `Si.pbe-n-rrkjus_psl.1.0.0.UPF`，`ecutwfc=30`、`ecutrho=240`，4×4×4 k 网格，`conv_thr=1.0d-8`）。2 个原子、截断和 k 网格都很小，计算量很小，适合作冒烟测试（本仓库没有运行过，不给出耗时）。赝势下载方法见 [第3章](03-pseudopotentials.md)。

```bash
cd QE/examples/01-install
mkdir -p pseudo        # put Si.pbe-n-rrkjus_psl.1.0.0.UPF here (see chapter 3)
bash check_install.sh            # only checks executables and input, does not run pw.x
RUN=1 NP=2 bash check_install.sh # also runs pw.x (you decide when to run it)
```

成功判据（逐条 grep）：

```bash
grep "convergence has been achieved" si_test.out   # SCF converged
grep "^!" si_test.out                              # final total energy line, starts with '!'
grep "highest occupied, lowest unoccupied level" si_test.out
grep "JOB DONE" si_test.out                        # clean termination
ls tmp/si_test.save/data-file-schema.xml           # XML data file written
```

（示意输出，数值仅供参考）

```text
     convergence has been achieved in   x iterations
!    total energy              =     -xx.xxxxxxxx Ry
     highest occupied, lowest unoccupied level (ev):     x.xxxx    x.xxxx
   JOB DONE.
```

同一台机器上，conda 版和源码编译版对同一输入、同一赝势得到的总能量通常只在 SCF 收敛阈值的量级上有差别（不同库、不同 FFT 网格也会带来微小差异，见手册 2.8.2）。差别明显时，先检查两边的赝势文件是否完全相同（`md5sum`），再看输出头部的 FFT 网格和版本号。

---

## 12. 学习路线与本节导航

| 阶段 | 章节 | 你会得到 |
|---|---|---|
| 入门 | 第1–4章：安装 → [输入结构](02-input-anatomy.md) → [赝势](03-pseudopotentials.md) → [收敛测试](04-convergence.md) | 能独立写出并检查一个 pw.x SCF |
| 基础性质 | [第5章 优化/EOS](05-relax-eos.md)、[第6章 能带/DOS](06-band-dos-pdos.md)、[第7章 磁性/SOC](07-magnetism-soc.md)、[第8章 DFT+U](08-dftu-hp.md) | 晶格常数、能带、磁矩 |
| 表面与界面 | [第9章 vdW/表面/功函数](09-vdw-surface-workfunction.md)、[第10章 电荷分析](10-charge-elf-bader.md)、[第11章 NEB](11-neb.md) | 吸附、功函数、反应能垒 |
| 线性响应 | [第12章 声子](12-phonon-dfpt.md)、[第13章 介电/IR/Raman](13-dielectric-ir-raman.md)、[第16章 电声/EPW](16-elph-epw.md) | 声子色散、光谱、超导 Tc |
| 进阶 | [第14章 杂化泛函](14-hybrid-metagga.md)、[第15章 Wannier90](15-wannier90.md)、[第17章 MD](17-md-cp.md) | 更准确的带隙、紧束缚模型、动力学 |
| 生产化 | [第18章 自动化/并行/排错](18-automation-parallel-troubleshooting.md) | 作业脚本、ASE/AiiDA、报错速查 |

配套脚本放在 `QE/scripts/`（如 `qe_check_input.py`、`qe_parse_out.py`、`job_qe_slurm.sh`），各章会在需要时提到。在仓库根目录下用 `python QE/scripts/<name>.py -h` 或 `bash QE/scripts/job_qe_slurm.sh -h` 查看用法。

---

## 参数选择与收敛建议

这里的“参数”指安装与运行参数，物理参数的收敛见 [第4章](04-convergence.md)。

| 场景 | 建议 |
|---|---|
| 选择安装方式 | 学习阶段用 conda；生产计算用集群 module 或自己编译（MKL/OpenBLAS + FFTW3 + ScaLAPACK） |
| MPI × OpenMP | CPU 上先用纯 MPI（`OMP_NUM_THREADS=1`）；内存紧张或进程过多时再引入 2–4 个线程 |
| `-nk`（k 点池） | 最有效的并行层级；`-nk` 应整除进程数，且不超过不可约 k 点数 |
| `-nd`（线性代数） | 大体系（几百个能带以上）且编译了 ScaLAPACK/ELPA 时才有意义 |
| GPU | 1 MPI/GPU，优先 `-nk`；架构参数 `QE_GPU_ARCHS` / `--with-cuda-cc` 要与硬件匹配 |
| `outdir` | 放在本地 scratch 或并行文件系统上，不要放在 NFS 家目录；可用环境变量 `ESPRESSO_TMPDIR` 设置默认值 |
| `pseudo_dir` | 每个输入里显式写出；也可以设置环境变量 `ESPRESSO_PSEUDO`（未设置时默认 `$HOME/espresso/pseudo/`） |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| CMake：`QE_ENABLE_CUDA, QE_ENABLE_ROCM, ... CMake options are removed. Use QE_GPU and/or QE_GPU_ARCHS instead.` | 7.6 删除了旧的 GPU 选项，命令行或 `CMakeCache.txt` 里还有残留 | 删除旧选项和 `CMakeCache.txt`（或新建构建目录），改用 `-DQE_GPU=... -DQE_GPU_ARCHS=...` |
| CMake：`Invalid QE_GPU selection ... "cuda", "rocm" and "oneapi" don't work together. Select one!` | `QE_GPU` 里同时写了多个厂商的后端 | 只保留一个，如 NVIDIA 用 `"openacc;cuda"` |
| CMake 版本过低（要求 3.20） | 系统自带的 cmake 太旧 | `pip install cmake` 或 `conda install -c conda-forge cmake` |
| configure 后 `DFLAGS` 里只有 `-D__FFTW`、`BLAS_LIBS` 为空 | 没找到外部 FFTW3/BLAS | 加载对应 module 或安装开发包，重新 configure；必要时手动编辑 `make.inc` 后 `make clean` |
| 编译时下载外部库失败 | 编译节点无法联网 | 在能联网的机器上下载，放入 `archive/`（不要解压），再 `make` |
| `pw.x: error while loading shared libraries: ... cannot open shared object file` | 运行时找不到动态库（手册 Troubleshooting 5.0.0.1） | 运行前加载与编译时相同的 module，或把库路径加入 `LD_LIBRARY_PATH` |
| 输出头部重复 N 遍，每遍的进程数都是 1 | `mpirun` 与 pw.x 链接的 MPI 不一致 | 统一 MPI：用 `ldd $(which pw.x)` 和 `which mpirun` 核对 |
| 并行运行时停在开头、不读输入（`Waiting for input...`） | 部分 MPI 不支持 `< input` 重定向 | 用 `pw.x -i x.in`（或 `-in`/`-inp`/`-input`），见手册 5.0.0.3 |
| test-suite 部分测试挂起（Open MPI 5.0.x） | 已知的组合问题报告 | 换 Open MPI 4.1.x 重新编译 |
| `apt` 装的 pw.x 读新输入报错 | Ubuntu 仓库的版本是 6.7 | 改用 conda-forge 的 7.6 或源码编译 |
| 计算奇慢，CPU 占用远超 100%×进程数 | OpenMP 线程超额订阅 | `export OMP_NUM_THREADS=1`（纯 MPI）或使进程数 × 线程数 = 物理核数 |
| `Error in routine mp_start_pools ... invalid number of pools, parent_nproc /= nproc_pool * npool` | `-nk` 不能整除 MPI 进程数 | 让 `-np` 是 `-nk` 的整数倍，例如 `-np 4 -nk 2` |
| `Error in routine read_namelists ... DFT+Hubbard input syntax has changed since v7.1`（前面还有 `WARNING!!! The input parameter lda_plus_u is obsolete.` 等） | 旧输入里有 `lda_plus_u`、`Hubbard_U(i)`、`U_projection_type` 等 7.1 前的关键词 | 删除这些关键词，改写为 `HUBBARD` 卡（[第8章](08-dftu-hp.md)） |
| `Error in routine iosys ... PPCG diagonalization not supported anymore (Dec. 2024)` | 输入里写了 `diagonalization='ppcg'` | 改用 `'david'`（默认）、`'cg'` 或 `'paro'` |

---

## 练习

1. 用 conda 安装 QE 7.6，运行 `which pw.x`，并用 `conda list` 确认它链接的是 Open MPI 还是 MPICH。
2. 阅读 [`build_qe_cmake.sh`](../examples/01-install/build_qe_cmake.sh)，把 `QE_GPU_ARCHS` 改成你所用 GPU 对应的值，并说明为什么切换到 7.6 时要删除 `CMakeCache.txt`。
3. 在 `make.inc` 中找出 `DFLAGS`，逐个解释每个 `-D__XXX` 的含义（参考源码中的 `include/defs.h.README`）。
4. 分别用 `-np 4 -nk 1`、`-np 4 -nk 2`、`-np 4 -nk 4` 运行 `si_test.in`，比较输出末尾的 `PWSCF :` 计时行（自己运行，本仓库不提供结果），并解释 `-nk` 的作用。
5. 在 WSL2 中用 Windows 记事本保存一个输入文件，用 `file` 命令检查换行符，再用 `dos2unix` 转换。

---

## 参考资料

**官方文档**

- [QE 官网](https://www.quantum-espresso.org) / [Learn 页（官方推荐的课程与学校）](https://www.quantum-espresso.org/learn/) / [Download 页](https://www.quantum-espresso.org/download-page/)
- [QE 总用户手册（安装、configure、CMake、GPU、并行）](https://www.quantum-espresso.org/Doc/user_guide/)，其中 [2.8 Running tests and examples](https://www.quantum-espresso.org/Doc/user_guide/node15.html)
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)（含 Troubleshooting）
- [GitLab Wiki：CMake 构建系统](https://gitlab.com/QEF/q-e/-/wikis/Developers/CMake-build-system)
- [QE 7.6 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.6) / [QE 7.5 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.5) / [release-notes（develop）](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes)
- [QE GitLab 主仓库](https://gitlab.com/QEF/q-e) / [users 邮件列表](https://lists.quantum-espresso.org/mailman/listinfo/users)
- [conda-forge qe 包](https://anaconda.org/conda-forge/qe) / [NVIDIA NGC QE 容器](https://catalog.ngc.nvidia.com/orgs/hpc/containers/quantum_espresso) / [Microsoft WSL 安装说明](https://learn.microsoft.com/en-us/windows/wsl/install)

**课程与教程（推荐，链接而非复制）**

- [Cornell PARADIM Summer School 2018](https://www.paradim.org/cu_ss_2018)：Giustino、Pamuk 的 QE 上机课（编译与 Si 示例、声子、能带），每讲有视频。原理讲解很好，但基于 QE 6.x 与 LDA，下载/编译部分已过时。
- [Cornell PARADIM Computational Tutorials](https://www.paradim.org/toolbox/computation/tutorials)：QE 配置、输入与收敛参数、平衡结构等系列。
- [Imperial College MSE404-MM Lab 1: Getting Started in Linux](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab01/readme/)：零基础的 Linux 与 QE 环境准备，持续维护。
- [pranabdas/espresso：安装](https://pranabdas.github.io/espresso/setup/install)：从安装到 EPW 的系统教程中的安装一节。
- [ENCCS: Efficient materials modelling on HPC](https://enccs.github.io/efficient-materials-modelling-on-hpc/)：`-nk`/`-nd`/`-nt` 与 GPU 运行，适合配合本章 §8–§9 阅读。
- [MaX e-School 2021（ICTP）](https://indico.ictp.it/event/9616/)：两周课程与视频。

**相关仓库（仅链接）**

- [QEF/q-e](https://github.com/QEF/q-e)：官方 GitHub 镜像（开发在 GitLab），`PW/examples`、`test-suite` 是权威示例。
- [hfp/xconfigure](https://github.com/hfp/xconfigure)：Intel 平台编译 QE 的配方。
- [pranabdas/espresso](https://github.com/pranabdas/espresso)：配套网站的源文件。
- [ImperialCollegeLondon/MSE404-MM](https://github.com/ImperialCollegeLondon/MSE404-MM)：课程仓库（CC BY-NC-SA 4.0）。

**关键文献**

- P. Giannozzi et al., *J. Phys.: Condens. Matter* 21, 395502 (2009). DOI: [10.1088/0953-8984/21/39/395502](https://doi.org/10.1088/0953-8984/21/39/395502)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)
- P. Giannozzi et al., *J. Chem. Phys.* 152, 154105 (2020)（QE toward the exascale）. DOI: [10.1063/5.0005082](https://doi.org/10.1063/5.0005082)
- I. Carnimeo et al., *J. Chem. Theory Comput.* 19, 6992 (2023)（GPU 与 exascale）. DOI: [10.1021/acs.jctc.3c00249](https://doi.org/10.1021/acs.jctc.3c00249)

---

[返回目录](../README.md) | [下一章 →](02-input-anatomy.md)
