# 示例 01：QE 安装脚本与 Si 冒烟测试（Installation Recipes & Si Smoke Test）

配套教程：[第1章 QE 概述、版本与安装](../../tutorials/01-install.md)

## 目的

- 提供 Quantum ESPRESSO 7.6 的两份**示例编译脚本**（CMake 与 configure + make，分 CPU/GPU 两种模式），供阅读和按本机情况修改；
- 用一个几秒钟就能完成的金刚石 Si SCF（PBE，超软赝势）确认 pw.x 可以运行：SCF 能收敛、MPI 配对正确、XML 数据文件能正常写出。

> [!WARNING]
> 本仓库只整理教程，**没有执行过这里的任何脚本，也没有运行过 pw.x**。编译脚本中的 module 名、编译器和库路径需要按你的机器修改。

## 文件说明

| 文件 | 说明 |
|---|---|
| `build_qe_cmake.sh` | CMake 编译示例：`MODE=cpu`（MPI + OpenMP + ScaLAPACK）或 `MODE=gpu`（`QE_GPU="openacc;cuda"` + `QE_GPU_ARCHS`，QE 7.6 写法）；每次都在新的构建目录中配置，避免旧 `CMakeCache.txt` 造成的问题 |
| `build_qe_configure.sh` | configure + make 编译示例（源码目录外构建，7.6 起支持）；`CONFIGURE_ONLY=1` 时只配置并打印 `make.inc` 关键行 |
| `check_install.sh` | 安装检查。默认**只做静态检查**（可执行文件、链接库、输入换行符、赝势是否存在）；`RUN=1` 时才运行 pw.x 并检查输出 |
| `si_test.in` | 冒烟测试输入：金刚石 Si，`ibrav=2`，`celldm(1)=10.26` Bohr，`ecutwfc=30`/`ecutrho=240` Ry，4×4×4 k 网格 |

## 赝势准备（不提交 UPF 文件）

本例只需要一个文件，来自 SSSP 1.3 PBE Efficiency：

| 元素 | 文件名 | 类型 | SSSP 推荐截断（Ry） |
|---|---|---|---|
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | 超软（pslibrary 1.0.0） | 30 / 240 |

```bash
cd QE/examples/01-install
mkdir -p pseudo
# Option A: single file from the QE pseudopotential site
wget -P pseudo https://pseudopotentials.quantum-espresso.org/upf_files/Si.pbe-n-rrkjus_psl.1.0.0.UPF
# Option B: take it from the SSSP_1.3.0_PBE_efficiency.tar.gz archive (see chapter 3)
```

两种来源的文件应当相同，可以用 SSSP JSON 中给出的 md5 核对，详见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/01-install

# 0) (optional) read and adapt a build script, then run it yourself
less build_qe_cmake.sh
# MODE=cpu bash build_qe_cmake.sh

# 1) static checks only (no calculation)
bash check_install.sh

# 2) when you are ready, run the smoke test
RUN=1 NP=2 bash check_install.sh
RUN=1 NP=4 NK=2 bash check_install.sh                 # with k-point pools
RUN=1 NP=4 NK=2 MPI_LAUNCHER=srun bash check_install.sh   # inside a Slurm allocation

# 3) or run by hand
mpirun -np 2 pw.x -i si_test.in > si_test.out
```

在集群上请在计算节点上运行（交互作业或批处理作业），作业脚本写法见 [公共篇第2章](../../../Common/tutorials/02-job-scheduler.md)。

## 预期结果 / 如何检查

`check_install.sh` 最后一行为 `RESULT: PASS` 表示通过。手动检查：

```bash
grep "Program PWSCF" si_test.out                     # version line
grep -c "Program PWSCF" si_test.out                  # must be 1 (N means N independent serial runs)
grep -E "Number of MPI processes|Parallel version" si_test.out   # process count must equal -np
grep "convergence has been achieved" si_test.out     # SCF converged
grep "^!" si_test.out                                # final total energy (Ry)
grep "highest occupied, lowest unoccupied level" si_test.out
grep "JOB DONE" si_test.out
ls tmp/si_test.save/data-file-schema.xml
```

（示意输出，数值仅供参考）

```text
!    total energy              =     -xx.xxxxxxxx Ry
     highest occupied, lowest unoccupied level (ev):     x.xxxx    x.xxxx
   JOB DONE.
```

- 绝缘体使用默认的 `occupations='fixed'`，所以输出的是 HOMO/LUMO，而不是费米能。
- 同一台机器上，不同安装（conda、源码）得到的总能量应当在 `conv_thr` 的量级内一致。

## 注意事项

- `outdir=./tmp` 在重跑时会被复用。为了得到干净的测试，重跑前删除 `tmp/` 和旧的 `si_test.out`（脚本在 `RUN=1` 时会自动删除）。
- 如果输出头部（`Program PWSCF` 行）出现多次、每次报告的进程数都是 1，说明启动用的 `mpirun` 与 pw.x 链接的 MPI 不一致。注意：只有 MPI + OpenMP 编译版才打印 `Number of MPI processes:`，纯 MPI 编译版打印的是 `Parallel version (MPI), running on N processors`。
- Windows/WSL2 用户：输入文件必须是 LF 换行，脚本会检查 CRLF。
- GPU 版本：每个 GPU 配 1 个 MPI 进程，`QE_GPU_ARCHS`（或 `--with-cuda-cc`）必须与硬件匹配。
