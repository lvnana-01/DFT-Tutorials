# 示例 04：Si 与 Al 的截断能、k 点与 smearing 收敛测试

配套教程：[第4章 收敛测试：ecutwfc/ecutrho、k 点与 smearing](../../tutorials/04-convergence.md)

## 目的

- 对半导体 Si（超软赝势）扫描 `ecutwfc`（固定 dual = `ecutrho/ecutwfc` = 8）和 Monkhorst-Pack k 网格，以 ΔE < 1 meV/atom 与压强变化 < 0.5 kbar 为判据；
- 对金属 Al（PAW）做 **k 网格 × `degauss`** 的二维扫描，理解为什么金属要把 k 点和 smearing 一起测试。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si_ecut_template.in` | Si 金刚石结构（`ibrav=2`，`celldm(1)=10.26` Bohr）模板，占位符 `@ECUT@`、`@ECUTRHO@`、`@K@` |
| `al_kpt_template.in` | fcc Al（`ibrav=2`，`celldm(1)=7.65` Bohr）模板，占位符 `@K@`、`@DEGAUSS@`、`@SMEAR@` |
| `run_conv.sh` | 生成各测试子目录与 `pw.in`；`RUN=1` 时调用 pw.x；最后汇总 `summary.dat` |

占位符写法与仓库脚本 `QE/scripts/qe_conv_test.sh` 一致（`@SMEAR@` 是本例额外的），两者都可以用。注意通用脚本的 dual 默认是 4，用于本例的超软/PAW 赝势时要加 `-d 8`；通用脚本用 `bash ../../scripts/qe_conv_test.sh -h` 查看用法。

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载以下两个文件，放到本目录的 `pseudo/` 子目录：

```bash
mkdir -p pseudo
# 下载 SSSP_1.3.0_PBE_efficiency.tar.gz 后解压，只保留需要的文件：
#   Si.pbe-n-rrkjus_psl.1.0.0.UPF   (SSSP 提示截断 30 / 240 Ry)
#   Al.pbe-n-kjpaw_psl.1.0.0.UPF    (SSSP 提示截断 30 / 240 Ry)
ls pseudo/
```

赝势来源与下载细节见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/04-convergence

# 1) 只生成输入，检查目录结构
bash run_conv.sh ecut
ls conv_ecut/            # ecut_25 ecut_30 ... labels.txt summary.dat

# 2) 真正运行（本机 4 核示例）
RUN=1 PW_CMD="mpirun -np 4 pw.x" bash run_conv.sh ecut
RUN=1 PW_CMD="mpirun -np 4 pw.x" ECUT_FIX=40 bash run_conv.sh kpt
RUN=1 PW_CMD="mpirun -np 4 pw.x" bash run_conv.sh alsmear

# 3) 只重新汇总（例如作业在集群上跑完后）
bash run_conv.sh collect conv_ecut
```

集群上可把第 2 步的命令写进 Slurm 脚本（模板见 `QE/scripts/job_qe_slurm.sh`，用法见 [Common 第2章](../../../Common/tutorials/02-job-scheduler.md)）。

## 预期结果 / 如何检查

`summary.dat` 的格式（**示意输出，数值仅供参考**，不是真实计算结果）：

```text
# label  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf  status
ecut_25        -15.8xxxxxxx        x.xxx      -xx.xx      8 done
...
ecut_80        -15.8xxxxxxx        0.000       -x.xx      8 done
```

检查要点：

1. 每个子目录的 `pw.out` 都有 `JOB DONE`，`status` 列为 `done`；
2. Si：`dE_vs_last` 随 `ecutwfc` 增大单调减小，取 |dE| < 1 meV/atom 且压强变化 < 0.5 kbar 的最小值并留余量；
3. Al：在每个 `degauss` 下找到 k 收敛点，再比较不同 `degauss` 之间的能量差与 `smearing contrib. (-TS)` 的大小；
4. 能量差要用**同一组赝势、同一 `ecutrho`** 的计算比较。

## 注意事项

- 所有输入为 ASCII、LF 换行；在 Windows 上编辑后请确认没有被改成 CRLF。
- `conv_thr = 1.0d-10` 远严于生产计算，目的是让 SCF 噪声远小于 1 meV/atom 的判据。
- Al 扫描默认 4 个 `degauss` × 8 个网格 = 32 个计算，k = 24 的计算较慢，可先删减 `KS`。
- 本仓库不提供任何真实计算结果，所有数值以你自己的测试为准。
