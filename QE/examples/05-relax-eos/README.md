# 示例 05：Si 变胞优化、CO 分子优化与 Cu 状态方程

配套教程：[第5章 结构优化与状态方程](../../tutorials/05-relax-eos.md)

## 目的

- `si.vc-relax.in`：从略微膨胀的晶胞出发，用 `vc-relax` 优化金刚石 Si 的晶格常数，观察末尾的 "final scf" 与 Pulay 应力；
- `co.relax.in`：固定晶胞的 `relax`，优化盒子中 CO 分子的键长，演示 `if_pos` 冻结原子；
- `cu_eos/`：fcc Cu 的 E(a) 扫描，用 `ev.x`（或 Python）拟合 Birch-Murnaghan 状态方程，得到 a0、B0、B0'。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.vc-relax.in` | Si，`ibrav=2`，`cell_dofree='ibrav'`，`press=0`，`press_conv_thr=0.2` kbar |
| `co.relax.in` | CO 分子，12 Å 立方盒子，`K_POINTS gamma`，C 原子 `if_pos = 0 0 0` |
| `job_relax_slurm.sh` | Slurm 作业脚本；`INPUT=co.relax.in NK=1 sbatch job_relax_slurm.sh` 可切换输入 |
| `cu_eos/cu_eos_template.in` | Cu 单点模板，占位符 `@ALAT@`（Å，用 `A` 关键词） |
| `cu_eos/run_eos.sh` | 生成 `a_<值>/pw.in`，`RUN=1` 时运行，汇总 `cu_eos.dat` 并调用 `ev.x` |

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 获取（文件名区分大小写）：

| 元素 | 文件名 | SSSP 提示截断（Ry） |
|---|---|---|
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | 30 / 240 |
| C | `C.pbe-n-kjpaw_psl.1.0.0.UPF` | 45 / 360 |
| O | `O.pbe-n-kjpaw_psl.0.1.UPF` | 50 / 400 |
| Cu | `Cu.paw.z_11.ld1.psl.v1.0.0-low.upf` | 90 / 720 |

放到 `05-relax-eos/pseudo/`（`cu_eos/` 下也建一个 `pseudo` 目录或软链接：`ln -s ../pseudo cu_eos/pseudo`）。

## 运行步骤

```bash
cd QE/examples/05-relax-eos
mkdir -p pseudo                        # 放入上表中的 UPF

# Si 变胞优化
mpirun -np 4 pw.x -nk 4 -in si.vc-relax.in > si.vc-relax.out
# CO 分子优化（Gamma 点，不用 -nk）
mpirun -np 4 pw.x -in co.relax.in > co.relax.out
# 或者提交 Slurm
sbatch job_relax_slurm.sh
INPUT=co.relax.in NK=1 sbatch job_relax_slurm.sh

# Cu 状态方程
cd cu_eos && ln -s ../pseudo pseudo
RUN=1 PW_CMD="mpirun -np 8 pw.x -nk 4" bash run_eos.sh
```

## 预期结果 / 如何检查

```bash
grep -e 'bfgs converged' -e 'End of BFGS' -e 'Final enthalpy' si.vc-relax.out
grep -A12 'Begin final coordinates' si.vc-relax.out
grep 'P=' si.vc-relax.out | tail -2        # 最后两次：最后一个 BFGS 步与 final scf
grep -e 'Final energy' co.relax.out
grep -A4 'ATOMIC_POSITIONS' co.relax.out | tail -4
head -3 cu_eos/cu_eos.ev.out               # ev.x 拟合结果（a0、k0、dk0）
```

成功标准：

1. 输出中有 `bfgs converged in ... scf cycles and ... bfgs steps` 与 `JOB DONE`；
2. Si：final scf 的压强与最后一个 BFGS 步差别不大（< 1 kbar 量级）；差别大说明截断不够（Pulay 应力），提高 `ecutwfc` 重做；
3. CO：O 原子受力分量都小于 `forc_conv_thr`；键长与实验值（约 1.128 Å）同量级（PBE 结果会有偏差，以你的计算为准）；
4. Cu：`cu_eos.ev.out` 中每个点的 `E_diff` 很小，拟合的 a0 落在扫描区间内部（若在边界，扩大扫描范围）。

本仓库不提供真实计算结果，所有数值以你的计算为准。

## 注意事项

- 所有输入为 ASCII、LF 换行。
- `etot_conv_thr`、`forc_conv_thr` 单位分别是 Ry 与 Ry/Bohr；`press`、`press_conv_thr` 单位是 kbar。
- `ev.x` 是交互式程序，`run_eos.sh` 用 here-document 依次回答：单位（`Ang`）、晶格类型（`fcc`）、状态方程类型（1–4）、输入文件、输出文件。
- Cu 的 k 网格和 `degauss` 只是示例，需按 [第4章](../../tutorials/04-convergence.md) 测试。
