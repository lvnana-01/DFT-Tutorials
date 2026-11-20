# 示例 12：Si 与 Al 的 DFPT 声子（ph.x → q2r.x → matdyn.x）

配套教程：[第12章 声子：ph.x、q2r.x 与 matdyn.x](../../tutorials/12-phonon-dfpt.md)

## 目的

- **Si（半导体，非极性）**：先做 Γ 点单 q 计算熟悉 ph.x 输出，再在 4×4×4 q 网格上做 DFPT，用 q2r.x 得到实空间力常数，用 matdyn.x 插值出 Γ–X–W–X'–K–Γ–L 色散和声子 DOS；
- **Al（金属）**：演示金属声子计算的要点：SCF 必须使用 smearing，k 点更密，声子对 k 点与 `degauss` 收敛较慢。

极性材料（GaAs、MgO）的 LO–TO 劈裂需要介电张量和 Born 有效电荷，放在 [第13章](../../tutorials/13-dielectric-ir-raman.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.scf.in` | Si 基态，`ibrav=2`，`conv_thr=1d-12`，8×8×8 k 点 |
| `si.ph_gamma.in` | 单个 q = Γ 的 ph.x 计算，`fildyn='si.dynG'` |
| `si.ph.in` | `ldisp=.true.`，`nq1=nq2=nq3=4`，`tr2_ph=1d-16`，`fildyn='si.dyn'` |
| `si.q2r.in` | `fildyn='si.dyn'` → `flfrc='si444.fc'`，`zasr='crystal'` |
| `si.matdyn.in` | `q_in_band_form=.true.`，6 个高对称点（2π/a 单位），输出 `si.freq` 和 `si.freq.gp` |
| `si.matdyn_dos.in` | `dos=.true.`，20×20×20 q 网格，输出 `si.phdos` |
| `al.scf.in` | fcc Al，`occupations='smearing'`，`smearing='mv'`，`degauss=0.02`，16×16×16 k 点 |
| `al.ph.in` / `al.q2r.in` / `al.matdyn.in` | Al 的 4×4×4 q 网格、力常数与色散 |
| `run_si.sh` | Si 全流程命令清单；默认只打印命令，`RUN=1` 时才执行 |
| `job_ph_slurm.sh` | Slurm 脚本：`SYS=si` 或 `SYS=al`；`NI>1` 时用 image 并行并自动追加 `recover=.true.` 的收集步骤 |

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载，放到本目录的 `pseudo/`：

| 元素 | 文件名 | SSSP 提示截断（Ry） |
|---|---|---|
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | 30 / 240 |
| Al | `Al.pbe-n-kjpaw_psl.1.0.0.UPF` | 30 / 240 |

SSSP 的截断提示已经考虑了声子频率的收敛（SSSP 验证项之一），但你的体系仍应按 [第4章](../../tutorials/04-convergence.md) 检查。下载方法见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/12-phonon-dfpt
mkdir -p pseudo                                   # 放入上表中的 UPF

# Si：Gamma 点单 q（几分钟量级的小计算）
mpirun -np 8 pw.x -nk 4 -in si.scf.in > si.scf.out
mpirun -np 8 ph.x -nk 4 -in si.ph_gamma.in > si.ph_gamma.out

# Si：完整色散
bash run_si.sh                                    # 先看命令（dry run）
RUN=1 bash run_si.sh                              # 真正执行
python ../../scripts/qe_phonon_plot.py -h         # 查看画图脚本用法

# Al：金属
mpirun -np 8 pw.x -nk 4 -in al.scf.in > al.scf.out
mpirun -np 8 ph.x -nk 4 -in al.ph.in > al.ph.out
q2r.x -in al.q2r.in > al.q2r.out
matdyn.x -in al.matdyn.in > al.matdyn.out

# Slurm
SYS=si sbatch job_ph_slurm.sh
NI=4 SYS=al sbatch job_ph_slurm.sh
```

## 预期结果 / 如何检查

```bash
grep -A10 'Diagonalizing the dynamical matrix' si.ph_gamma.out    # 6 个 freq 行
grep -A12 'uniform grid of q-points' si.ph.out                     # 4x4x4 网格的不可约 q 点列表
grep -c 'Convergence has been achieved' si.ph.out                  # 每个不可约表示一次
grep 'No convergence has been achieved' si.ph.out                  # 应为空
ls si.dyn*                                                          # si.dyn0 ... si.dynN
head -3 si.freq.gp                                                  # 第 1 列为路径长度，其后 6 列为频率 (cm^-1)
awk '{for(i=2;i<=NF;i++) if ($i < -1.0) print NR": "$0}' si.freq.gp | head   # 负频率（虚频）检查，应为空
```

成功标准：

1. 所有输出以 `JOB DONE` 结束，ph.x 输出中没有 `No convergence has been achieved`；
2. Si 的 Γ 点：3 个声学模频率接近 0（DFPT 中通常在几个 cm⁻¹ 以内），3 个光学模简并；
3. 色散中没有虚频（`si.freq.gp` 中没有明显的负值），Γ 点附近声学支线性趋于 0；
4. 声子 DOS（`si.phdos`）对频率积分等于 3 × 原子数（matdyn.x 文档：DOS 归一化到 3·nat）；
5. Al 的色散只有 3 支声学支，没有虚频；改变 k 点或 `degauss` 时频率变化小于你设定的容差（例如 2–3%）。

本仓库不提供参考输出；频率数值请以你自己的计算为准，可与实验中子散射数据对比。

## 注意事项

- ph.x 读取的是 SCF 写在 `outdir/prefix.save` 中的数据：`prefix`、`outdir` 必须一致；SCF 不能用 `K_POINTS gamma`（Γ 点技巧不支持 ph.x）。
- ph.x 不能直接接在 `vc-relax` 后面运行（报错 `The phonon code is not working after vc-relax`）：用优化后的晶格常数重新跑一次 `scf`。
- 尽量保留 `ibrav` 和 Wyckoff 位置，不要把高对称结构写成 `ibrav=0` 加带噪声的坐标，否则可能出现对称性相关的报错（PHonon 用户手册 Troubleshooting）。
- 金属声子对 k 点和 `degauss` 收敛很慢，Al 的 16×16×16 + 0.02 Ry 只是演示值。
- 中断后在 `&INPUTPH` 中加 `recover = .true.` 用同一输入重新运行即可续算。
- 所有输入为 ASCII、LF 换行。
