# 示例 11：H 原子在 Al(100) 表面相邻四重空位之间扩散的 CI-NEB

配套教程：[第11章 过渡态搜索：neb.x](../../tutorials/11-neb.md)

## 目的

用一个小体系走完 neb.x 的完整流程：

1. 用**同一套** `&SYSTEM`/`&ELECTRONS`/k 点分别弛豫初态（H 位于空位 A）和末态（H 位于相邻空位 B）；
2. 把弛豫后的坐标填入 `h_al100.neb.in`，运行 climbing-image NEB（`CI_scheme='auto'`）；
3. 从 `h_al100_neb.dat`、`.int`、`.xyz`/`.axsf` 读取能垒和鞍点结构；
4. 在鞍点处做只位移 H 的 Γ 点 ph.x（`nat_todo=1`）+ dynmat.x，确认恰好有一个虚频（快速近似检查）。

体系：Al(100) 2×2 表面晶胞、3 层 12 个 Al（底层 4 个固定）+ 1 个 H，PBE，SSSP 1.3 赝势，6×6×1 k 点，Marzari–Vanderbilt smearing。这是教学规模的模型，正式研究需要更厚的 slab、更大的表面晶胞和收敛测试。

## 文件说明

| 文件 | 说明 |
|---|---|
| `h_al100.initial.relax.in` | 初态弛豫：H 在 (1/4, 1/4) 空位，底层 Al `if_pos = 0 0 0` |
| `h_al100.final.relax.in` | 末态弛豫：H 在 (3/4, 1/4) 空位，其余设置与初态逐字相同 |
| `h_al100.neb.in` | neb.x 输入：`BEGIN ... END` 超级卡片，7 个 image，`opt_scheme='broyden'`，`CI_scheme='auto'` |
| `h_ts.scf.in` | 鞍点结构的 SCF（`conv_thr=1d-12`），坐标需替换为爬坡 image 的坐标 |
| `h_ts.ph_gamma.in` | ph.x，q = Γ，`nat_todo=1`，只位移第 13 个原子（H） |
| `h_ts.dynmat.in` | dynmat.x，`remove_interaction_blocks=.true.` 去掉与未位移原子的耦合；频率写入 `filout='h_ts.modes'` |
| `job_neb_slurm.sh` | Slurm 脚本：`STEP=relax` 弛豫两个端点，`STEP=neb` 运行 neb.x（`-ni` 为 image 并行组数） |

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载，放到本目录的 `pseudo/`：

| 元素 | 文件名 | SSSP 提示截断（Ry） |
|---|---|---|
| Al | `Al.pbe-n-kjpaw_psl.1.0.0.UPF` | 30 / 240 |
| H | `H.pbe-rrkjus_psl.1.0.0.UPF` | 60 / 480 |

输入中取两者的最大值 `ecutwfc=60`、`ecutrho=480`。整包下载方法见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/11-neb
mkdir -p pseudo                                   # 放入上表中的 UPF

# 1. 端点弛豫（同一套设置）
mpirun -np 16 pw.x -nk 4 -in h_al100.initial.relax.in > h_al100.initial.relax.out
mpirun -np 16 pw.x -nk 4 -in h_al100.final.relax.in   > h_al100.final.relax.out
grep -A15 'Begin final coordinates' h_al100.initial.relax.out
grep -A15 'Begin final coordinates' h_al100.final.relax.out
#    -> 把两组 ATOMIC_POSITIONS 分别粘贴到 h_al100.neb.in 的 FIRST_IMAGE / LAST_IMAGE

# 2. CI-NEB：neb.x 不读标准输入，必须用 -inp
mpirun -np 40 neb.x -ni 5 -inp h_al100.neb.in > h_al100.neb.out
#    或 Slurm：
STEP=relax sbatch job_neb_slurm.sh
STEP=neb NI=5 sbatch job_neb_slurm.sh

# 3. 鞍点的虚频检查（把爬坡 image 的坐标填进 h_ts.scf.in 之后）
mpirun -np 16 pw.x -nk 4 -in h_ts.scf.in > h_ts.scf.out
mpirun -np 16 ph.x -nk 4 -in h_ts.ph_gamma.in > h_ts.ph_gamma.out
dynmat.x < h_ts.dynmat.in > h_ts.dynmat.out
```

## 预期结果 / 如何检查

```bash
grep -e 'bfgs converged' -e '^!' h_al100.initial.relax.out | tail -2
grep -e 'bfgs converged' -e '^!' h_al100.final.relax.out   | tail -2
grep 'activation energy' h_al100.neb.out | tail -2          # 正向 (->) 与逆向 (<-) 能垒，eV
grep -e 'neb: convergence achieved' -e 'reached the maximum number' h_al100.neb.out
grep 'climbing image' h_al100.neb.out | tail -1
cat h_al100_neb.dat                                          # 反应坐标(归一化)  E-E(1) (eV)  error (eV/A)
ls tmp/h_al100_neb_*/PW.out                                  # 每个 image 的 pw.x 输出
grep 'freq' h_ts.modes                                       # 39 行：36 个 0（被置零的 Al），H 的 3 个模式中应有 1 个为负
```

成功标准：

1. 两个端点都 `bfgs converged`；由于两个空位等价，两端能量应几乎相同；
2. `neb: convergence achieved in ... iterations`，最后一步所有 image 的 `error (eV/A)` 小于 `path_thr`；
3. 能量剖面只有一个峰，最高点在路径中间（image 4），对应 H 越过两个空位之间的桥位；
4. 正向、逆向能垒相近（对称路径）；
5. dynmat.x 的 `h_ts.modes` 中共有 3·nat = 39 个频率：36 个为 0（被置零的 Al 行列，无物理意义），H 的 3 个模式中恰好 1 个为负（虚频，方向沿 A→B），其余 2 个为正。

本仓库不提供参考输出；能垒数值取决于 slab 厚度、k 点、截断能等，以你自己的计算为准。

## 注意事项

- 端点与 NEB 必须使用**完全相同**的 `&SYSTEM`、`&ELECTRONS`、k 点、晶胞和固定原子；原子顺序必须一致。
- `if_pos` 只在 `FIRST_IMAGE` 中读取（`LAST_IMAGE` 中可以不写），NEB 中所有 image 使用同样的固定约束。
- `nstep_path` 默认值是 1，不写的话 neb.x 只走一步就结束。
- 中断后续算：把 `restart_mode` 改为 `'restart'`，保留 `h_al100_neb.path` 和 `tmp/` 目录重新提交。
- `nat_todo` 只位移部分原子是近似（INPUT_PH 原话："Use with care, at your own risk"），只能用来确认虚频个数和方向；正式的振动/零点能修正需要更完整的计算。
- 所有输入为 ASCII、LF 换行。
