# 示例 10：电荷密度、差分电荷、ELF、STM 与 Bader 分析

配套教程：[第10章 电荷密度、差分电荷、ELF、STM 与 Bader 分析](../../tutorials/10-charge-elf-bader.md)

## 目的

- `si.scf.in` + `si_charge.pp.in` / `si_charge_line.pp.in`：体相 Si 的价电子密度，导出 3D cube 和沿 Si–Si 键的 1D 曲线，熟悉 pp.x 的"两步"结构（`&INPUTPP` 提取 → `&PLOT` 输出）；
- `h2o.scf.in` + `h2o_elf.pp.in` / `h2o_rho.pp.in`：用模守恒（NC）赝势计算孤立 H2O 分子的 ELF（`plot_num=8`）和密度；
- `co_on_cu_diff/`：CO 吸附在 Cu(111) 2×2 三层 slab 上，演示差分电荷 Δρ = ρ(AB) − ρ(slab) − ρ(CO)、PAW 全电子密度（`plot_num=21`）+ Henkelman `bader` 程序、Tersoff–Hamann STM（`plot_num=5`）和 projwfc.x 的 Löwdin 电荷。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.scf.in` | Si 金刚石结构 SCF（`ibrav=2`，US 赝势） |
| `si_charge.pp.in` | `plot_num=0` → `si.rho.pp` → `si_rho.cube`（`iflag=3`，`output_format=6`） |
| `si_charge_line.pp.in` | `&INPUTPP` 留空，只读 `si.rho.pp`，沿 [-1 1 1] 键方向输出 1D 曲线 `si_rho_bond.dat` |
| `h2o.scf.in` | 12 Å 立方盒子中的 H2O，`K_POINTS gamma`，`assume_isolated='mt'`，SG15 ONCV 赝势 |
| `h2o_elf.pp.in` | `plot_num=8`（ELF）→ `h2o_elf.cube` |
| `h2o_rho.pp.in` | `plot_num=0` → `h2o_rho.cube`（可直接交给 `bader`，NC 赝势只有价电子密度，见注意事项） |
| `co_on_cu_diff/ab.scf.in` | CO/Cu(111)（14 原子，未弛豫的初猜结构） |
| `co_on_cu_diff/slab.scf.in` | 只有 Cu slab，坐标与 AB 中**完全相同** |
| `co_on_cu_diff/co.scf.in` | 只有 CO，晶胞、坐标、截断、k 点与 AB 完全相同 |
| `co_on_cu_diff/rho_{ab,slab,co}.pp.in` | 三次 `plot_num=0`，只做第一步，写出 `*.rho.pp` |
| `co_on_cu_diff/diff_plot.pp.in` | `nfile=3`，`weight = 1, -1, -1` 线性组合 → `drho.cube` |
| `co_on_cu_diff/ab_rho.pp.in` | 把 `ab.rho.pp` 写成 `ab_rho.cube`（Bader 的 chargefile） |
| `co_on_cu_diff/ab_ae.pp.in` | `plot_num=21`：PAW 全电子（价 + 芯）密度 → `ab_ae.cube`（Bader 参考密度） |
| `co_on_cu_diff/stm.pp.in` | `plot_num=5`，`sample_bias=-0.0735` Ry（约 −1 eV）→ `ab_stm.cube` |
| `co_on_cu_diff/ab.projwfc.in` | projwfc.x，输出 Löwdin 电荷和 spilling parameter |
| `co_on_cu_diff/run_diff.sh` | 整个流程的命令清单；默认只打印命令，`RUN=1` 时才执行 |

## 赝势准备（不提交 UPF 文件）

| 元素 | 文件名 | 来源 | 用于 |
|---|---|---|---|
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency)（提示 30 / 240 Ry） | `si.scf.in` |
| Cu | `Cu.paw.z_11.ld1.psl.v1.0.0-low.upf` | SSSP 1.3 PBE Efficiency（90 / 720 Ry） | `co_on_cu_diff/` |
| C | `C.pbe-n-kjpaw_psl.1.0.0.UPF` | SSSP 1.3 PBE Efficiency（45 / 360 Ry） | `co_on_cu_diff/` |
| O | `O.pbe-n-kjpaw_psl.0.1.UPF` | SSSP 1.3 PBE Efficiency（50 / 400 Ry） | `co_on_cu_diff/` |
| O | `O_ONCV_PBE-1.2.upf` | [SG15 ONCV](http://www.quantum-simulation.org/potentials/sg15_oncv/)（`sg15_oncv_upf_2020-02-06.tar.gz`） | `h2o.scf.in` |
| H | `H_ONCV_PBE-1.2.upf` | SG15 ONCV | `h2o.scf.in` |

全部放到本目录的 `pseudo/` 子目录（`co_on_cu_diff/` 中的输入使用 `../pseudo`）。SSSP 的整包下载方法见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/10-charge-elf-bader
mkdir -p pseudo                                   # 放入上表中的 UPF

# (1) Si 电荷密度
mpirun -np 4 pw.x -nk 4 -in si.scf.in > si.scf.out
mpirun -np 4 pp.x -in si_charge.pp.in > si_charge.pp.out
pp.x -in si_charge_line.pp.in > si_charge_line.pp.out

# (2) H2O 的 ELF 与密度（Gamma 点，不用 -nk）
mpirun -np 4 pw.x -in h2o.scf.in > h2o.scf.out
mpirun -np 4 pp.x -in h2o_elf.pp.in > h2o_elf.pp.out
mpirun -np 4 pp.x -in h2o_rho.pp.in > h2o_rho.pp.out

# (3) CO/Cu(111)：先检查命令，再执行
cd co_on_cu_diff
bash run_diff.sh                                   # 只打印命令（dry run）
RUN=1 PW="mpirun -np 16 pw.x -nk 4" PP="mpirun -np 16 pp.x" bash run_diff.sh
```

在 Slurm 集群上可以把上述命令放进 [`job_qe_slurm.sh`](../../scripts/job_qe_slurm.sh) 模板（用 `bash ../../scripts/job_qe_slurm.sh -h` 查看用法）。

## 预期结果 / 如何检查

```bash
grep -e '^!' -e 'JOB DONE' si.scf.out h2o.scf.out co_on_cu_diff/*.scf.out
grep 'Dense  grid' co_on_cu_diff/*.scf.out     # 三行的 FFT dimensions 必须完全一致
head -6 si_rho.cube                            # cube 头：注释 2 行、原子数+原点、三行网格向量
grep -i 'not fully' h2o_elf.pp.out   # 此提示总会打印；对 NC 赝势无影响，对 US/PAW 才意味着 ELF 不完整
head co_on_cu_diff/ACF.dat                     # bader 输出：每个原子的 CHARGE
grep -A20 'Lowdin Charges' co_on_cu_diff/ab.projwfc.out | head -25
```

成功标准：

1. 所有 pw.x/pp.x 输出末尾有 `JOB DONE`；`si_rho.cube`、`h2o_elf.cube`、`drho.cube`、`ab_ae.cube`、`ab_stm.cube` 都能在 VESTA 或 XCrySDen 中打开；
2. ELF 取值在 0–1 之间，H2O 中 O 的两个孤对电子区域和 O–H 键处出现高 ELF；
3. `drho.cube` 的体积分接近 0（三者电子数守恒：N(AB) = N(slab) + N(CO)），等值面上 C 端与表面 Cu 之间出现电荷重新分布（正负两种颜色）；
4. `ACF.dat` 中 `CHARGE` 列的总和与 `ab.scf.out` 中的 `number of electrons` 一致（使用 `-ref` 全电子参考时，`CHARGE` 列积分的是价电子密度）；
5. Löwdin 输出中 `Spilling Parameter` 较小（几个百分点以内说明原子轨道基组能很好地展开 KS 态）。

本仓库不提供真实计算结果，所有数值以你自己的计算为准。

## 注意事项

- **差分电荷的三个计算必须共享同一个 FFT 网格**：同晶胞、同 `ecutrho`，必要时在三个输入中显式写相同的 `nr1, nr2, nr3`。
- `diff_plot.pp.in` 中第一个文件必须是原子数最多的那个（pp.x 从第一个文件读取原子坐标）。
- `co_on_cu_diff/ab.scf.in` 的几何是**未弛豫的初猜**：先按 [第5章](../../tutorials/05-relax-eos.md) 和 [第9章](../../tutorials/09-vdw-surface-workfunction.md) 弛豫，再把弛豫后的坐标同步到三个输入。
- `plot_num=21` 只适用于 PAW，并且需要很密的实空间网格；网格不够时 Bader 电荷会随网格变化，正式结果请做网格收敛（提高 `ecutrho` 或显式加密 `nr1/nr2/nr3`）。
- `plot_num=8`（ELF）：pp.x 每次都会打印 `elf/kin + US not fully implemented`（源码 `PP/src/elf.f90` 无条件输出该提示）。对 US/PAW，动能密度只含平滑部分，ELF 不完整；对 NC 赝势则没有这个问题，所以 ELF 示例使用 NC 赝势。
- 本例的 k 网格、`degauss`、截断能仅为演示，需按 [第4章](../../tutorials/04-convergence.md) 测试。
- 所有输入为 ASCII、LF 换行。
