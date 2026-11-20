# 示例 14：Si/MgO 的 HSE、PBE0 与 r2SCAN，以及 HSE 能带的 Wannier 插值

对应教程：[第14章 杂化泛函与 meta-GGA](../../tutorials/14-hybrid-metagga.md)

## 目的

- 在相同截断与 k 网格下比较 Si 的 PBE 与 HSE 带隙（`input_dft = 'hse'`，`nqx1/2/3`）；
- 对 MgO 做 PBE0 自洽，熟悉外层 EXX 循环的输出（`est. exchange err (dexx)`、`EXX self-consistency reached`）；
- 在启用 Libxc 的 pw.x 上运行 r2SCAN（`input_dft = 'r2scan'`）；
- 用显式全网格 HSE 自洽 + Wannier90 得到 HSE 能带（杂化泛函在 QE 中画能带的推荐路线）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.pbe.scf.in` | Si，PBE，SG15 ONCV，4×4×4 网格，`nbnd = 8`（输出 HOMO/LUMO） |
| `si.hse.scf.in` | 同上，`input_dft='hse'`，`nqx1=nqx2=nqx3=2` |
| `mgo.pbe0.scf.in` | MgO 岩盐结构，`input_dft='pbe0'`，`nqx = 2` |
| `si.r2scan.scf.in` | Si，`input_dft='r2scan'`，需要 Libxc |
| `si.hse.fullgrid.in` | HSE 自洽，`K_POINTS crystal` 显式列出 4×4×4 全网格（64 点，等权重），`nosym`/`noinv`，`nbnd = 12` |
| `si_hse.win` | Wannier90 输入：8 个 sp3 WF，解纠缠窗口，`bands_plot`，k 路径 L–Γ–X–U/K–Γ |
| `si_hse.pw2wan.in` | `pw2wannier90.x` 输入：`write_mmn`、`write_amn` |
| `job.slurm` | 用 `STEPS` 选择 `pbe hse pbe0 r2scan wannier` |

## 赝势（不随仓库提供）

使用 SG15 ONCV 模守恒赝势（PBE，**无**非线性芯修正，避免 pw.x 对杂化泛函和 meta-GGA 给出 NLCC 不一致的警告）：

```bash
mkdir -p pseudo && cd pseudo
wget http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2020-02-06.tar.gz
tar -xzf sg15_oncv_upf_2020-02-06.tar.gz Si_ONCV_PBE-1.2.upf Mg_ONCV_PBE-1.2.upf O_ONCV_PBE-1.2.upf
grep -m1 'core_correction=' *.upf      # all "F"
cd ..
```

注意：

- pw.x 中杂化泛函配 US/PAW 赝势时不能做 relax/md（`forces for hybrid functionals + US/PAW not implemented`），`ecutfock` 也被忽略；meta-GGA 完全不支持 US/PAW（`Meta-GGA not implemented with USPP/PAW`）。所以本章统一用模守恒赝势。
- 截断值（Si 40 Ry、MgO 80 Ry、r2SCAN 50 Ry）只是起点，请按 [第4章](../../tutorials/04-convergence.md) 的方法自行收敛。

## 运行步骤

```bash
cd QE/examples/14-hybrid-metagga
sbatch --export=ALL,STEPS="pbe hse" job.slurm      # Si PBE vs HSE
sbatch --export=ALL,STEPS="pbe0" job.slurm         # MgO PBE0
sbatch --export=ALL,STEPS="r2scan" job.slurm       # only with a Libxc-enabled pw.x
sbatch --export=ALL,STEPS="wannier" job.slurm      # HSE bands via Wannier90
```

Wannier 部分的手动命令：

```bash
mpirun -np 16 pw.x -nk 4 -in si.hse.fullgrid.in > si.hse.fullgrid.out
wannier90.x -pp si_hse                              # writes si_hse.nnkp
mpirun -np 16 pw2wannier90.x -in si_hse.pw2wan.in > si_hse.pw2wan.out
wannier90.x si_hse                                  # writes si_hse.wout, si_hse_band.dat
```

## 预期结果 / 如何检查

```bash
grep "highest occupied" si.pbe.scf.out
grep "highest occupied" si.hse.scf.out | tail -1          # last occurrence = converged HSE
grep -e "est. exchange err" -e "EXX self-consistency reached" si.hse.scf.out
grep "EXX: q-point mesh" si.hse.scf.out                   # should echo 2 2 2
grep -A6 "Final State" si_hse.wout
```

- HSE 的 `highest occupied, lowest unoccupied level` 差值应明显大于 PBE（Si 的实验带隙约 1.17 eV，PBE 典型地低估一半左右）。注意：这里给出的是 **4×4×4 网格点上**的最小间隙，Si 的导带底不在网格点上，真正的带隙要看插值能带。
- `est. exchange err (dexx)` 应逐次下降，最后出现 `EXX self-consistency reached`，最终能量行以 `!!` 开头。
- Wannier90：`Final State` 中总 spread 稳定，`si_hse_band.dat` 中的价带与 `si.hse.fullgrid.out` 中网格点上的本征值一致（在冻结窗口内应几乎重合）。

## 注意事项

- 显式 k 列表（`K_POINTS crystal`）时没有默认的 EXX q 网格，必须写 `nqx1/2/3`；q 网格必须与 k 网格相容（k+q 必须在 k 列表中），否则 pw.x 会打印 `EXX: WARNING: q-point mesh has been updated!` 并自动改小 q 网格。
- `si_hse.win` 中的 `dis_froz_max` 是绝对能量（eV），要根据你自己的 `si.hse.fullgrid.out` 中价带顶的位置调整（冻结窗口放在价带顶以上约 1 eV 以内）。
- `pw2wannier90.x` 不支持 `-nb`（带组并行）。
- QE 7.6.1（开发分支，截至 2026-11 尚未发布）的 release notes 写明：此前版本杂化泛函的应力有误；从 7.6.1 起 US/PAW 下不再计算杂化泛函应力，vc-relax 会停止。在 7.6 及以前不要用杂化泛函做 vc-relax。
- 所有带隙数值都需要你自己算；教程中的输出片段只是格式示意。
