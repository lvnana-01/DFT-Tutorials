# 示例 02：pw.x 输入结构与第一个 SCF（Si、Al、H2O）

配套教程：[第2章 pw.x 输入文件结构与第一个 SCF](../../tutorials/02-input-anatomy.md)

## 目的

- 用三个最小算例走通“写输入 → 运行 pw.x → grep 输出”的流程：
  - 半导体：金刚石 Si（`occupations` 默认为 fixed，输出 HOMO/LUMO）；
  - 金属：fcc Al（`occupations='smearing'`，输出费米能和 `-TS` 项）；
  - 分子：12 Å 立方盒子中的 H2O（`K_POINTS gamma` + `assume_isolated='mt'`）。
- 对比 `ibrav=2` + `celldm(1)` 与 `ibrav=0` + `CELL_PARAMETERS angstrom` 两种等价写法。

> [!NOTE]
> 本仓库没有运行过这些输入，下面的“预期结果”只描述应该检查什么，不提供数值。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.scf.in` | Si，`ibrav=2`，`celldm(1)=10.26` Bohr，30/240 Ry，8×8×8 k 网格，打印力和应力 |
| `si.scf.ibrav0.in` | 同一个 Si 晶胞，用 `ibrav=0` + `CELL_PARAMETERS angstrom`（a = 5.43 Å）+ `ATOMIC_POSITIONS crystal` |
| `al.scf.in` | Al，`ibrav=2`，`celldm(1)=7.65` Bohr，Marzari–Vanderbilt 冷展宽，`degauss=0.02` Ry，12×12×12 k 网格 |
| `h2o.scf.in` | H2O，`ibrav=1`，`A=12.0` Å，60/480 Ry（H 与 O 推荐值取最大），`nbnd=8`，Γ 点，Martyna–Tuckerman |
| `run_all.sh` | 默认只检查赝势与换行符；`RUN=1` 时依次运行四个输入并 grep 关键行 |

## 赝势准备（不提交 UPF 文件）

全部来自 SSSP 1.3 PBE Efficiency（下载与校验见 [第3章](../../tutorials/03-pseudopotentials.md)）：

| 元素 | 文件名 | 类型 | SSSP 推荐截断 ecutwfc / ecutrho（Ry） |
|---|---|---|---|
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | US | 30 / 240 |
| Al | `Al.pbe-n-kjpaw_psl.1.0.0.UPF` | PAW | 30 / 240 |
| O | `O.pbe-n-kjpaw_psl.0.1.UPF` | PAW | 50 / 400 |
| H | `H.pbe-rrkjus_psl.1.0.0.UPF` | US | 60 / 480 |

```bash
cd QE/examples/02-input-anatomy
mkdir -p pseudo
# copy the four files from the extracted SSSP_1.3.0_PBE_efficiency archive, e.g.
# cp ~/pseudo/SSSP_1.3.0_PBE_efficiency/{Si.pbe-n-rrkjus_psl.1.0.0.UPF,Al.pbe-n-kjpaw_psl.1.0.0.UPF,O.pbe-n-kjpaw_psl.0.1.UPF,H.pbe-rrkjus_psl.1.0.0.UPF} pseudo/
```

## 运行步骤

```bash
cd QE/examples/02-input-anatomy

# 1) static checks (no calculation)
bash run_all.sh

# 2) run everything yourself when ready
RUN=1 NP=4 bash run_all.sh

# 3) or run one input by hand
mpirun -np 4 pw.x -nk 2 -i si.scf.in > si.scf.out
```

集群上请写成作业脚本提交，见仓库的 `QE/scripts/job_qe_slurm.sh` 与 [公共篇第2章](../../../Common/tutorials/02-job-scheduler.md)。

## 预期结果 / 如何检查

```bash
grep -A12 "bravais-lattice index" si.scf.out   # cell, electrons, cutoffs, XC as understood by pw.x
grep "^!" *.out                                # final total energies (Ry)
grep -e "highest occupied" -e "Fermi energy" *.out
grep "convergence has been achieved" *.out
grep -A3 "total   stress" si.scf.out al.scf.out
grep "JOB DONE" *.out
```

检查要点：

| 算例 | 应当看到 |
|---|---|
| `si.scf.in` | `number of electrons = 8.00`；只有 `highest occupied level`（默认 `nbnd`=4，没有空带）；两个原子的力都为 0 |
| `si.scf.ibrav0.in` | 对称操作数、k 点数和 `!` 总能量与 `si.scf.in` 一致（数值精度内） |
| `al.scf.in` | `the Fermi energy is ... ev`；`smearing contrib. (-TS)` 很小 |
| `h2o.scf.in` | `highest occupied, lowest unoccupied level (ev)` 两个数；力不为 0（几何结构未优化） |

（示意输出，数值仅供参考）

```text
si.scf.out:!    total energy              =     -xx.xxxxxxxx Ry
al.scf.out:!    total energy              =     -xx.xxxxxxxx Ry
h2o.scf.out:!    total energy              =     -xx.xxxxxxxx Ry
```

## 注意事项

- `celldm(1)` 的单位是 Bohr，`A` 的单位是 Å，两者不能同时出现。
- 四个输入共用 `outdir=./tmp`，但 `prefix` 各不相同，所以不会互相覆盖。
- H2O 的盒子大小需要收敛测试（例如 10、12、14、16 Å）；`ecutwfc`/`ecutrho` 取所有元素推荐值的最大值。
- 绝对总能量依赖于赝势，只比较同一套设置下的能量差。
