# 示例 15：Si 的最大局域 Wannier 函数与插值能带（pw.x + pw2wannier90.x + Wannier90）

对应教程：[第15章 Wannier90 接口与插值能带](../../tutorials/15-wannier90.md)

## 目的

- 走通 QE + Wannier90 的标准五步流程：`pw.x` scf → `pw.x` nscf（全网格显式 k 列表）→ `wannier90.x -pp` → `pw2wannier90.x` → `wannier90.x`；
- 价带孤立流形：4 个以 Si–Si 键中心为初始投影的 MLWF（不需要解纠缠）；
- 价带 + 低导带：8 个 sp3 投影，使用外窗口/冻结窗口做解纠缠（disentanglement）；
- 输出插值能带 `*_band.dat`、Wannier 中心 `*_centres.xyz`、实空间哈密顿量 `si_hr.dat`。

## 文件说明

| 文件 | 说明 |
|---|---|
| `si.scf.in` | Si 自洽，PBE，SSSP 1.3 Efficiency 超软赝势，8×8×8 网格 |
| `si.nscf.in` | 非自洽：`K_POINTS crystal` 显式列出 4×4×4 全网格（64 点），`nosym`、`noinv`，`nbnd = 12` |
| `si_val.win` | 4 个 MLWF（价带），`exclude_bands = 5-12`，键中心 s 轨道投影 |
| `si.win` | 8 个 MLWF（价带 + 导带），`Si:sp3` 投影，`dis_win_max`、`dis_froz_max`，`write_hr` |
| `si_val.pw2wan.in`, `si.pw2wan.in` | `pw2wannier90.x` 输入（`seedname` 分别为 `si_val`、`si`） |
| `job.slurm` | 依次运行五步；`SEEDS` 变量选择 `si_val`、`si` |

## 赝势（不随仓库提供）

从 Materials Cloud 下载 SSSP 1.3 PBE Efficiency 压缩包，取出 Si 的文件：

- 文件名：`Si.pbe-n-rrkjus_psl.1.0.0.UPF`（pslibrary 1.0.0 超软赝势）
- SSSP 推荐截断：`ecutwfc = 30` Ry，`ecutrho = 240` Ry

```bash
mkdir -p pseudo
# download the SSSP 1.3 PBE efficiency archive from
#   https://www.materialscloud.org/discover/sssp/table/efficiency
# then extract only the Si file into ./pseudo
ls pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF
```

`pw2wannier90.x` 支持超软赝势与 PAW；也可以换成任意 PBE 模守恒赝势（相应修改 `ATOMIC_SPECIES` 与截断）。

## 运行步骤

```bash
cd QE/examples/15-wannier90
sbatch job.slurm                                 # both seeds
sbatch --export=ALL,SEEDS="si_val" job.slurm     # valence only
```

手动运行：

```bash
mpirun -np 8 pw.x -nk 4 -in si.scf.in  > si.scf.out
mpirun -np 8 pw.x -nk 4 -in si.nscf.in > si.nscf.out
wannier90.x -pp si_val                                   # writes si_val.nnkp
mpirun -np 8 pw2wannier90.x -in si_val.pw2wan.in > si_val.pw2wan.out
wannier90.x si_val                                       # writes si_val.wout
```

`wannier90.x` 是串行程序（另有 MPI 版本时可并行）；`pw2wannier90.x` 可以 MPI 并行，但不支持 `-nb`。

## 预期结果 / 如何检查

```bash
grep "highest occupied" si.scf.out                     # VBM: place dis_froz_max above it
grep -c "k =" si.nscf.out                              # k blocks printed with verbosity='high'
grep -e "K-points are ok" -e "Reciprocal lattice is ok" si_val.pw2wan.out
grep -A6 "Final State" si_val.wout
grep "<-- CONV" si_val.wout | tail -3
grep "All done" si_val.wout
ls si_val_band.dat si_val_band.gnu si_val_centres.xyz
```

- `si_val`：4 个 Wannier 中心应位于 4 个 Si–Si 键中心附近，4 个 spread 相等（对称性要求）；收敛过程中总 spread 单调下降并稳定。
- `si`：8 个 sp3 型 WF 中心位于两个 Si 原子附近；冻结窗口内的插值能带应与 DFT 本征值重合。
- 插值能带的价带应与 [第6章](../../tutorials/06-band-dos-pdos.md) 用 `calculation='bands'` 得到的 PBE 能带一致；在窗口以外（高导带）出现偏差是正常的。

## 注意事项

- nscf 的 k 列表必须与 `.win` 中 `kpoints` 块**完全相同且顺序一致**，否则 `pw2wannier90.x` 报 `Wrong number of k-points` 或 `problems with k-points`。本例两者都由同一程序生成。
- `dis_froz_max`、`dis_win_max` 是绝对能量（eV），不随费米能/价带顶平移；更换赝势、晶格常数后要根据 `si.scf.out` 重新设置。冻结窗口内的能带数不能超过 `num_wann`，否则报 `More states in the frozen window than target WFs`。
- `exclude_bands` 写在 `.win` 中，`wannier90.x -pp` 会把它写进 `.nnkp`，`pw2wannier90.x` 自动读取。
- k 路径坐标是相对于本例 `unit_cell_cart` 那组晶格矢量的倒格子分数坐标；换了晶胞矢量就要重新换算。
- 本目录所有数值都需要你自己算；教程中的输出片段只是格式示意。
