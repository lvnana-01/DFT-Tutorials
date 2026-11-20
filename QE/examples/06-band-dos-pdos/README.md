# 示例 06：Si 的能带、态密度与投影态密度

配套教程：[第6章 能带结构、态密度与投影态密度](../../tutorials/06-band-dos-pdos.md)

## 目的

用金刚石 Si 走完 QE 电子结构后处理的标准流程：

```text
pw.x scf -> pw.x bands (crystal_b path) -> bands.x
         -> pw.x nscf (dense grid, tetrahedra_opt) -> dos.x
                                                   -> projwfc.x -> sumpdos.x
```

## 文件说明

| 文件 | 程序 | 说明 |
|---|---|---|
| `si.scf.in` | pw.x | 自洽计算，8×8×8 网格，`nbnd=8` 以便输出 HOMO/LUMO |
| `si.bands.in` | pw.x | `calculation='bands'`，`K_POINTS crystal_b`：L–Γ–X–U\|K–Γ，`nbnd=12` |
| `si.bands_pp.in` | bands.x | `&BANDS`：`filband='si.bands.dat'`，`lsym=.true.` |
| `si.nscf.in` | pw.x | `calculation='nscf'`，16×16×16，`occupations='tetrahedra_opt'` |
| `si.dos.in` | dos.x | `&DOS`：`fildos='si.dos'`，`Emin/Emax/DeltaE`（eV） |
| `si.projwfc.in` | projwfc.x | `&PROJWFC`：`filpdos='si'`，输出 `si.pdos_tot` 与 `si.pdos_atm#N(Si)_wfc#M(l)` |
| `run_all.sh` | bash | 按顺序运行全部步骤并 grep 关键行 |

## 赝势准备（不提交 UPF 文件）

从 [SSSP 1.3 PBE Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency) 下载 `Si.pbe-n-rrkjus_psl.1.0.0.UPF`（SSSP 提示截断 30/240 Ry），放到 `pseudo/`。获取方法见 [第3章](../../tutorials/03-pseudopotentials.md)。

## 运行步骤

```bash
cd QE/examples/06-band-dos-pdos
mkdir -p pseudo          # 放入 Si.pbe-n-rrkjus_psl.1.0.0.UPF

# 一次跑完
bash run_all.sh
# 或指定并行方式
PW_CMD="mpirun -np 8 pw.x -nk 4" PP_CMD="mpirun -np 8" bash run_all.sh

# 分步运行（与 run_all.sh 等价）
mpirun -np 4 pw.x -in si.scf.in   > si.scf.out
mpirun -np 4 pw.x -in si.bands.in > si.bands.out
mpirun -np 4 bands.x -in si.bands_pp.in > si.bands_pp.out
mpirun -np 4 pw.x -in si.nscf.in  > si.nscf.out
mpirun -np 4 dos.x -in si.dos.in > si.dos.out
mpirun -np 4 projwfc.x -in si.projwfc.in > si.projwfc.out
```

作图：仓库脚本 `QE/scripts/qe_plot_bands.py`（读 `si.bands.dat.gnu`）和 `qe_plot_dos.py`（读 `si.dos` 与 `si.pdos_*`），分别用 `python ../../scripts/qe_plot_bands.py -h`、`python ../../scripts/qe_plot_dos.py -h` 查看用法。

## 预期结果 / 如何检查

| 检查项 | 命令 | 成功标准 |
|---|---|---|
| scf 正常结束 | `grep -e '^!' -e 'JOB DONE' si.scf.out` | 有总能与 `JOB DONE` |
| 带边 | `grep 'highest occupied, lowest unoccupied level' si.scf.out` | 两个数值（eV）；它们是 **8×8×8 网格上** 的带边，不一定是真实 CBM |
| 能带路径 | `grep 'high-symmetry point' si.bands_pp.out` | 列出 6 个高对称点及其 x 坐标 |
| 能带文件 | `ls si.bands.dat si.bands.dat.gnu si.bands.dat.rap` | 三个文件都存在 |
| DOS | `head -1 si.dos` | 表头含 `EFermi = ... eV` |
| PDOS | `ls si.pdos_*` | `si.pdos_tot`、`si.pdos_atm#1(Si)_wfc#1(s)`、`si.pdos_atm#1(Si)_wfc#2(p)` 等 |
| 投影质量 | `grep 'Spilling Parameter' si.projwfc.out` | 数值较小（越小说明原子轨道基越完备） |

定性预期：Si 为间接带隙半导体，价带顶在 Γ，导带底在 Γ–X 线上靠近 X 处；PBE 带隙明显低估实验值。具体数值以你的计算为准，本仓库不提供真实计算结果。

## 注意事项

- 四个 pw.x 输入的 `prefix`、`outdir`、`pseudo_dir`、`ecutwfc`、`ecutrho`、结构必须一致；bands/nscf 会读取 scf 写下的电荷密度。
- 注意 **nscf 会覆盖** `tmp/si.save` 中的波函数与能级：先跑 bands + bands.x，再跑 nscf + dos.x + projwfc.x（`run_all.sh` 就是这个顺序）。若要重画能带，重新运行 bands 步骤即可。
- `K_POINTS crystal_b` 中各点坐标基于 QE `ibrav=2` 的倒格矢；换用其它晶胞设置时需要重新生成路径（SeeK-path）。
- `run_all.sh` 中 `sumpdos.x` 的通配符需要用反斜杠转义括号。
- 所有输入为 ASCII、LF 换行。
