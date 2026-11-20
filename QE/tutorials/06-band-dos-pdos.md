# 第6章 能带结构、态密度与投影态密度（Band Structure, DOS & PDOS: bands.x, dos.x, projwfc.x）

本章讲 Quantum ESPRESSO 的电子结构后处理：先用 pw.x 做 scf，再用 `calculation='bands'` 加 `K_POINTS crystal_b` 计算高对称路径上的能带，交给 bands.x（`&BANDS`）整理；然后用 `calculation='nscf'` 在更密的均匀网格上配合 `occupations='tetrahedra_opt'` 计算本征值，用 dos.x（`&DOS`）得到总态密度（DOS），用 projwfc.x（`&PROJWFC`）得到投影态密度（PDOS）、Löwdin 电荷和 spilling 参数，并用 sumpdos.x 合并。我们以间接带隙半导体 Si 为主线，讨论 SeeK-path 路径、`nbnd`、费米能与带隙的正确读取，并对比金属 Cu 与 GaAs 的注意事项。

> [!NOTE]
> 本章全部输入在示例目录 [06-band-dos-pdos](../examples/06-band-dos-pdos/README.md)。截断能与 k 网格来自 [第4章](04-convergence.md)，晶格常数应先按 [第5章](05-relax-eos.md) 优化。文中"示意输出"只演示格式，**数值仅供参考，不是真实计算结果**。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 工作流 | scf → bands → bands.x；scf → nscf → dos.x / projwfc.x，理解每一步读写什么 |
| k 路径 | 用 SeeK-path 或教材表格生成高对称点；分清 `crystal_b` 与 `tpiba_b`；理解"权重 = 段内点数" |
| 空带 | 何时以及如何设置 `nbnd`，`diago_full_acc` 的作用 |
| DOS | `tetrahedra_opt` vs smearing；`Emin`/`Emax`/`DeltaE`；DOS 文件格式 |
| PDOS | `pdos_tot` 与 `pdos_atm#N(X)_wfc#M(l)` 文件的列含义；m 分量顺序；sumpdos.x 合并 |
| 能量零点 | 费米能、`highest occupied, lowest unoccupied level` 应以 scf/nscf 为准 |
| 作图 | gnuplot、Python（`qe_plot_bands.py`、`qe_plot_dos.py`） |

与 VASP 用户的对应关系：

| VASP | QE | 说明 |
|---|---|---|
| `ICHARG=11` + line-mode KPOINTS | `calculation='bands'` + `K_POINTS crystal_b` | 都是读自洽电荷的非自洽计算 |
| `ICHARG=11` + 均匀密网格 | `calculation='nscf'` + `K_POINTS automatic` | DOS 用 |
| `NBANDS` | `nbnd` | QE 绝缘体默认只算占据带 |
| `ISMEAR=-5` | `occupations='tetrahedra_opt'`（或 `'tetrahedra'`） | 平滑 DOS |
| `NEDOS`、`EMIN`/`EMAX` | dos.x 的 `DeltaE`、`Emin`/`Emax` | QE 在后处理里设定，单位 eV |
| `LORBIT=11`（PROCAR/DOSCAR 分波） | projwfc.x | QE 投影到正交化原子轨道（Löwdin） |
| `E-fermi`（OUTCAR） | `the Fermi energy is` / `highest occupied, lowest unoccupied level` | 以 scf/nscf 为准 |
| vaspkit / sumo / py4vasp 作图 | bands.x 的 `.gnu`、plotband.x、Python | — |

VASP 的同类内容见 [VASP 第5章](../../VASP/tutorials/05-band-dos.md)，CP2K 的见 [CP2K 第10章](../../CP2K/tutorials/10-dos-pdos-band.md)。

---

## 2. 原理速览

### 2.1 为什么要分两步

Kohn-Sham 方程的自洽只需要在一个均匀 k 网格上求得电荷密度 $n(\mathbf r)$。一旦 $n(\mathbf r)$ 收敛，哈密顿量就确定了，任何 k 点上的本征值都可以通过**非自洽**对角化得到：

- `calculation='bands'`：在任意 k 点列表（通常是高对称路径）上求本征值，**不计算费米能**，也不要求 k 点权重有意义；
- `calculation='nscf'`：在均匀网格上求本征值，并计算占据数与费米能，供 DOS/PDOS 使用。

两者都读 scf 写入 `outdir/prefix.save/` 的电荷密度，因此 `prefix`、`outdir`、赝势、截断能、结构**必须与 scf 完全一致**。

### 2.2 态密度与四面体方法

$$
g(E) = \sum_n \int_\text{BZ} \frac{d^3k}{\Omega_\text{BZ}}\, \delta(E-\varepsilon_{n\mathbf k})
$$

δ 函数需要离散化处理：

- **展宽（smearing）**：把 δ 换成宽度为 `degauss` 的高斯等函数，简单但峰形被抹平，k 点稀时出现锯齿或假峰；
- **四面体方法**：在 k 网格的四面体内线性插值能带再解析积分。`'tetrahedra'` 是 Blöchl 版本；`'tetrahedra_opt'` 是 Kawamura 的优化四面体方法，INPUT_PW 中注明它比线性四面体更高效（`'tetrahedra_lin'` 仅作参考）。四面体方法需要 pw.x 自动生成的均匀网格（`K_POINTS automatic`）。

官方文档同时提醒：四面体方法**不是变分的**，不适合力/优化/动力学；所以只在 nscf（DOS 专用）中使用，scf 仍按第4章的设置。

### 2.3 投影态密度（PDOS）

projwfc.x 把每个 Kohn-Sham 态投影到赝势文件中给出的**正交化原子波函数**（Löwdin 正交化）上，得到每个原子、每个角动量 $l$、每个 $m$ 分量的权重，进而计算 PDOS 与 Löwdin 电荷。原子轨道基不完备，投影权重之和小于 1，差值就是 **spilling parameter**（`Spilling Parameter:` 一行）。因此：

- PDOS 之和一般小于总 DOS；
- 能用于判断"哪个原子、哪个轨道主导价带顶/导带底"，但不要把 Löwdin 电荷当作精确的离子价态；
- 投影基来自 UPF 中的 `PP_CHI`；若赝势没有某个轨道（如缺少 d 态投影），PDOS 中就不会出现它。

PAW 赝势还可以设 `pawproj=.true.`，用 PAW 投影子与全电子分波计算权重（非共线情况未实现）。

---

## 3. 流程总览

```text
            si.scf.in                         si.bands.in              si.bands_pp.in
pw.x  -----------------> tmp/si.save  ---pw.x (bands)-----> eigenvalues ---bands.x---> si.bands.dat(.gnu/.rap)
                              |
                              +------pw.x (nscf, si.nscf.in)--> eigenvalues on 16x16x16
                                                                   |-- dos.x (si.dos.in) ------> si.dos
                                                                   |-- projwfc.x (si.projwfc.in) -> si.pdos_*
```

> [!WARNING]
> bands 和 nscf 都会把新的本征值写回同一个 `tmp/si.save`。所以顺序应为 **bands → bands.x**，然后 **nscf → dos.x → projwfc.x**。若在 nscf 之后再运行 bands.x，读到的就是 nscf 网格上的数据。需要两者并存时，可以把 scf 的 `tmp` 目录复制一份分别使用。

---

## 4. 第一步：scf

[`si.scf.in`](../examples/06-band-dos-pdos/si.scf.in)：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.26           ! replace by your relaxed value (Chapter 5)
  nat       = 2
  ntyp      = 1
  ecutwfc   = 40.0
  ecutrho   = 320.0
  occupations = 'fixed'
  nbnd      = 8               ! 4 occupied + 4 empty: also prints the LUMO
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.7
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

`nbnd` 的默认值：绝缘体 = 价电子数/2（Si 两个原子共 8 个价电子 → 4 条带），金属多 20%（至少多 4 条）。scf 中多算 4 条空带的好处是输出会同时给出 HOMO 和 LUMO：

```bash
mpirun -np 4 pw.x -in si.scf.in > si.scf.out
grep -e '^!' -e 'highest occupied' -e 'JOB DONE' si.scf.out
```

（示意输出，数值仅供参考）

```text
     highest occupied, lowest unoccupied level (ev):     x.xxxx    x.xxxx
!    total energy              =     -15.8xxxxxxx Ry
```

注意：这里的 LUMO 是 **8×8×8 网格上** 的最低空能级。Si 的导带底不在高对称点上（在 Γ–X 线上），网格未必恰好包含它，所以这个"带隙"只是估计值；精确的带边要结合能带路径或更密网格。

---

## 5. 第二步：能带（bands + bands.x）

### 5.1 高对称路径从哪里来

推荐顺序：

1. **SeeK-path**（[在线工具](https://www.materialscloud.org/work/tools/seekpath)，或 Python 包 `seekpath`）：上传结构，可直接导出 QE 的 `K_POINTS crystal_b` 块。注意 SeeK-path 可能会把结构转换成**它的标准原胞**，此时 scf 也要改用它给出的晶胞，否则坐标对不上；
2. **教材/文献表格**（Setyawan-Curtarolo 2010、Hinuma 2017）；
3. 仓库脚本 `QE/scripts/qe_inputgen.py` 可以结合 seekpath 生成路径，在示例目录下用 `python ../../scripts/qe_inputgen.py -h` 查看用法。

### 5.2 `crystal_b` 与 `tpiba_b`

| 格式 | 坐标单位 | 适用 |
|---|---|---|
| `tpiba_b` | 笛卡尔坐标，单位 $2\pi/a$ | 教科书上的立方晶系高对称点（如 X = (0,1,0)）可直接抄 |
| `crystal_b` | 倒格矢 $\mathbf b_1,\mathbf b_2,\mathbf b_3$ 的分数坐标 | 任意晶系；SeeK-path 输出常用 |

两者语法相同：第一行是点数 `nks`，之后每行 `kx ky kz N`。pw.x 在每一对相邻点之间生成 N 个点，**N 是起点那一行的权重**；最后一个点的权重不使用（写 1 即可）。若某点权重为 1，则它与下一点之间只有一个点，相当于"跳"到下一点——这就是在路径中写 `U|K` 不连续段的方法。

> [!IMPORTANT]
> 同一个点的 `crystal_b` 坐标依赖于晶胞矢量的约定。QE `ibrav=2` 的 fcc 原胞矢量是 $\mathbf a_1=\frac a2(-1,0,1)$、$\mathbf a_2=\frac a2(0,1,1)$、$\mathbf a_3=\frac a2(-1,1,0)$，与 VASP/Materials Project 常用的约定不同，所以别处抄来的 `crystal` 坐标不能直接用。下表的坐标是针对 QE `ibrav=2` 推导的。

| 点 | 笛卡尔（$2\pi/a$，`tpiba_b`） | QE `ibrav=2` 的 `crystal_b` |
|---|---|---|
| Γ | (0, 0, 0) | (0, 0, 0) |
| X | (0, 1, 0) | (0, 0.5, 0.5) |
| L | (0.5, 0.5, 0.5) | (0, 0.5, 0) |
| K | (0.75, 0.75, 0) | (−0.375, 0.375, 0) |
| U | (0.25, 1, 0.25) | (0, 0.625, 0.375) |
| W | (0.5, 1, 0) | (−0.25, 0.5, 0.25) |

（坐标由 $\mathbf k = \sum_i c_i \mathbf b_i$ 换算得到；如有疑问，在 bands 输出中用 `verbosity='high'` 打印的 k 点笛卡尔坐标核对。）

### 5.3 bands 输入

[`si.bands.in`](../examples/06-band-dos-pdos/si.bands.in)：

```fortran
&CONTROL
  calculation = 'bands'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  verbosity   = 'high'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.26
  nat       = 2
  ntyp      = 1
  ecutwfc   = 40.0
  ecutrho   = 320.0
  occupations = 'fixed'
  nbnd      = 12              ! more empty bands to see the conduction band
/
&ELECTRONS
  conv_thr       = 1.0d-10
  diago_full_acc = .true.     ! empty states as accurate as occupied ones
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS crystal_b
6
  0.000  0.500  0.000  35   ! L
  0.000  0.000  0.000  40   ! G
  0.000  0.500  0.500  14   ! X
  0.000  0.625  0.375   1   ! U
 -0.375  0.375  0.000  42   ! K
  0.000  0.000  0.000   1   ! G
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `calculation` | `'bands'` | 非自洽，读 scf 电荷；不计算费米能 |
| `nbnd` | 12 | 默认只有 4 条占据带，必须显式增加才能看到导带 |
| `diago_full_acc` | `.true.` | 默认空带用较松阈值对角化；画导带时建议打开 |
| `conv_thr` | 1e-10 | nscf/bands 中用于设定对角化阈值 |
| `verbosity` | `'high'` | 打印每个 k 点的坐标与本征值，便于核对 |
| 权重 35/40/14/42 | 段内点数 | 大致按段长分配，使点密度均匀 |

运行：

```bash
mpirun -np 4 pw.x -in si.bands.in > si.bands.out
grep -e 'Band Structure Calculation' -e 'End of band structure calculation' si.bands.out
grep -c '          k =' si.bands.out          # 粗略统计输出的 k 点块数
```

### 5.4 bands.x

[`si.bands_pp.in`](../examples/06-band-dos-pdos/si.bands_pp.in)：

```fortran
&BANDS
  prefix  = 'si'
  outdir  = './tmp'
  filband = 'si.bands.dat'
  lsym    = .true.            ! symmetry analysis -> si.bands.dat.rap
/
```

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `prefix` / `outdir` | `'pwscf'` / `ESPRESSO_TMPDIR` 或 `./` | 与 pw.x 一致 |
| `filband` | `'bands.out'` | 输出文件名；同时生成 `filband.gnu`（eV，可直接用 gnuplot 画）和（`lsym` 时）`filband.rap` |
| `lsym` | `.true.` | 按 k 点小群不可约表示分类能带，写入 `.rap` 文件 |
| `no_overlap` | `.true.` | 为 `.false.` 且 `lsym=.false.` 时，按相邻 k 点波函数重叠最大来排序能带 |
| `spin_component` | — | 自旋极化时选 1（up）或 2（down），见 [第7章](07-magnetism-soc.md) |
| `lsigma(i)` | — | 非共线时计算自旋算符期望值 |
| `lp` / `filp` | `.false.` / `'p_avg.dat'` | 计算导带–价带间动量矩阵元 |

```bash
mpirun -np 4 bands.x -in si.bands_pp.in > si.bands_pp.out
grep 'high-symmetry point' si.bands_pp.out
```

（示意输出，数值仅供参考）

```text
     high-symmetry point:  0.5000 0.5000 0.5000   x coordinate   0.0000
     high-symmetry point:  0.0000 0.0000 0.0000   x coordinate   0.8660
     high-symmetry point:  0.0000 1.0000 0.0000   x coordinate   1.8660
     ...
```

`high-symmetry point` 后的坐标是笛卡尔坐标（$2\pi/a$），`x coordinate` 是沿路径的累计长度，正是作图时高对称点竖线的位置。注意 bands.x 是靠路径方向变化来"猜"高对称点的，U|K 跳跃处的标签需要你自己按输入顺序标注。

### 5.5 作图

`si.bands.dat.gnu` 是两列（路径坐标、能量 eV）、能带之间以空行分隔的文件。最简单的 gnuplot：

```bash
gnuplot -persist << 'EOF'
set ylabel "E (eV)"
set xtics ("L" 0.0, "G" 0.8660, "X" 1.8660)   # complete with your x coordinates
plot "si.bands.dat.gnu" u 1:2 w l lc rgb "black" notitle
EOF
```

能量零点应取价带顶（绝缘体/半导体）或费米能（金属），数值从 **scf 或 nscf 输出**中读取（第 8 节）。仓库脚本 `QE/scripts/qe_plot_bands.py` 可以读取 `.gnu` 文件和 bands.x 输出，自动对齐能量零点、加标签并报告带隙，在示例目录下用 `python ../../scripts/qe_plot_bands.py -h` 查看用法。QE 自带的交互式 plotband.x 也能把 `filband` 转成 gnuplot/PostScript 格式。

---

## 6. 第三步：nscf + dos.x

### 6.1 nscf 输入

[`si.nscf.in`](../examples/06-band-dos-pdos/si.nscf.in) 与 scf 相比改了 `calculation`、`occupations`、`nbnd`、`K_POINTS`，并加了 `diago_full_acc`（其余必须保持一致）：

```fortran
&CONTROL
  calculation = 'nscf'
  ...
/
&SYSTEM
  ...
  occupations = 'tetrahedra_opt'
  nbnd      = 12
/
...
K_POINTS automatic
  16 16 16 0 0 0
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `calculation` | `'nscf'` | 均匀网格非自洽；会计算费米能 |
| `occupations` | `'tetrahedra_opt'` | 平滑 DOS；需要 `K_POINTS automatic`，否则报 `tetrahedra need automatic k-point grid` |
| `nbnd` | 12 | 覆盖你想看的导带能量范围 |
| `K_POINTS` | 16×16×16 | DOS 通常要比 scf 密 1.5–2 倍或更多；检查 DOS 不再随网格变化 |

```bash
mpirun -np 4 pw.x -in si.nscf.in > si.nscf.out
grep -e 'the Fermi energy is' -e 'highest occupied' si.nscf.out
grep 'number of k points' si.nscf.out
```

使用四面体或展宽时，pw.x 会打印 `the Fermi energy is ... ev`（对半导体它落在带隙中的某处）；`occupations='fixed'` 时才打印 `highest occupied, lowest unoccupied level`。

### 6.2 dos.x 输入

[`si.dos.in`](../examples/06-band-dos-pdos/si.dos.in)：

```fortran
&DOS
  prefix = 'si'
  outdir = './tmp'
  fildos = 'si.dos'
  Emin   = -15.0              ! eV (absolute energies, not shifted to E_F)
  Emax   =  15.0
  DeltaE =  0.01              ! eV grid step
/
```

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `fildos` | `'<prefix>.dos'` | 输出文件 |
| `Emin`, `Emax` | 能带极值（有展宽时外扩 3 倍展宽） | 单位 **eV**，绝对能量 |
| `DeltaE` | — | 能量网格步长（eV） |
| `bz_sum` | 若给 `degauss` 则 `'smearing'`，否则从 XML 读取 | 可选 `'smearing'`、`'tetrahedra'`、`'tetrahedra_lin'`、`'tetrahedra_opt'` |
| `degauss` | — | 单位 **Ry**（不是 eV）；**一旦给出就改用高斯展宽** |
| `ngauss` | 0 | 0 高斯，1 MP，−1 冷展宽，−99 Fermi-Dirac |

官方文档对展宽选择的规则：只有当 pw.x 数据使用了四面体占据、**并且** `&DOS` 中没有给 `degauss` 时才用四面体方法；其他情况都用高斯展宽（给了 `degauss` 就用它，否则从数据文件读 scf 的 `degauss`/`ngauss`，都没有则用 `degauss = DeltaE`）。

> [!NOTE]
> dos.x 从 QE 5 起 namelist 名为 `&DOS`，不再是 `&INPUTPP`。旧教程中的 `&INPUTPP` 写法会报错。

```bash
mpirun -np 4 dos.x -in si.dos.in > si.dos.out
head -3 si.dos
```

（示意输出，数值仅供参考）

```text
#  E (eV)   dos(E)     Int dos(E) EFermi =    x.xxx eV
 -15.000  0.0000E+00  0.0000E+00
 -14.990  0.0000E+00  0.0000E+00
```

三列分别是能量（eV）、DOS（states/eV）、积分 DOS；自旋极化时为 `E dosup(E) dosdw(E) Int dos(E)`。积分 DOS 在价带顶处应等于价电子数（Si 原胞为 8），这是一个很好的自检。

---

## 7. 第四步：projwfc.x 与 PDOS

### 7.1 输入

[`si.projwfc.in`](../examples/06-band-dos-pdos/si.projwfc.in)：

```fortran
&PROJWFC
  prefix  = 'si'
  outdir  = './tmp'
  filpdos = 'si'
  Emin    = -15.0
  Emax    =  15.0
  DeltaE  =  0.01
/
```

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `filpdos` | `prefix` 的值 | PDOS 文件名前缀 |
| `filproj` | 标准输出 | 投影系数另存的文件 |
| `Emin`/`Emax`/`DeltaE` | 能带极值 / — | eV |
| `degauss`/`ngauss` | 0.0 / 0 | Ry；展宽规则与 dos.x 相同，不给 `degauss` 且 nscf 为四面体时用四面体 |
| `lsym` | `.true.` | 对投影做对称化并计算 PDOS |
| `diag_basis` | `.false.` | `.true.` 时投影到使占据矩阵对角化的局域坐标系 |
| `pawproj` | `.false.` | PAW 专用，用 PAW 投影子计算权重 |
| `kresolveddos` | `.false.` | k 分辨 DOS（每个 k 点单独输出） |
| `tdosinboxes` | `.false.` | 在实空间盒子内积分的局域 DOS |

投影需要波函数：projwfc.x 必须在 nscf **之后**、同一 `outdir` 下运行（数据文件说明里写着 "wavefunctions are needed"）。

```bash
mpirun -np 4 projwfc.x -in si.projwfc.in > si.projwfc.out
ls si.pdos_*
grep -A6 'Lowdin Charges' si.projwfc.out
grep 'Spilling Parameter' si.projwfc.out
```

### 7.2 输出文件

| 文件 | 列（非自旋极化） | 说明 |
|---|---|---|
| `si.pdos_tot` | `E  DOS(E)  PDOS(E)` | 总 DOS 与所有投影之和；两列之差反映 spilling |
| `si.pdos_atm#1(Si)_wfc#1(s)` | `E  LDOS(E)  PDOS_1(E)` | 原子 1 的 s 轨道 |
| `si.pdos_atm#1(Si)_wfc#2(p)` | `E  LDOS(E)  PDOS_1..3(E)` | 原子 1 的 p 轨道，LDOS = 各 m 分量之和 |
| `si.pdos_atm#2(Si)_wfc#...` | 同上 | 原子 2 |

文件名规则：`<filpdos>.pdos_atm#N(X)_wfc#M(l)`，N 为原子序号，X 为元素符号（ATOMIC_SPECIES 标签），M 为该原子赝势中的原子波函数序号，l 为 s/p/d/f。自旋极化时每个量拆成 up/dw 两列；含全相对论赝势（SOC）时文件名变为 `..._wfc#M(l_j)`，投影到总角动量本征态（见 [第7章](07-magnetism-soc.md)）。

m 分量的顺序（官方文档 "Orbital Order"）：

| l | 第 1 列 | 第 2 列 | 第 3 列 | 第 4 列 | 第 5 列 |
|---|---|---|---|---|---|
| p（l=1） | pz | px | py | — | — |
| d（l=2） | dz2 | dzx | dzy | dx2-y2 | dxy |

> [!WARNING]
> p 轨道的顺序是 **pz, px, py**，不是 px, py, pz。从 VASP 迁移过来的用户特别容易弄错。

### 7.3 用 sumpdos.x 合并

sumpdos.x 把多个 PDOS 文件求和后写到**标准输出**（`sumpdos.x -h` 显示用法：`sumpdos [-h] [-f <filein>] [-Ef <value>] [<file1> ... <fileN>]`）：

```bash
# 所有 Si 原子的 s 和 p 分别求和（括号要转义）
sumpdos.x si.pdos_atm#*\(Si\)_wfc#*\(s\) > si.pdos_Si_s.dat
sumpdos.x si.pdos_atm#*\(Si\)_wfc#*\(p\) > si.pdos_Si_p.dat
# 把能量平移到费米能（eV）
sumpdos.x -Ef 6.xx si.pdos_atm#*\(Si\)_wfc#*\(p\) > si.pdos_Si_p_shift.dat
```

（`-Ef` 后的数值请替换为你自己 nscf 输出中的值。）

仓库脚本 `QE/scripts/qe_plot_dos.py` 可以直接读 `si.dos` 与 `si.pdos_*`，按元素/轨道求和、对齐能量零点并作图，在示例目录下用 `python ../../scripts/qe_plot_dos.py -h` 查看用法。

### 7.4 Löwdin 电荷

`si.projwfc.out` 中的 `Lowdin Charges:` 块给出每个原子的总电荷与 s/p/d 分解，末尾有 `Spilling Parameter:`。（示意输出，数值仅供参考）

```text
Lowdin Charges:

     Atom #   1: total charge =   x.xxxx, s =  x.xxxx, p =  x.xxxx,
     Atom #   2: total charge =   x.xxxx, s =  x.xxxx, p =  x.xxxx,
     Spilling Parameter:   0.0xxx
```

对 Si 这样的单质，两个原子的 Löwdin 电荷应相同；对 GaAs 等化合物，可用来定性比较电荷转移方向。定量电荷分析请用 Bader（[第10章](10-charge-elf-bader.md)）。

---

## 8. 费米能、带隙与能量零点

| 情形 | 去哪里读 | 输出行 |
|---|---|---|
| 绝缘体，`occupations='fixed'`，有空带 | scf 输出 | `highest occupied, lowest unoccupied level (ev):` |
| 绝缘体，`'fixed'`，无空带 | scf 输出 | `highest occupied level (ev):` |
| 金属或使用展宽/四面体 | scf 或 nscf 输出 | `the Fermi energy is ... ev` |
| DOS 文件 | `si.dos` 表头 | `EFermi = ... eV` |
| `calculation='bands'` 输出 | **不要用** | bands 运行只算本征值，不求占据、不打印费米能或 HOMO/LUMO；路径上的点也没有积分权重 |

求 Si 的带隙：

1. 价带顶：Si 在 Γ 点，从能带数据中取占据带最大值；
2. 导带底：在 Γ–X 线上，从能带数据中取空带最小值；
3. 带隙 = 二者之差；与 scf 网格给出的 `lowest unoccupied level` 比较，体会网格是否包含导带底。

> [!IMPORTANT]
> PBE 等半局域泛函系统性低估带隙。需要可信的带隙数值，请使用杂化泛函（[第14章](14-hybrid-metagga.md)）或 GW（Yambo）。本章的重点是流程与读取方法，而不是带隙的绝对值。

---

## 9. 其它体系的注意事项

### 9.1 金属 Cu（d 带）

- scf 用 `occupations='smearing'`（如 `'mv'`，`degauss` 按 [第4章](04-convergence.md) 测试），bands 计算中也保留同样的 smearing 设置；
- `nbnd` 默认值对金属已多出 20%（至少 4 条），画能带时常需再加；
- DOS 用 nscf + `tetrahedra_opt` + 密网格（金属比半导体更需要密网格）；
- 能量零点取 nscf 的 `the Fermi energy is`；
- Cu 3d 带位于费米能以下几 eV 处，可用 `pdos_atm#1(Cu)_wfc#2(d)`（序号以实际文件名为准）查看。SSSP 1.3 中 Cu 使用 `Cu.paw.z_11.ld1.psl.v1.0.0-low.upf`（提示截断 90/720 Ry）。

### 9.2 GaAs（直接带隙、半芯态）

- SSSP 1.3 Efficiency 中 Ga 为 `Ga.pbe-dn-kjpaw_psl.1.0.0.UPF`（包含 3d 半芯态，提示 70/560 Ry）、As 为 `As.pbe-n-rrkjus_psl.0.2.UPF`（35/280 Ry）；截断取两者最大值；
- Ga 3d 半芯带在价带底以下较深处，`Emin` 要设得足够低才能在 DOS 中看到；
- 闪锌矿与金刚石同为 fcc 晶格（`ibrav=2`），本章的 `crystal_b` 路径可以直接复用；
- 价带顶附近的 SOC 劈裂见 [第7章](07-magnetism-soc.md)。

### 9.3 自旋极化与 SOC

自旋极化时 bands.x 用 `spin_component = 1` 和 `2` 分别输出两个自旋通道，DOS/PDOS 文件变为 up/dw 两列；非共线/SOC 情况见 [第7章](07-magnetism-soc.md)。

---

## 参数选择与收敛建议

| 参数 | scf | bands | nscf（DOS/PDOS） |
|---|---|---|---|
| `K_POINTS` | 收敛测试值（Si：8×8×8） | `crystal_b`，每段 20–50 点，按段长分配 | 比 scf 更密（Si：16×16×16 起），直到 DOS 不变 |
| `occupations` | 绝缘体 `'fixed'`；金属 `'smearing'` | 与 scf 相同 | `'tetrahedra_opt'` |
| `nbnd` | 默认，或加几条以输出 LUMO | 占据带数 + 足够空带 | 覆盖 `Emax` |
| `diago_full_acc` | 默认 | `.true.` | `.true.` |
| `ecutwfc`/`ecutrho`、结构、赝势 | 收敛值 | **必须与 scf 相同** | **必须与 scf 相同** |
| dos.x `DeltaE` | — | — | 0.01–0.05 eV |
| dos.x/projwfc.x `degauss` | — | — | 不给（用四面体）；或 0.005–0.02 Ry 高斯 |

收敛检查：把 scf 网格加密一级，能带在路径上的变化应在 meV 量级；把 nscf 网格加密一级，DOS 主要峰位与峰形应基本不变。

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| 能带图只有价带，看不到导带 | `nbnd` 用了默认值（绝缘体只算占据带） | bands/nscf 中显式设 `nbnd` |
| `Error in routine setup: tetrahedra need automatic k-point grid` | 四面体占据配了显式 k 点列表 | nscf 使用 `K_POINTS automatic` |
| bands/nscf 报找不到数据文件，或结果与 scf 不符 | `prefix`/`outdir` 与 scf 不一致，或 scf 未成功 | 统一 `prefix`/`outdir`；确认 scf 有 `JOB DONE` |
| 能带形状怪异、高对称点位置不对 | `crystal_b` 坐标按别的晶胞约定写的；`crystal_b` 与 `tpiba_b` 混用 | 按本章表格或 SeeK-path 重新生成；用 `verbosity='high'` 核对笛卡尔坐标 |
| 用 SeeK-path 路径后能带完全不对 | SeeK-path 把结构转换成了它的标准原胞，但 scf 仍用旧晶胞 | scf/bands 都改用 SeeK-path 输出的晶胞 |
| 能带图中导带"毛刺"或交叉处乱跳 | 空带对角化阈值太松；能带排序问题 | `diago_full_acc=.true.`；bands.x 中试 `lsym=.false.`、`no_overlap=.false.` |
| DOS 锯齿严重、出现假峰 | nscf k 点太稀，或用了很小的高斯展宽 | 加密网格，用 `tetrahedra_opt` |
| 设了 `tetrahedra_opt` 但 DOS 仍是高斯展宽 | `&DOS` 或 `&PROJWFC` 中给了 `degauss` | 删除 `degauss`，让程序使用四面体 |
| dos.x 读 namelist 报错 | 用了旧的 `&INPUTPP` | 改为 `&DOS` |
| `sumpdos.x` 报 `No file to sum` 或找不到文件 | shell 没有正确展开含括号的文件名 | 用 `\(`、`\)` 转义，或把文件名写入列表用 `-f` 传入 |
| PDOS 之和明显小于总 DOS | 原子轨道基不完备（spilling） | 正常现象；查看 `Spilling Parameter`，必要时换带更多投影轨道的赝势 |
| 能带图中费米能位置不对 | 用了 bands 计算输出的费米能 | 从 scf/nscf 输出或 `si.dos` 表头读取 |
| 在 nscf 之后运行 bands.x 得到奇怪的图 | nscf 覆盖了 bands 写入的本征值 | 重新运行 bands 再运行 bands.x，或复制 `tmp` 目录分开做 |

---

## 练习

1. 运行完整流程，从能带数据中读出 Si 的价带顶位置、导带底在 Γ–X 线上的相对位置和间接带隙；与 scf 的 `highest occupied, lowest unoccupied level` 比较并解释差异。
2. 把 nscf 网格分别设为 8×8×8、16×16×16、24×24×24，比较 DOS 曲线；再把 `occupations` 换成 `'smearing'` 并在 `&DOS` 中给 `degauss=0.01`，比较峰形。
3. 检查 `si.dos` 第三列（积分 DOS）在价带顶处的数值是否等于 8。
4. 用 `tpiba_b` 重写 `si.bands.in` 的 K_POINTS（坐标见第 5.2 节表格第二列），确认两种写法得到相同的能带。
5. 对 fcc Cu 完成 scf（smearing）→ bands → nscf（`tetrahedra_opt`）→ projwfc.x，指出 d 带的能量范围与费米能附近的主要轨道成分。
6. 对 GaAs 计算 PDOS，用 sumpdos.x 分别得到 Ga 与 As 的 s、p、d 贡献，判断价带顶与导带底分别由哪种元素主导。

---

## 参考资料

**官方文档**

- [INPUT_PW（pw.x）](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — `calculation='bands'/'nscf'`、`nbnd`、`occupations`、`K_POINTS crystal_b/tpiba_b`、`diago_full_acc`
- [INPUT_BANDS（bands.x）](https://www.quantum-espresso.org/Doc/INPUT_BANDS.html) — `filband`、`lsym`、`spin_component`、`lp`
- [INPUT_DOS（dos.x）](https://www.quantum-espresso.org/Doc/INPUT_DOS.html) — `bz_sum`、`degauss`、`Emin`/`Emax`/`DeltaE`、展宽规则
- [INPUT_PROJWFC（projwfc.x）](https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html) — 文件格式、轨道顺序、`pawproj`、`kresolveddos`
- [PWscf User's Guide](https://www.quantum-espresso.org/Doc/pw_user_guide/)
- [SeeK-path 在线工具](https://www.materialscloud.org/work/tools/seekpath) — 标准原胞与高对称路径，可导出 QE 格式

**课程与教程（推荐，链接而非复制）**

- [Imperial College MSE404-MM Lab 4: Crystals and the Electronic Band Structure](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab04/readme/) — QE 能带计算的清晰入门
- [Imperial College MSE404-MM Lab 5: Metals and the Density of States](https://imperialcollegelondon.github.io/MSE404-MM/labs/lab05/readme/) — 金属 DOS 与展宽
- [Cornell PARADIM Summer School 2016](https://www.paradim.org/cu_ss_2016) — 能带与 Kohn-Sham 波函数可视化上机（QE 5.x/6.x 时代，原理部分仍适用）
- [Cornell PARADIM Summer School 2018](https://www.paradim.org/cu_ss_2018) — Si 能带与 GaAs UV-Vis，每讲有视频
- [pranabdas/espresso: Band structure](https://pranabdas.github.io/espresso/hands-on/bands)、[DOS](https://pranabdas.github.io/espresso/hands-on/dos)、[PDOS](https://pranabdas.github.io/espresso/hands-on/pdos) — 附 Python 绘图
- [MaX e-School 2021（ICTP）](https://indico.ictp.it/event/9616/) — QE 官方开发者授课，含视频

**相关仓库（仅链接）**

- [materialscloud-org/seekpath](https://github.com/materialscloud-org/seekpath) — 高对称路径生成
- [yyyu200/QEbandplot](https://github.com/yyyu200/QEbandplot) — 读取 bands.x 输出作图
- [shubics/QEPlotter](https://github.com/shubics/QEPlotter) — 能带、PDOS、fatband 绘图
- [Griffin-Group/pymatgen-io-espresso](https://github.com/Griffin-Group/pymatgen-io-espresso) — 把 QE XML 转成类 Vasprun 对象，便于复用 pymatgen 绘图
- [QEF/postqe](https://github.com/QEF/postqe) — 官方 Python 后处理（能带、DOS）
- [band-unfolding/bandup](https://github.com/band-unfolding/bandup) — 超胞能带反折叠

**关键文献**

- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* 49, 299 (2010). DOI: [10.1016/j.commatsci.2010.05.010](https://doi.org/10.1016/j.commatsci.2010.05.010)
- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* 128, 140 (2017)（SeeK-path）. DOI: [10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015)
- P. E. Blöchl, O. Jepsen, O. K. Andersen, *Phys. Rev. B* 49, 16223 (1994)（四面体方法）. DOI: [10.1103/PhysRevB.49.16223](https://doi.org/10.1103/PhysRevB.49.16223)
- M. Kawamura, Y. Gohda, S. Tsuneyuki, *Phys. Rev. B* 89, 094515 (2014)（优化四面体方法）. DOI: [10.1103/PhysRevB.89.094515](https://doi.org/10.1103/PhysRevB.89.094515)
- P. E. Blöchl, *Phys. Rev. B* 50, 17953 (1994)（PAW，`pawproj` 依据）. DOI: [10.1103/PhysRevB.50.17953](https://doi.org/10.1103/PhysRevB.50.17953)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017). DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)

---

[← 上一章](05-relax-eos.md) | [返回目录](../README.md) | [下一章 →](07-magnetism-soc.md)
