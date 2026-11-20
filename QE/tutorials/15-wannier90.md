# 第15章 Wannier90 接口与插值能带（Maximally Localized Wannier Functions with Wannier90）

本章讲 Quantum ESPRESSO 与 Wannier90 的标准接口：pw.x scf → pw.x nscf（显式全网格 `K_POINTS crystal`、`nosym`/`noinv`）→ `wannier90.x -pp` 生成 `seedname.nnkp` → `pw2wannier90.x` 计算重叠矩阵 Mmn 与投影 Amn → `wannier90.x` 求最大局域 Wannier 函数（MLWF）。重点是 `seedname.win` 中的 `num_wann`、`num_bands`、`projections`、`exclude_bands`、解纠缠窗口（`dis_win_max`、`dis_froz_max`）、`mp_grid`/`kpoints`、`bands_plot`，以及如何从 `seedname.wout` 判断 spread 收敛（`Final State`）。示例为 Si 价带（4 个键中心 MLWF）和 Si 价带 + 导带（8 个 sp3 MLWF，解纠缠），并简要介绍 Cu/石墨烯的解纠缠思路、SCDM 自动投影和 WannierTools 拓扑分析。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 Wannier 化所需的 nscf 输入，用 `kmesh.pl` 生成全网格 k 列表，并让它与 `.win` 的 `kpoints` 块一一对应；
- 写出 `seedname.win`：选择 `num_wann`/`num_bands`、写 `projections` 块、设置解纠缠窗口；
- 正确执行 `wannier90.x -pp` → `pw2wannier90.x` → `wannier90.x` 三步，知道每一步读写哪些文件；
- 从 `seedname.wout` 读 Wannier 中心、spread、`Final State`，判断结果是否合理；
- 用插值能带与 DFT 能带对比来验证 Wannier 化质量；
- 知道 SCDM 自动投影、`write_hr`（紧束缚哈密顿量）、WannierTools/postw90 等后续用法。

| 目标 | 做法 | 本章位置 |
|---|---|---|
| 化学成键分析（键中心、孤对） | 孤立价带流形 MLWF | §4 |
| 任意 k 路径上的平滑能带 | `bands_plot = .true.` | §4、§5 |
| 杂化泛函能带 | HSE 全网格自洽后 Wannier 插值 | [第14章 §6](14-hybrid-metagga.md) |
| 紧束缚模型 / 输运 / 拓扑 | `write_hr`、`write_tb` → WannierTools 等 | §7 |
| 金属费米面附近能带 | 解纠缠（`dis_win_*`、`dis_froz_*`） | §5、§6 |
| 无法猜测投影 | SCDM（`auto_projections` + `scdm_proj`） | §6.3 |
| 电声耦合插值 | EPW（内部调用 Wannier90） | [第16章](16-elph-epw.md) |

---

## 2. 原理速览

### 2.1 Wannier 函数与规范自由度

对一组（孤立的）能带 n = 1…N，Bloch 态 ψ_{nk} 可以在每个 k 点做任意幺正变换 U^{(k)} 而不改变物理。Wannier 函数定义为

$$
w_{n\mathbf R}(\mathbf r) = \frac{V}{(2\pi)^3}\int_{\mathrm{BZ}} d\mathbf k\; e^{-i\mathbf k\cdot\mathbf R}\sum_m U^{(\mathbf k)}_{mn}\,\psi_{m\mathbf k}(\mathbf r)
$$

Marzari–Vanderbilt 方法选择 U^{(k)} 使总 spread 最小：

$$
\Omega = \sum_n \left[\langle w_{n0}|r^2|w_{n0}\rangle - |\langle w_{n0}|\mathbf r|w_{n0}\rangle|^2\right] = \Omega_I + \widetilde\Omega,
\qquad \widetilde\Omega = \Omega_D + \Omega_{OD}
$$

Ω_I 与规范无关；最小化只作用于 Ω̃。在离散 k 网格上，所需的全部信息是相邻 k 点周期部分的重叠矩阵

$$
M^{(\mathbf k,\mathbf b)}_{mn} = \langle u_{m\mathbf k}|u_{n,\mathbf k+\mathbf b}\rangle
$$

以及初始投影

$$
A^{(\mathbf k)}_{mn} = \langle \psi_{m\mathbf k}|g_n\rangle
$$

（g_n 为用户给定的试探轨道）。`pw2wannier90.x` 的工作就是计算 M（写入 `seedname.mmn`）、A（`seedname.amn`）和本征值（`seedname.eig`）。

### 2.2 解纠缠

当目标能带与其他能带交叉（金属、导带）时，需要先从更多的能带（`num_bands`）中挑出一个 `num_wann` 维的"最优连接子空间"，即最小化 Ω_I（Souza–Marzari–Vanderbilt）。两个能量窗口控制这一过程：

- **外窗口**（`dis_win_min`/`dis_win_max`）：只有窗口内的态参与；
- **冻结窗口**（`dis_froz_min`/`dis_froz_max`）：窗口内的态原样保留，插值能带在此范围内与 DFT 精确一致。

Wannier90 文档说明：`num_wann < num_bands` 时才启用解纠缠；冻结窗口内的态数在任何 k 点都不能超过 `num_wann`。

### 2.3 为什么能做插值

一旦得到局域的 w_{nR}，实空间哈密顿量 H_{mn}(R) = ⟨w_{m0}|H|w_{nR}⟩ 随 |R| 快速衰减，于是任意 k 点的哈密顿量可以由傅里叶求和得到，再对角化得到插值能带。这就是 `bands_plot` 的原理，也是 EPW、BoltzWann、WannierTools 等工具的基础。

---

## 3. 准备工作

示例目录：[../examples/15-wannier90/README.md](../examples/15-wannier90/README.md)。

### 3.1 程序

| 程序 | 来源 | 说明 |
|---|---|---|
| `pw.x` | QE | scf 与 nscf |
| `pw2wannier90.x` | QE 的 PP 包（`make pp`） | 读 `.nnkp`，写 `.mmn`、`.amn`、`.eig`（以及可选的 `UNK*`） |
| `wannier90.x` | Wannier90 | 预处理（`-pp`）与 Wannier 化 |
| `postw90.x` | Wannier90 | 可选：DOS、Berry 曲率、BoltzWann 等 |

安装方式：

- QE 源码树中 `make w90` 会编译随附的 Wannier90（QE 7.6 的 `external/wannier90` 子模块指向 Wannier90 **v3.1.0**）；
- 独立安装最新版 Wannier90（v4.0.x：v4.0.1 于 2026-07 发布，截至 2026-11 最新为 v4.0.3；conda-forge 上的 `wannier90` 包也是 4.0.x）。本章的 `.win` 关键词在 v3.1 与 v4 中通用；`dis_froz_proj` 等可投影性解纠缠关键词在 QE 随附的 v3.1.0 中不存在；
- `pw2wannier90.x` 与 `wannier90.x` 之间通过文件交互（standalone 模式），两者版本不必严格一致。QE 开发版（8.0 预告）中 `wan_mode='library'` 将改为调用 Wannier90 v4 库并行 Wannier 化，正式发布前不要依赖。

```bash
which pw.x pw2wannier90.x wannier90.x
wannier90.x -v     # prints "Wannier90: <version>" (-v/--version, also in v3.1.0)
```

### 3.2 赝势

本章示例用 SSSP 1.3 PBE Efficiency 中的 Si 超软赝势 `Si.pbe-n-rrkjus_psl.1.0.0.UPF`（推荐 `ecutwfc = 30`、`ecutrho = 240` Ry）。`pw2wannier90.x` 支持 NC、US、PAW；源码头部注明了 Koretsune 与 Thoele 对非共线和 USPP 的扩展。下载方法见 [第3章](03-pseudopotentials.md)。

---

## 4. Si 价带：孤立流形的 4 个 MLWF

### 4.1 第 1 步：scf

```fortran
! Si (diamond), PBE, SSSP Efficiency 1.3 pseudopotential (ultrasoft).
! Step 1 of the Wannier90 workflow: self-consistent charge density.
&CONTROL
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 30.0      ! SSSP 1.3 Efficiency recommendation for Si
  ecutrho   = 240.0     ! dual 8 (ultrasoft)
/
&ELECTRONS
  conv_thr  = 1.0d-10
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  -2.7150   0.0000   2.7150
   0.0000   2.7150   2.7150
  -2.7150   2.7150   0.0000
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

```bash
mpirun -np 8 pw.x -nk 4 -in si.scf.in > si.scf.out
grep "highest occupied" si.scf.out
```

记下价带顶（VBM）的绝对能量（eV），§5 设置解纠缠窗口时要用。

### 4.2 第 2 步：nscf（全网格、显式 k 列表）

Wannier90 需要**整个布里渊区**的均匀网格（Γ 中心 Monkhorst–Pack），并且 k 点顺序必须与 `.win` 中的 `kpoints` 块一致。用 Wannier90 自带的 `utility/kmesh.pl` 生成：

```bash
perl kmesh.pl 4 4 4          > kpts_qe.txt    # "K_POINTS crystal" + 64 lines with weights -> for pw.x
perl kmesh.pl 4 4 4 wannier  > kpts_w90.txt   # 64 lines without weights -> for the .win kpoints block
```

`kmesh.pl` 的用法（据源码）：`n1 n2 n3 [wan]`，给第 4 个参数时省略权重列。示例目录中两份列表已经写好。

```fortran
! Step 2: non-SCF on the FULL 4x4x4 grid, listed explicitly (generated as by
! Wannier90 utility/kmesh.pl 4 4 4). The order must match the kpoints block in si.win.
&CONTROL
  calculation = 'nscf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  verbosity   = 'high'     ! print eigenvalues at every k point
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 30.0
  ecutrho   = 240.0
  nbnd      = 12         ! 4 valence + 8 conduction bands
  nosym     = .true.     ! do not reduce or reorder the k list
  noinv     = .true.     ! do not use time reversal to reduce the k list
/
&ELECTRONS
  conv_thr        = 1.0d-10
  diago_full_acc  = .true.   ! empty states converged as tightly as occupied ones
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  -2.7150   0.0000   2.7150
   0.0000   2.7150   2.7150
  -2.7150   2.7150   0.0000
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS crystal
64
  0.00000000  0.00000000  0.00000000  1.562500e-02
  0.00000000  0.00000000  0.25000000  1.562500e-02
  ...（共 64 行，完整列表见示例目录 si.nscf.in）
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `calculation` | `'nscf'` | 在 scf 电荷密度上求解更多能带；也可用 `'bands'`，但 nscf 更常见 |
| `K_POINTS crystal` | `kmesh.pl` 生成的全网格 | 不要用 `automatic`：pw.x 会按对称性约化 |
| `nosym` | `.true.` | INPUT_PW：显式列表"as is"使用，不生成对称等价点 |
| `noinv` | `.true.` | 不用 k ↔ −k 时间反演约化 |
| `nbnd` | ≥ `num_bands` + `exclude_bands` 的数目 | 本例 12 |
| `diago_full_acc` | `.true.` | INPUT_PW：空带与占据带同等精度对角化 |
| `verbosity` | `'high'` | 打印所有 k 点本征值，便于与插值能带对比 |

```bash
mpirun -np 8 pw.x -nk 4 -in si.nscf.in > si.nscf.out
```

> [!IMPORTANT]
> nscf 必须与 scf 使用相同的 `prefix`、`outdir`、赝势与截断；nscf 会覆盖 save 目录中的波函数（电荷密度不变）。nscf 网格（这里 4×4×4）与 scf 网格（8×8×8）可以不同：Wannier 网格决定插值质量和计算量，scf 网格只决定电荷密度。

### 4.3 第 3 步：seedname.win

```text
! si_val.win -- 4 MLWFs for the isolated valence manifold of Si (no disentanglement)
num_wann      = 4
num_bands     = 4
exclude_bands = 5-12      ! the nscf run has 12 bands; keep only the 4 valence bands
num_iter      = 100
conv_window   = 3
conv_tol      = 1.0e-10

bands_plot       = .true.
bands_num_points = 100
write_xyz        = .true.

! one s-like trial orbital at the centre of each Si-Si bond (fractional coordinates)
begin projections
f= 0.125, 0.125, 0.125:s
f=-0.375, 0.125, 0.125:s
f= 0.125,-0.375, 0.125:s
f= 0.125, 0.125,-0.375:s
end projections

begin unit_cell_cart
ang
 -2.7150   0.0000   2.7150
  0.0000   2.7150   2.7150
 -2.7150   2.7150   0.0000
end unit_cell_cart

begin atoms_frac
Si  0.00  0.00  0.00
Si  0.25  0.25  0.25
end atoms_frac

! fractional coordinates w.r.t. the reciprocal vectors of the cell above
begin kpoint_path
L  0.000  0.500  0.000   G  0.000  0.000  0.000
G  0.000  0.000  0.000   X  0.000  0.500  0.500
X  0.000  0.500  0.500   U  0.000  0.625  0.375
K -0.375  0.375  0.000   G  0.000  0.000  0.000
end kpoint_path

mp_grid = 4 4 4

begin kpoints
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.25000000
  ...（共 64 行，与 si.nscf.in 顺序相同）
end kpoints
```

关键词（据 Wannier90 用户指南 Parameters 一节）：

| 关键词 | 本例 | 说明 |
|---|---|---|
| `num_wann` | 4 | Wannier 函数个数 |
| `num_bands` | 4 | 传给 Wannier90 的能带数（排除 `exclude_bands` 之后）；等于 `num_wann` 时不做解纠缠 |
| `exclude_bands` | `5-12` | 不参与 Wannier 化的能带编号（支持 `2, 6-8, 12` 这种写法） |
| `num_iter` | 100（默认 100） | 最小化 Ω̃ 的迭代次数；设为 0 得到投影 WF（不做最大局域化） |
| `conv_window` / `conv_tol` | 3 / 1e-10 | `conv_window > 1` 时，若连续 `conv_window` 步 Ω 变化都小于 `conv_tol`（Å²）即认为收敛；默认 `conv_window = -1`，跑满 `num_iter` |
| `projections` 块 | 键中心 s 轨道 | 初始投影；`f=` 为分数坐标，`c=` 为笛卡尔坐标（单位由块首行 `Ang`/`Bohr` 指定，默认 Å） |
| `unit_cell_cart` | 与 pw.x 相同 | 晶格矢量；首行可写单位 |
| `atoms_frac` / `atoms_cart` | — | 原子位置（二者只能给一个） |
| `mp_grid` | `4 4 4` | 必须与 nscf 网格一致 |
| `kpoints` 块 | 64 行 | 与 nscf 的 k 列表**同序** |
| `bands_plot` | `.true.` | 输出插值能带 |
| `kpoint_path` | L–Γ–X–U / K–Γ | 每行一段：起点标签 + 坐标、终点标签 + 坐标 |
| `bands_num_points` | 100（默认 100） | 第一段上的点数，其余段按长度等密度 |
| `write_xyz` | `.true.` | 写 `seedname_centres.xyz`（原子 + WF 中心） |

键中心坐标的由来：Si 原子位于 0 和 τ = (¼,¼,¼)（晶体坐标），四个最近邻分别是 τ − 0、τ − a1、τ − a2、τ − a3，键中心是 (τ − a_i)/2，换算回这组原胞矢量的分数坐标就是上面四个点。Wannier90 官方 tutorial 11（Si 价带与导带）的价带部分也得到位于 Si–Si 键中心的 MLWF，但它的输入文件使用自己的晶胞约定，坐标不能直接照抄。

> [!TIP]
> 也可以用原子中心投影 `Si:sp3`（每个 Si 4 个 sp3，共 8 个）再配合 `select_projections` 选其中 4 个；或者干脆用 SCDM（§6.3）自动生成初始投影。

### 4.4 第 4 步：wannier90.x -pp

```bash
wannier90.x -pp si_val          # reads si_val.win, writes si_val.nnkp
grep -A3 "begin projections" si_val.nnkp
grep -A2 "begin exclude_bands" si_val.nnkp
```

`seedname.nnkp` 包含：晶格与倒格矢、k 点列表、每个 k 点的近邻 k+b（`nnkpts` 块）、投影信息（`projections` 块或 `auto_projections` 块）、被排除的能带。`pw2wannier90.x` 完全依据这个文件工作。

### 4.5 第 5 步：pw2wannier90.x

```fortran
&INPUTPP
  outdir    = './tmp'
  prefix    = 'si'
  seedname  = 'si_val'
  write_mmn = .true.
  write_amn = .true.
  write_unk = .false.     ! set .true. only to plot Wannier functions (large files)
/
```

| 关键词（`&inputpp`） | 默认值（QE 7.6 源码） | 说明 |
|---|---|---|
| `outdir` / `prefix` | `$ESPRESSO_TMPDIR` 或 `./` / 空 | 与 nscf 相同 |
| `seedname` | `'wannier'` | 与 `.win`/`.nnkp` 文件名一致 |
| `write_mmn` | `.true.` | 写重叠矩阵 `seedname.mmn` |
| `write_amn` | `.true.` | 写投影 `seedname.amn` |
| `write_eig` | `.true.` | 写本征值 `seedname.eig` |
| `write_unk` | `.false.` | 写周期部分 `UNKxxxxx.1`，画 WF 时需要（文件很大） |
| `reduce_unk` | `.false.` | 降低 UNK 文件的实空间分辨率 |
| `wvfn_formatted` | `.false.` | UNK 是否用格式化文件 |
| `spin_component` | `'none'` | 自旋极化时选择 `'up'` 或 `'down'` |
| `write_spn` | `.false.` | 写自旋矩阵（SOC + postw90 时用） |
| `wan_mode` | `'standalone'` | 与独立 wannier90.x 通过文件交互 |
| `scdm_proj` | `.false.` | 用 SCDM 生成投影（§6.3） |
| `scdm_entanglement` | `'isolated'` | `'isolated'`、`'erfc'`、`'gaussian'` |
| `scdm_mu` / `scdm_sigma` | 0.0 / 1.0 | `erfc`/`gaussian` 时的参数（eV） |
| `atom_proj` | `.false.` | 用赝势中的原子轨道作投影（配合 `auto_projections`） |

```bash
mpirun -np 8 pw2wannier90.x -in si_val.pw2wan.in > si_val.pw2wan.out
grep -e "Reciprocal lattice is ok" -e "K-points are ok" si_val.pw2wan.out
ls si_val.mmn si_val.amn si_val.eig
```

`pw2wannier90.x` 会先核对 `.nnkp` 中的倒格矢和 k 点与 save 目录是否一致（输出 `- Reciprocal lattice is ok`、`- K-points are ok`）。它可以 MPI 并行（pools 也可），但不支持 `-nb`。

### 4.6 第 6 步：wannier90.x

```bash
wannier90.x si_val
grep "<-- CONV" si_val.wout | tail -5
grep -A8 "Final State" si_val.wout
grep "All done" si_val.wout
```

`seedname.wout` 中的关键内容（示意输出，数值仅供参考）：

```text
 +--------------------------------------------------------------------+<-- CONV
 | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
 +--------------------------------------------------------------------+<-- CONV
      0     0.xxxE+01     0.0000000000        x.xxxxxxxxxx       0.xx  <-- CONV
      ...
 Final State
  WF centre and spread    1  (  x.xxxxxx,  x.xxxxxx,  x.xxxxxx )     x.xxxxxxxx
  WF centre and spread    2  (  x.xxxxxx,  x.xxxxxx,  x.xxxxxx )     x.xxxxxxxx
  ...
  Sum of centres and spreads ( x.xxxxxx, x.xxxxxx, x.xxxxxx )     x.xxxxxxxx

         Spreads (Ang^2)       Omega I      =     x.xxxxxxxxx
        ================       Omega D      =     0.xxxxxxxxx
                               Omega OD     =     x.xxxxxxxxx
    Final Spread (Ang^2)       Omega Total  =     x.xxxxxxxxx
 ...
 All done: wannier90 exiting
```

成功判据：

- `<-- CONV` 行中 Delta Spread 逐步减小，Spread 稳定；
- `Final State` 中 4 个 WF 中心位于 4 个键中心附近（可在 `si_val_centres.xyz` 中用 VESTA/Jmol 查看），4 个 spread 相等；
- Ω_D 很小（孤立流形中 Ω_D 可被优化到接近 0）；
- 插值能带 `si_val_band.dat` 与 DFT 价带（[第6章](06-band-dos-pdos.md) 的 `bands` 计算，或 nscf 网格点上的本征值）重合。

画图：

```bash
gnuplot -persist si_val_band.gnu                      # quick look (generated by wannier90)
cat si_val_band.labelinfo.dat                         # positions of the high-symmetry labels
```

也可以把 `si_val_band.dat`（两列：路径坐标、能量 eV，每条带之间空行分隔）与第 6 章的 DFT 能带画在一起；`qe_plot_bands.py` 的用法请在示例目录下用 `python ../../scripts/qe_plot_bands.py -h` 查看。注意 Wannier90 与 bands.x 的路径横坐标归一化可能不同，叠加前要统一。

---

## 5. Si 价带 + 导带：解纠缠

### 5.1 输入

导带与更高能带纠缠在一起。取每个 Si 4 个 sp3（共 8 个 WF），用 12 条能带解纠缠：

```text
! si.win -- 8 MLWFs (sp3 on both Si) for valence + low conduction bands, with disentanglement
num_wann      = 8
num_bands     = 12
num_iter      = 200
conv_window   = 3
conv_tol      = 1.0e-10

! energy windows in eV, ABSOLUTE energies of the pw.x run.
! Set dis_froz_max ~1-2 eV above the VBM printed by YOUR si.scf.out
dis_win_max   = 17.0
dis_froz_max  = 7.0
dis_num_iter  = 200

bands_plot       = .true.
bands_num_points = 100
write_hr         = .true.
write_xyz        = .true.

begin projections
Si:sp3
end projections
（unit_cell_cart、atoms_frac、kpoint_path、mp_grid、kpoints 与 si_val.win 相同）
```

| 关键词 | 本例 | 说明（据 Wannier90 用户指南） |
|---|---|---|
| `dis_win_min` | 不设 | 外窗口下限，默认为最低本征值 |
| `dis_win_max` | 17.0 | 外窗口上限（eV），默认为最高本征值 |
| `dis_froz_min` | 不设 | 冻结窗口下限；给了 `dis_froz_max` 时默认等于 `dis_win_min` |
| `dis_froz_max` | VBM + 1~2 eV | 冻结窗口上限；不设则没有冻结态 |
| `dis_num_iter` | 200（默认 200） | Ω_I 最小化的迭代次数 |
| `dis_mix_ratio` | 默认 0.5 | 解纠缠混合比 |
| `dis_conv_tol` / `dis_conv_window` | 默认 1e-10 / 3 | Ω_I 收敛判据 |
| `write_hr` | `.true.` | 写 `seedname_hr.dat`（Wannier 基下的实空间哈密顿量） |

> [!WARNING]
> 窗口能量是 pw.x 输出的**绝对**本征值（eV），不是相对费米能或价带顶。换赝势、晶格常数或泛函后，VBM 的绝对位置会变，窗口必须重新设置。可以用 `grep "highest occupied" si.scf.out` 或看 nscf 输出中的本征值确定。

### 5.2 运行与检查

```bash
wannier90.x -pp si
mpirun -np 8 pw2wannier90.x -in si.pw2wan.in > si.pw2wan.out
wannier90.x si
grep "<-- DIS" si.wout | tail -3          # disentanglement convergence (Omega_I)
grep "Final Omega_I" si.wout
grep -A12 "Final State" si.wout
ls si_hr.dat si_band.dat si_centres.xyz
```

判据：

- `<-- DIS` 表中 Ω_I 收敛（Delta 很小）；
- 8 个 WF 中心位于两个 Si 原子附近，呈 sp3 形状；
- 冻结窗口以下（价带与最低导带）的插值能带与 DFT 精确重合；冻结窗口以上到外窗口上限之间允许小偏差；
- 如果插值能带出现 DFT 中没有的"摆动"或尖峰，通常是窗口设置不当、投影选择不当或 k 网格太粗。

---

## 6. 进阶：金属、二维材料与自动投影

### 6.1 Cu：d 带与费米面附近

思路（不提供完整输入，见参考资料中 Wannier90 tutorial 4/6）：fcc Cu 可取每个 Cu 的 `s;d` 或 `s;p;d` 投影（6 或 9 个 WF），`num_bands` 取更多（例如 12–20），冻结窗口上限设在 E_F 以上 1–2 eV，使费米面附近能带被精确保留。nscf 与 scf 都要用 smearing（金属）。检查费米能附近插值能带是否与 DFT 一致，以及 `Final State` spread 是否合理。

### 6.2 石墨烯：二维与 pz 轨道

石墨烯的 π 带由 C 的 pz 轨道构成，最小模型是每个 C 一个 `C:pz`（2 个 WF），但它们与 σ 带纠缠，需要解纠缠；更稳健的选择是 `C:sp2;pz`（每个 C 4 个）。二维体系 nscf 网格取 N×N×1（`mp_grid` 第三维为 1）。Wannier90 官方 tutorial 34 用"可投影性解纠缠"（projectability disentanglement，`dis_froz_proj` + `dis_proj_min`/`dis_proj_max`，需要 `pw2wannier90.x` 用 `atom_proj` 生成归一化的 `.amn`）处理石墨烯，比手动调窗口更自动化。QE 7.6 随附的 Wannier90 v3.1.0 没有这些关键词，需要使用独立安装的新版 Wannier90。

### 6.3 SCDM 自动投影

不知道该用什么投影时，可以用 SCDM（selected columns of the density matrix）由密度矩阵自动构造初始投影：

1. `.win` 中**不写** `projections` 块，加上 `auto_projections = .true.`；
2. `wannier90.x -pp` 会在 `.nnkp` 中写 `auto_projections` 块；
3. `pw2wannier90.x` 输入中设 `scdm_proj = .true.`，按体系选择 `scdm_entanglement`：

```fortran
&INPUTPP
  outdir   = './tmp'
  prefix   = 'si'
  seedname = 'si_scdm'
  scdm_proj         = .true.
  scdm_entanglement = 'erfc'     ! 'isolated' for an isolated manifold
  scdm_mu           = 10.0       ! eV, illustrative: choose from your band structure
  scdm_sigma        = 2.0        ! eV, illustrative
/
```

`scdm_entanglement` 的合法值（QE 7.6 源码检查）：`'isolated'`、`'erfc'`、`'gaussian'`。`scdm_mu`/`scdm_sigma` 的选择方法见 Wannier90 tutorial 27 与 32（后者讲如何从可投影性确定参数），以及 Vitale 等人的高通量 Wannier 化论文。

> [!NOTE]
> `pw2wannier90.x` 的一致性检查：同时出现 `projections` 块、`auto_projections` 块和 `scdm_proj = .true.` 会报 `Inconsistent options for projections.`；`scdm_proj = .true.` 但 `.nnkp` 中没有 `auto_projections` 块会报 `scdm_proj = T but cannot find an auto_projections block`。

### 6.4 自旋极化与 SOC

- 共线自旋（`nspin = 2`）：两个自旋通道分别 Wannier 化，各用一个 seedname，`pw2wannier90.x` 中分别设 `spin_component = 'up'`/`'down'`；`.win` 中的 `spin` 关键词（`up`/`down`）只影响 Wannier90 读取哪一个自旋的能带（用于 postw90 等），见 Wannier90 用户指南；
- 非共线/SOC：`.win` 中 `spinors = .true.`，投影可写成带自旋的形式（例如 `Fe:d(u,d)`），`num_wann` 翻倍。详见 Wannier90 tutorial 17（Fe 的 SOC 能带与费米面）。

---

## 7. Wannier 函数之后：常见后续用法

| 用途 | 需要的输出 | 工具 |
|---|---|---|
| 紧束缚哈密顿量 | `seedname_hr.dat`（`write_hr`）、`seedname_tb.dat`（`write_tb`） | 自写脚本、PythTB、TBmodels 等 |
| 表面态、Wilson loop、Weyl 点、Z2 | `seedname_hr.dat` | [WannierTools](https://github.com/quanshengwu/wannier_tools) |
| DOS、Berry 曲率、反常霍尔、BoltzWann 输运 | `seedname.chk` 等 | `postw90.x` |
| 画 Wannier 函数 | `write_unk = .true.` + `wannier_plot = .true.` | VESTA、XCrySDen（tutorial 28 讲 cube 格式） |
| 电声耦合 | EPW 内部调用 Wannier90 | [第16章](16-elph-epw.md) |
| 杂化泛函能带 | HSE 全网格自洽 + Wannier 插值 | [第14章 §6](14-hybrid-metagga.md) |
| Wannier 函数作为 Hubbard 投影 | QE 7.5 起支持 | [第8章](08-dftu-hp.md) |

---

## 8. 与 VASP 的对应关系

VASP 一侧见 [VASP 第17章 高级电子结构](../../VASP/tutorials/17-advanced-electronic.md)。

| 步骤 | VASP | Quantum ESPRESSO |
|---|---|---|
| 生成 Wannier90 输入文件 | `LWANNIER90 = .TRUE.` 自动写 `wannier90.win` | 手写 `seedname.win` |
| 投影 | `WANNIER90_WIN` 字符串或 `.win` | `.win` 的 `projections` 块 |
| k 点 | VASP 需 `ISYM = -1` 的全网格 | nscf 显式全网格 + `nosym`/`noinv` |
| Mmn/Amn 计算 | VASP 内部（库模式） | `wannier90.x -pp` + `pw2wannier90.x` |
| Wannier 化 | VASP 内部调用或独立 `wannier90.x` | 独立 `wannier90.x` |
| 能带数 | `NBANDS`、`NUM_WANN` | `nbnd`（pw.x）、`num_bands`/`num_wann`（.win） |
| SCDM | `LSCDM` | `auto_projections` + `scdm_proj` |

---

## 参数选择与收敛建议

| 参数 | 建议 | 判据 |
|---|---|---|
| nscf 网格（`mp_grid`） | 绝缘体 4³–6³；金属/费米面 8³ 以上 | 插值能带在网格点之间不再变化；`seedname_hr.dat` 远端矩阵元足够小 |
| `nbnd`（nscf） | ≥ `num_bands` + 排除的能带 | 外窗口上限以下的能带完整 |
| `num_bands` | 孤立流形 = `num_wann`；解纠缠时多出若干条 | 外窗口内有足够的态 |
| `dis_froz_max` | VBM 或 E_F 以上 1–2 eV | 冻结窗口内态数 ≤ `num_wann`（每个 k 点） |
| `dis_win_max` | 能覆盖目标能带的上沿 | 插值能带无伪振荡 |
| `num_iter` | 100–500 | `<-- CONV` 中 spread 稳定 |
| `conv_window` | 3–5 | 提前结束已收敛的最小化 |
| 投影 | 化学直觉（键中心、sp3、d）；不确定时 SCDM | 初始 spread 合理、最终中心位置合理 |
| 截断与赝势 | 与 scf 相同 | — |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 处理 |
|---|---|---|
| `pw2wannier90`：`Wrong number of k-points` | nscf 的 k 点数与 `.nnkp` 不一致（多半是 nscf 用了对称性约化或 `automatic`） | nscf 用 `kmesh.pl` 生成的 `K_POINTS crystal` + `nosym`/`noinv` |
| `pw2wannier90`：`problems with k-points` | k 点顺序或坐标与 `.nnkp` 不同 | `.win` 的 `kpoints` 与 nscf 列表用同一工具、同一顺序生成 |
| `Could not find the file seedname.nnkp` | 没有先运行 `wannier90.x -pp`，或 `seedname` 与 `.win` 文件名不一致 | 先跑 `-pp`；检查 `seedname` |
| `Could not find kpoints block in seedname.nnkp` | `.nnkp` 文件不完整（例如 `-pp` 中途出错） | 检查 `-pp` 的 `.wout`，重新生成 `.nnkp` |
| `wannier90`：`dis_windows: More states in the frozen window than target WFs` | 某个 k 点上冻结窗口内能带数 > `num_wann` | 降低 `dis_froz_max` 或增加 `num_wann` |
| `pw2wannier90`：`bands (-nb) not implemented` | 用了 `-nb` | 去掉 `-nb` |
| `pw2wannier90`：`pools not implemented for library mode` | `wan_mode='library'` + `-nk` | standalone 模式或不用 pools |
| `Inconsistent options for projections.` | 投影块、auto_projections、scdm_proj 同时出现 | 只保留一种投影方式 |
| `scdm_proj = T but cannot find an auto_projections block` | `.win` 中没有 `auto_projections = .true.` | 在 `.win` 中加上并重新 `-pp` |
| `Can not recognize the choice for scdm_entanglement` | 拼写错误 | 只能是 `isolated`、`erfc`、`gaussian` |
| spread 很大且不收敛 | 投影选择不当（中心远离电子密集区）、窗口错误、k 网格太粗 | 换投影或用 SCDM；检查窗口；加密网格 |
| 插值能带在高对称点附近有尖峰/振荡 | 解纠缠不充分；冻结窗口过高或过低 | 调整窗口；增加 `dis_num_iter`；增大 `num_bands` |
| 插值能带整体平移 | 与 DFT 对比时一方减了费米能 | Wannier90 输出的是绝对能量，统一参考零点 |
| 杂化泛函下 nscf 失败 | 杂化泛函不能在新 k 点 nscf | 直接在自洽中用显式全网格，见第 14 章 |

---

## 练习

1. 把 `si_val` 的投影改为 `Si:sp3` 并用 `select_projections` 选择其中 4 个，比较 `Final State` 的中心与 spread。
2. 把 nscf 网格从 4×4×4 增加到 6×6×6（同时改 `mp_grid` 和 `kpoints`），比较插值价带的变化。
3. 对 `si.win` 扫描 `dis_froz_max`（VBM + 0.5、1.0、2.0 eV），观察冻结窗口对插值导带的影响。
4. 用 SCDM（`auto_projections` + `scdm_proj`，`scdm_entanglement='isolated'`）重做价带 Wannier 化，与键中心投影比较。
5. 打开 `write_unk` 与 `wannier_plot`，用 VESTA 画出一个键中心 WF。
6. 把第 14 章的 HSE 全网格计算接到本章流程，画出同一路径上 PBE 与 HSE 的插值能带，读出两者的间接带隙。
7. 用 `si_hr.dat` 统计 |H(R)| 随 |R| 的衰减，判断 4×4×4 网格是否足够。

---

## 参考资料

**官方文档**

- [Wannier90 文档](https://wannier90.readthedocs.io/) — 用户指南（Parameters、Projections、Post-processing）与 tutorials
- [Wannier90 support 页面](https://wannier.org/support/) — 教程、邮件列表与学校资料入口
- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — `nosym`、`noinv`、`nbnd`、`diago_full_acc`
- QE 源码 `PP/src/pw2wannier90.f90`（[GitHub 镜像](https://github.com/QEF/q-e)）— `&inputpp` 全部变量与默认值（本章表格据 QE 7.6 版本整理）
- [QE 7.6 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.6) / [release-notes（develop）](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes) — pw2wannier90 修复与 8.0 中 `wan_mode='library'` 预告

**课程与教程（只链接）**

- Wannier90 官方 tutorials（[文档首页左栏](https://wannier90.readthedocs.io/)）：tutorial 5（金刚石价带）、6（Cu 费米面）、10（石墨）、11（Si 价带与导带）、17（Fe SOC）、27（Si SCDM）、28（cube 格式画 WF）、34（石墨烯可投影性解纠缠）
- [QuantumESPRESSO-school-2023](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023) — EPFL 高级学校，含 Wannier90 上机
- [MaX e-School 2021](https://indico.ictp.it/event/9616/) — 两周课程，含 Wannier 函数讲座与视频
- [pranabdas/espresso](https://pranabdas.github.io/espresso/) — 从安装到 EPW 的系统教程，含 Wannier90
- [PARADIM Summer School 2018（Cornell）](https://www.paradim.org/cu_ss_2018) — Giustino 主讲的 DFT 课程（能带与 Wannier 背景），旧版 QE

**相关 GitHub 仓库（只链接）**

- [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90) — Wannier90 源码、`utility/kmesh.pl`、tutorials
- [quanshengwu/wannier_tools](https://github.com/quanshengwu/wannier_tools) — 表面态、Wilson loop、Weyl 点
- [materialscloud-org/seekpath](https://github.com/materialscloud-org/seekpath) — 标准原胞与 k 路径
- [aiidateam/aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso) — 自动化工作流（含 Wannier 相关插件生态）

**关键文献**

- N. Marzari, D. Vanderbilt, *Phys. Rev. B* 56, 12847 (1997)（最大局域 Wannier 函数）. DOI: [10.1103/PhysRevB.56.12847](https://doi.org/10.1103/PhysRevB.56.12847)
- I. Souza, N. Marzari, D. Vanderbilt, *Phys. Rev. B* 65, 035109 (2001)（解纠缠）. DOI: [10.1103/PhysRevB.65.035109](https://doi.org/10.1103/PhysRevB.65.035109)
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, D. Vanderbilt, *Rev. Mod. Phys.* 84, 1419 (2012)（综述）. DOI: [10.1103/RevModPhys.84.1419](https://doi.org/10.1103/RevModPhys.84.1419)
- G. Pizzi et al., *J. Phys.: Condens. Matter* 32, 165902 (2020)（Wannier90 v3）. DOI: [10.1088/1361-648X/ab51ff](https://doi.org/10.1088/1361-648X/ab51ff)
- A. A. Mostofi et al., *Comput. Phys. Commun.* 185, 2309 (2014)（Wannier90 v2）. DOI: [10.1016/j.cpc.2014.05.003](https://doi.org/10.1016/j.cpc.2014.05.003)
- A. Damle, L. Lin, L. Ying, *J. Chem. Theory Comput.* 11, 1463 (2015)（SCDM）. DOI: [10.1021/ct500985f](https://doi.org/10.1021/ct500985f)
- V. Vitale et al., *npj Comput. Mater.* 6, 66 (2020)（高通量自动 Wannier 化）. DOI: [10.1038/s41524-020-0312-y](https://doi.org/10.1038/s41524-020-0312-y)
- Q. Wu, S. Zhang, H.-F. Song, M. Troyer, A. A. Soluyanov, *Comput. Phys. Commun.* 224, 405 (2018)（WannierTools）. DOI: [10.1016/j.cpc.2017.09.033](https://doi.org/10.1016/j.cpc.2017.09.033)

---

[← 上一章](14-hybrid-metagga.md) | [返回目录](../README.md) | [下一章 →](16-elph-epw.md)
