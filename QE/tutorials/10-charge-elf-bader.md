# 第10章 电荷密度、差分电荷、ELF、STM 与 Bader 分析（Charge Density, Difference, ELF, STM & Bader）

本章讲如何把 pw.x 的自洽结果变成可以看、可以积分的实空间数据。主要工具是后处理程序 pp.x：用 `&INPUTPP` 的 `plot_num` 选物理量（电荷密度 0、自旋密度 6、单个波函数 7、ELF 8、STM 5、PAW 全电子密度 21 等），再用 `&PLOT` 的 `iflag`/`output_format` 输出 Gaussian cube 或 XCrySDen xsf 文件。在此基础上依次完成差分电荷密度（charge density difference，同一 FFT 网格上的 ρ(AB) − ρ(A) − ρ(B)）、电子局域函数（ELF）、Tersoff–Hamann STM 模拟、Henkelman `bader` 程序的 Bader 电荷，以及 projwfc.x 的 Löwdin 电荷和 LOBSTER COHP。

> [!NOTE]
> 本章全部输入在示例目录 [10-charge-elf-bader](../examples/10-charge-elf-bader/README.md)。关键词均对照 QE 7.5/7.6 的 [INPUT_PP](https://www.quantum-espresso.org/Doc/INPUT_PP.html) 核实。仓库不提交 UPF 文件；文中输出片段均为**示意输出，数值仅供参考**，不是真实计算结果。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 pp.x 的"两步"结构：第一步从 `outdir/prefix.save` 读数据、算出物理量并写入中间文件 `filplot`；第二步把一个或多个中间文件按权重组合后写成 1D/2D/3D 图形文件；
- 按需选择 `plot_num`，导出 cube 文件并在 VESTA/XCrySDen 中画等值面和截面；
- 正确地做差分电荷：三次计算共用晶胞、截断能和原子坐标，保证 FFT 网格一致，再用 `nfile`/`weight` 在 pp.x 内部相减；
- 理解 ELF 在 US/PAW 下的局限，知道什么时候换 NC 赝势；
- 用 PAW 全电子密度（`plot_num=21`）作为参考，跑 Henkelman `bader` 程序得到 Bader 电荷；
- 用 `plot_num=5` 模拟恒流 STM 图像，理解 `sample_bias` 的单位和符号；
- 从 projwfc.x 输出中读取 Löwdin 电荷和 spilling parameter，并知道各种"原子电荷"之间不能直接比较。

| 想做的事 | 工具与关键设置 | 本章小节 | 示例文件 |
|---|---|---|---|
| 画价电子密度等值面、截面、沿键曲线 | pp.x `plot_num=0`，`iflag=3/2/1` | §4 | `si_charge.pp.in`、`si_charge_line.pp.in` |
| 吸附、成键引起的电荷转移 | 三次 SCF + `nfile=3`、`weight` | §5 | `co_on_cu_diff/diff_plot.pp.in` |
| 孤对电子、共价键、金属键的可视化 | `plot_num=8`（ELF），NC 赝势 | §6 | `h2o_elf.pp.in` |
| 原子电荷、氧化态趋势 | `plot_num=21` + `bader -ref` | §7 | `co_on_cu_diff/ab_ae.pp.in` |
| 模拟 STM 图像 | `plot_num=5`，`sample_bias`（Ry） | §8 | `co_on_cu_diff/stm.pp.in` |
| 快速原子电荷、基组完备性 | projwfc.x 的 Löwdin 电荷 | §9 | `co_on_cu_diff/ab.projwfc.in` |
| 成键/反键分析 | LOBSTER（外部程序） | §9.3 | — |

**前置知识**：pw.x 输入结构（[第2章](02-input-anatomy.md)）、赝势类型（[第3章](03-pseudopotentials.md)）、slab 建模与真空层（[第9章](09-vdw-surface-workfunction.md)）。功函数用到的 `plot_num=11` 和 average.x 在第9章已经讲过，本章不重复。

### 与 VASP 的对应关系

从 VASP 转过来的读者可以先看这张表（概念对应，文件格式和归一化方式不同）：

| VASP | QE | 说明 |
|---|---|---|
| `CHGCAR`（自动写出） | pp.x `plot_num=0` | QE 不自动写体数据文件，必须用 pp.x 从 `prefix.save` 中提取 |
| `AECCAR0` + `AECCAR2`（`LAECHG=.TRUE.`） | pp.x `plot_num=21`（PAW） | 全电子（价 + 芯）密度，Bader 的参考密度 |
| `ELFCAR`（`LELF=.TRUE.`） | pp.x `plot_num=8` | QE 对 US/PAW 的 ELF 不完整，见 §6 |
| `PARCHG`（`LPARD`、`IBAND`、`EINT`） | pp.x `plot_num=7`（指定 k 点与能带）/ `plot_num=23`（能量窗口内的电荷）/ `plot_num=5`（STM） | 见 §8 |
| `LOCPOT` | pp.x `plot_num=1`（总势）/ `plot_num=11`（V_bare + V_H） | 功函数见第9章 |
| 自旋密度（`CHGCAR` 第二块） | pp.x `plot_num=6` | ρ↑ − ρ↓ |
| `bader CHGCAR -ref CHGCAR_sum` | `bader rho.cube -ref ae.cube` | 同一个 Henkelman 程序，读 cube 格式 |

VASP 侧的完整流程见 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md)；CP2K 的 cube 输出与电荷分析见 [CP2K 第11章](../../CP2K/tutorials/11-cube-charge-workfunction.md)。

---

## 2. 原理速览

### 2.1 pp.x 输出的"电荷密度"是什么

平面波 + 赝势计算中，pw.x 保存的是**价电子的赝电荷密度** ρ̃(r)：

- **NC 赝势**：ρ̃ = Σ f_i |ψ_i|²，在芯区内是平滑的赝密度，不等于真实的全电子密度；
- **US 赝势**：ρ = Σ f_i |ψ_i|² + 增强项（augmentation charges），pp.x `plot_num=0` 输出的是包含增强项的价电子密度；
- **PAW**：除了上面的价电子密度，还可以重建全电子密度。pp.x 提供 `plot_num=17`（全电子**价**电子密度）和 `plot_num=21`（**价 + 芯**全电子密度），官方说明这两项只适用于 PAW，并且需要非常密的实空间网格。

INPUT_PP 文档还明确了两点：pp.x 输出的量均为 Rydberg 原子单位（除非另行说明）；**所有电荷密度的积分等于电子数**（不是总电荷）。所以 cube 中的密度单位是 e/bohr³，乘以体素体积求和应得到价电子数，这是检查文件是否正确的第一步。

### 2.2 差分电荷密度

把体系分成 A（如 slab）和 B（如吸附分子），定义

$$
\Delta\rho(\mathbf r) = \rho_{AB}(\mathbf r) - \rho_{A}(\mathbf r) - \rho_{B}(\mathbf r)
$$

其中 ρ_A、ρ_B 必须是 **A、B 在复合体系 AB 中的原子位置上**分别计算得到的密度。Δρ > 0 的区域电子聚集，Δρ < 0 的区域电子流失。沿表面法向做平面平均可以得到

$$
\Delta\bar\rho(z) = \frac{1}{S}\iint \Delta\rho\,dx\,dy,\qquad
\Delta Q(z) = \int_0^z \Delta\bar\rho(z')\,S\,dz'
$$

ΔQ(z) 在两个子体系分界处的值可以粗略估计电荷转移量（这个数依赖分界面的选取，只适合比较趋势）。

pp.x 还提供 `plot_num=9`：ρ − Σ ρ_atom（总密度减去孤立原子密度的叠加），也叫"变形密度"（deformation density）。它只需要一次计算，适合看共价键的形成，但它回答的是"原子 → 固体"的电荷重排，与上面"分子 + 表面 → 吸附体系"的差分电荷是两个不同的问题。

### 2.3 电子局域函数（ELF）

ELF 由 Becke 和 Edgecombe 提出，pp.x 源码 `PP/src/elf.f90` 中的实现是

$$
\mathrm{ELF} = \frac{1}{1+d^2},\qquad
d = \frac{t(\mathbf r) - t_W(\mathbf r)}{t_{TF}(\mathbf r)}
$$

- $t = \frac{1}{2}\sum_i f_i|\nabla\psi_i|^2$：KS 动能密度；
- $t_W = \frac{1}{8}|\nabla\rho|^2/\rho$：von Weizsäcker 动能密度（玻色子极限）；
- $t_{TF} = \frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}$：均匀电子气的 Thomas–Fermi 动能密度。

ELF 介于 0 和 1 之间：接近 1 表示电子高度局域（共价键、孤对电子、芯层），约 0.5 表示类似均匀电子气（简单金属），接近 0 表示电子稀少或离域。

**局限**：源码注释写明，对 US 赝势，动能密度只包含波函数平滑部分的贡献，增强区的贡献缺失；pp.x 运行时会打印 `elf/kin + US not fully implemented`。因此做 ELF 分析时应使用 **NC 赝势**（PseudoDojo 或 SG15 ONCV）。另外 `plot_num=8` 不支持非共线计算（报错 `elf+noncolin not yet implemented`）。

### 2.4 Bader 电荷

Bader（QTAIM）用电荷密度的**零通量面**（∇ρ·n = 0）把空间划分为原子盆（basin），对每个盆积分得到原子电荷。Henkelman 组的 `bader` 程序在网格上实现这一划分，读取 VASP CHGCAR 或 Gaussian cube 格式。

赝势计算的价电子密度在核附近没有极大值，零通量面的位置会受到影响，所以推荐的做法是：**用包含芯电子的全电子密度确定盆的边界（`-ref` 参考密度），再对价电子密度在这些盆中积分**。在 QE 中，参考密度来自 PAW 的 `plot_num=21`。

原子的净 Bader 电荷 = 该元素赝势的价电子数 Z_val（UPF 文件头中的 `z_valence`）− 该原子盆内的价电子数。

### 2.5 Tersoff–Hamann STM

Tersoff–Hamann 近似把隧穿电流写成针尖位置处、费米能附近一个能量窗口内的局域态密度（LDOS）积分：

$$
I(\mathbf r_0, V) \propto \int_{E_F}^{E_F+eV} \rho(\mathbf r_0, E)\,dE
$$

pp.x `plot_num=5` 计算的就是这个积分：源码 `PP/src/stm.f90` 说明，**费米能 E_F 到 E_F + `sample_bias` 之间的态**都被计入。`sample_bias` 的单位是 **Ry**：负值对应样品的占据态（电子从样品隧穿到针尖），正值对应空态。源码还要求 slab 的表面法向沿第三个晶格矢量，并提醒两侧对称的 slab 可能处理不正确。

---

## 3. pp.x 的输入结构

### 3.1 两个 namelist

```fortran
&INPUTPP              ! step 1: read prefix.save, compute the quantity, write filplot
  prefix   = 'si'
  outdir   = './tmp'
  filplot  = 'si.rho.pp'
  plot_num = 0
/
&PLOT                 ! step 2: read filepp(i), combine with weight(i), write fileout
  nfile         = 1
  filepp(1)     = 'si.rho.pp'
  iflag         = 3
  output_format = 6
  fileout       = 'si_rho.cube'
/
```

INPUT_PP 规定：只想做第二步时，把 `&INPUTPP` **留空**（只写 `&INPUTPP` 和 `/`）；只想做第一步时，**不写** `&PLOT`。中间文件 `filplot` 是格式化文本，可以拷到别的机器上再画图。差分电荷正是利用这一点：先分别提取三个 `filplot`，最后一次性组合。

### 3.2 `&INPUTPP` 关键词

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `prefix` | 与 pw.x 相同 | 读取 `outdir/prefix.save` |
| `outdir` | 与 pw.x 相同 | 默认取环境变量 `ESPRESSO_TMPDIR`，否则为当前目录 |
| `filplot` | 有意义的文件名，如 `ab.rho.pp` | 第一步的中间文件，默认 `'tmp.pp'` |
| `plot_num` | 见下表 | 默认 `-1`（跳过第一步） |
| `spin_component` | `0` | `plot_num=0/9/17` 时：0 = 总密度，1 = 自旋向上，2 = 自旋向下 |
| `sample_bias` | 如 `-0.0735`（约 −1 eV） | `plot_num=5` 时必需，单位 **Ry** |
| `kpoint`、`kband`、`lsign` | 按需 | `plot_num=7` 时选择 k 点和能带；`lsign=.true.` 且为 Γ 点时输出带符号的 \|ψ\|² |
| `emin`、`emax` | eV | `plot_num=3/10/23` 的能量窗口 |

常用的 `plot_num`（完整列表见 INPUT_PP）：

| `plot_num` | 物理量 | 备注 |
|---|---|---|
| 0 | 价电子（赝）电荷密度 | 最常用 |
| 1 | 总势 V_bare + V_H + V_xc | |
| 5 | STM 图像（Tersoff–Hamann） | 需 `sample_bias` |
| 6 | 自旋极化 ρ↑ − ρ↓ | 需 `nspin=2` |
| 7 | 选定波函数对密度的贡献 | NC 时即 \|ψ\|² |
| 8 | ELF | US/PAW 不完整，不支持非共线 |
| 9 | ρ − 原子密度叠加 | 变形密度 |
| 10 | ILDOS（`emin` 到 `emax`） | |
| 11 | V_bare + V_H | 功函数（第9章） |
| 17 | 全电子价电子密度 | 仅 PAW，需要很密的网格 |
| 19 / 20 | 约化密度梯度 / sign(λ₂)ρ | NCI 分析 |
| 21 | 全电子（价 + 芯）密度 | 仅 PAW，需要很密的网格 |
| 22 | 动能密度 | |
| 23 | `emin` 与 `emax` 之间的态的电荷密度 | 类似 VASP 的 `EINT` |
| 123 | DORI | |

### 3.3 `&PLOT` 关键词

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `nfile` | 1；差分电荷用 3 | 读取的中间文件数 |
| `filepp(i)` | 第一步的 `filplot` | 原子坐标从**第一个**文件读取；原子数不同时，第一个文件的原子数必须最多 |
| `weight(i)` | 1.0 / −1.0 | 输出 Σ weight(i) × ρ(i) |
| `iflag` | 3 | 0 = 球平均 1D，1 = 1D 线，2 = 2D 面，3 = 3D，4 = 球面 2D |
| `output_format` | 6（cube）或 5（xsf） | 0 = gnuplot 1D，2 = plotrho 2D，3 = XCrySDen 2D/自定义 3D 区域，5 = XCrySDen 3D（整个 FFT 网格），6 = Gaussian cube，7 = gnuplot 2D |
| `fileout` | 如 `si_rho.cube` | 默认写到标准输出 |
| `x0`、`e1`、`e2`、`e3` | alat 单位 | 1D/2D 必需；3D 可选，给出时必须相互正交 |
| `nx`、`ny`、`nz` | 如 201 | 线/面/体上的点数 |

> [!TIP]
> 3D 输出优先用 `output_format=6`（cube）。它被 VESTA、VMD、Multiwfn、`bader` 以及本仓库的 `qe_cube_diff.py` 直接读取。`output_format=5` 写出整个 FFT 网格的 xsf，适合 XCrySDen，对任意晶轴都很快。

---

## 4. 实战一：Si 的价电子密度

### 4.1 SCF

示例 `si.scf.in` 是标准的金刚石 Si（`ibrav=2`，`celldm(1)=10.26` Bohr，SSSP US 赝势 `Si.pbe-n-rrkjus_psl.1.0.0.UPF`，8×8×8 k 点）。与第2章相比，只把 `conv_thr` 收紧到 `1.0d-10`，以减少密度中的数值噪声：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.26
  nat       = 2
  ntyp      = 1
  ecutwfc   = 30.0
  ecutrho   = 240.0
/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

### 4.2 导出 3D cube 和沿键的 1D 曲线

```bash
cd QE/examples/10-charge-elf-bader
mpirun -np 4 pw.x -nk 4 -in si.scf.in > si.scf.out
mpirun -np 4 pp.x -in si_charge.pp.in > si_charge.pp.out        # -> si_rho.cube
pp.x -in si_charge_line.pp.in > si_charge_line.pp.out           # -> si_rho_bond.dat
```

`si_charge_line.pp.in` 的 `&INPUTPP` 留空，直接复用上一步的 `si.rho.pp`，只做一维输出：

```fortran
&INPUTPP
/
&PLOT
  nfile         = 1
  filepp(1)     = 'si.rho.pp'
  iflag         = 1
  output_format = 0
  x0(1) =  0.0, x0(2) = 0.0, x0(3) = 0.0
  e1(1) = -0.5, e1(2) = 0.5, e1(3) = 0.5
  nx            = 201
  fileout       = 'si_rho_bond.dat'
/
```

`ibrav=2` 时 QE 的 fcc 基矢为 a/2·(−1,0,1)、a/2·(0,1,1)、a/2·(−1,1,0)，晶体坐标 (1/4,1/4,1/4) 对应笛卡尔坐标 (−1/4,1/4,1/4)·a。所以这条线从第一个 Si 出发，经过键中点到达第二个 Si，再延伸到 (−1/2,1/2,1/2)·a 处的空位。

### 4.3 检查与可视化

```bash
grep -e '^!' -e 'JOB DONE' si.scf.out
head -8 si_rho.cube
python ../../scripts/qe_cube_diff.py info si_rho.cube     # usage: python ../../scripts/qe_cube_diff.py -h
```

cube 文件头的格式（示意输出，数值仅供参考；两行注释的具体文字取决于 QE 版本）：

```text
<comment line 1>
<comment line 2>
    2    0.000000    0.000000    0.000000
   45   -0.113...    0.000000    0.113...
   45    0.000000    0.113...    0.113...
   45   -0.113...    0.113...    0.000000
   14    ...         0.000000    0.000000    0.000000
   14    ...        -2.565000    2.565000    2.565000
```

- 第 3 行：原子数和原点；
- 第 4–6 行：每个方向的网格点数和体素向量（单位 bohr，等于晶格矢量除以网格点数）；
- 原子行：原子序数、（第二列的数值依写出程序而定）、坐标（bohr）。

**成功标准**：密度对整个晶胞的积分等于价电子数（Si 两个原子共 8 个价电子），`qe_cube_diff.py info` 会打印这个积分；VESTA 中等值面在 Si–Si 键中点有明显的电荷积聚（共价键）。

> [!NOTE]
> US/PAW 的 `plot_num=0` 密度在芯区附近可能出现很小的负值（增强电荷的数值误差），画等值面时不必在意；但做 Bader 分析时这会干扰盆的划分，这也是要用全电子参考密度的原因之一（§7）。

---

## 5. 实战二：CO/Cu(111) 的差分电荷密度

### 5.1 三次计算的设计原则

差分电荷的每个网格点都在做减法，所以三次计算在**网格层面**必须完全一致：

1. **同一个晶胞**：`CELL_PARAMETERS` 逐位相同；
2. **相同的 `ecutwfc`/`ecutrho`**：FFT 网格由 `ecutrho` 和晶胞决定，两者相同则网格相同；不放心时在三个输入中显式写相同的 `nr1, nr2, nr3`；
3. **原子位置完全相同**：slab 和 CO 的坐标从（已弛豫的）AB 结构中原样拷出，**不要**分别弛豫；
4. **相同的 k 网格、smearing、泛函和赝势**：孤立 CO 本可以只算 Γ 点，但为了和 AB 一致，这里也用 6×6×1 和同样的 smearing；
5. **自旋**：若 AB 无磁性而某个子体系有磁性（例如吸附的自由基），需要想清楚参考态是什么，并在输出中注明。

三次 SCF 跑完后，第一件事是确认网格一致：

```bash
grep 'Dense  grid' ab.scf.out slab.scf.out co.scf.out
```

```text
ab.scf.out:     Dense  grid:   ...... G-vectors     FFT dimensions: (  48,  48, 216)
slab.scf.out:   Dense  grid:   ...... G-vectors     FFT dimensions: (  48,  48, 216)
co.scf.out:     Dense  grid:   ...... G-vectors     FFT dimensions: (  48,  48, 216)
```

（示意输出，数值仅供参考）三行的 `FFT dimensions` 必须完全一致。

### 5.2 AB 输入

`co_on_cu_diff/ab.scf.in` 是 2×2 Cu(111) 三层 slab（a = 3.63 Å，层间距 a/√3）加上顶位 CO；O 原子到上方周期像中底层 Cu 的距离约 15 Å：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'ab'
  outdir      = './tmp'
  pseudo_dir  = '../pseudo'
/
&SYSTEM
  ibrav       = 0
  nat         = 14
  ntyp        = 3
  ecutwfc     = 90.0
  ecutrho     = 720.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr    = 1.0d-8
  mixing_beta = 0.3
  mixing_mode = 'local-TF'
/
ATOMIC_SPECIES
  Cu  63.546  Cu.paw.z_11.ld1.psl.v1.0.0-low.upf
  C   12.011  C.pbe-n-kjpaw_psl.1.0.0.UPF
  O   15.999  O.pbe-n-kjpaw_psl.0.1.UPF
CELL_PARAMETERS angstrom
  5.13360  0.00000   0.00000
  2.56680  4.44582   0.00000
  0.00000  0.00000  22.00000
ATOMIC_POSITIONS angstrom
  Cu  0.00000  0.00000  2.00000
  Cu  1.28340  2.22291  2.00000
  Cu  2.56680  0.00000  2.00000
  Cu  3.85020  2.22291  2.00000
  Cu  1.28340  0.74097  4.09578
  Cu  2.56680  2.96388  4.09578
  Cu  3.85020  0.74097  4.09578
  Cu  5.13360  2.96388  4.09578
  Cu  2.56680  1.48194  6.19156
  Cu  3.85020  3.70485  6.19156
  Cu  5.13360  1.48194  6.19156
  Cu  6.41699  3.70485  6.19156
  C   2.56680  1.48194  8.04156
  O   2.56680  1.48194  9.19156
K_POINTS automatic
  6 6 1 0 0 0
```

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `ecutwfc` / `ecutrho` | 90 / 720 Ry | 三种元素 SSSP 1.3 提示值中的最大值（Cu 决定） |
| `occupations` / `smearing` / `degauss` | `'smearing'` / `'mv'` / 0.02 Ry | 金属 slab；三个输入保持一致 |
| `mixing_mode` | `'local-TF'` | 有真空的非均匀体系更稳定 |
| `K_POINTS` | 6 6 1 0 0 0 | 演示值，需按第4章测试；z 方向只取 1 |

> [!WARNING]
> 示例中的坐标是按理想 fcc 堆垛和经验 C–Cu、C–O 键长摆出的**未弛豫初猜**。正式分析前先用第5章的 `relax`（底层 Cu 用 `if_pos` 固定，见第9章）优化 AB，再把**弛豫后的全部坐标**同步到 `slab.scf.in` 和 `co.scf.in`。

`slab.scf.in` 和 `co.scf.in` 只是从 AB 中删去对应原子（并相应修改 `nat`、`ntyp`、`ATOMIC_SPECIES`），其余逐字相同。

### 5.3 提取三个密度并在 pp.x 中相减

先对三次计算分别只做第一步（不写 `&PLOT`）：

```fortran
&INPUTPP
  prefix   = 'ab'
  outdir   = './tmp'
  filplot  = 'ab.rho.pp'
  plot_num = 0
/
```

（`rho_slab.pp.in`、`rho_co.pp.in` 只改 `prefix` 和 `filplot`。）然后用 `nfile=3` 在第二步做线性组合：

```fortran
&INPUTPP
/
&PLOT
  nfile     = 3
  filepp(1) = 'ab.rho.pp',   weight(1) =  1.0
  filepp(2) = 'slab.rho.pp', weight(2) = -1.0
  filepp(3) = 'co.rho.pp',   weight(3) = -1.0
  iflag         = 3
  output_format = 6
  fileout       = 'drho.cube'
/
```

完整的命令序列写在 `run_diff.sh` 里（默认只打印命令，`RUN=1` 才执行）：

```bash
cd QE/examples/10-charge-elf-bader/co_on_cu_diff
bash run_diff.sh                       # dry run: review the commands first
RUN=1 PW="mpirun -np 16 pw.x -nk 4" PP="mpirun -np 16 pp.x" bash run_diff.sh
```

另一种做法是三个密度各自写成 cube，再在 Python 中相减：

```bash
python ../../../scripts/qe_cube_diff.py -h        # 查看 diff / planar / info 子命令
```

两种方法结果相同，前提都是网格一致；pp.x 内部组合省去了三个大 cube 文件的读写，`qe_cube_diff.py` 的好处是可以顺便做平面平均和 ΔQ(z)。

### 5.4 可视化与定量

- **VESTA**：打开 `drho.cube`，在 Properties → Isosurfaces 中设置正负两个等值面（例如 ±0.002 e/bohr³，具体数值要按体系调整），用两种颜色分别表示电子聚集和流失；
- **截面**：VESTA 的 2D 截面工具选取穿过 Cu–C–O 轴的竖直平面，比等值面更容易看出 σ 给予/π 反馈的特征形状；
- **平面平均**：用 `qe_cube_diff.py planar`（在 `co_on_cu_diff/` 目录下用 `python ../../../scripts/qe_cube_diff.py -h` 查看用法）沿 z 方向得到 Δρ̄(z) 和 ΔQ(z)。

**成功标准**：

1. `drho.cube` 的全晶胞积分接近 0（电子数守恒：N(AB) = N(slab) + N(CO)）；
2. 远离吸附位点的 Cu 原子周围 Δρ 很小，否则说明 slab 的坐标与 AB 不一致；
3. CO 与表面 Cu 之间出现成对的正负区域，而不是噪声状的散点（散点通常意味着 `conv_thr` 太松或网格不一致）。

---

## 6. 实战三：H2O 分子的 ELF

### 6.1 为什么换 NC 赝势

§2.3 已经说明：pp.x 的 ELF 对 US/PAW 只包含平滑部分的动能密度。所以本例改用 SG15 ONCV 模守恒赝势（`O_ONCV_PBE-1.2.upf`、`H_ONCV_PBE-1.2.upf`），`ecutrho` 取 4 × `ecutwfc`：

```fortran
&CONTROL
  calculation = 'scf'
  prefix      = 'h2o'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 1
  A         = 12.0
  nat       = 3
  ntyp      = 2
  ecutwfc   = 80.0
  ecutrho   = 320.0
  assume_isolated = 'mt'
/
&ELECTRONS
  conv_thr = 1.0d-9
/
ATOMIC_SPECIES
  O  15.999  O_ONCV_PBE-1.2.upf
  H   1.008  H_ONCV_PBE-1.2.upf
ATOMIC_POSITIONS angstrom
  O   6.0000  6.0000  6.0000
  H   6.7569  6.5860  6.0000
  H   5.2431  6.5860  6.0000
K_POINTS gamma
```

`ecutwfc = 80` Ry 是 NC 赝势的常见量级，请按第4章自行测试。`assume_isolated='mt'`（Martyna–Tuckerman）消除周期像之间的静电相互作用，第2章有介绍。

### 6.2 ELF 输入与运行

```fortran
&INPUTPP
  prefix   = 'h2o'
  outdir   = './tmp'
  filplot  = 'h2o.elf.pp'
  plot_num = 8
/
&PLOT
  nfile         = 1
  filepp(1)     = 'h2o.elf.pp'
  iflag         = 3
  output_format = 6
  fileout       = 'h2o_elf.cube'
/
```

```bash
mpirun -np 4 pw.x -in h2o.scf.in > h2o.scf.out
mpirun -np 4 pp.x -in h2o_elf.pp.in > h2o_elf.pp.out
```

> [!IMPORTANT]
> 当前版本的 pp.x 每次计算 ELF 时都会打印 `elf/kin + US not fully implemented`（源码中无条件调用，与赝势类型无关）。用 NC 赝势时可以忽略；用 US/PAW 时这条提示意味着结果不完整。

> [!NOTE]
> `K_POINTS gamma` 会启用 Γ 点技巧（只存一半平面波系数）。pp.x 的 ELF 和密度输出支持这种情况；但第12章的 ph.x **不能**读取 Γ 点技巧的数据（报错 `cannot start from pw.x data file using Gamma-point tricks`），做分子振动时要改用 `K_POINTS automatic 1 1 1 0 0 0`。

### 6.3 解读

- 在 VESTA 中把等值面设在 0.8 左右，可以看到 O 原子两侧的两个孤对电子"兔耳"和两个 O–H 成键区域；
- 用 2D 截面（分子平面）显示 0–1 的颜色图，比等值面更直观；
- ELF 是**拓扑**描述，不是电荷：它告诉你电子对在哪里局域，但不能直接给出转移了多少电子；
- 金属表面的 ELF 在真空区和间隙区通常在 0.5 附近起伏，解释时要与均匀电子气的参考值对比。

---

## 7. 实战四：Bader 电荷

### 7.1 获取 bader 程序

Henkelman 组的 Bader 程序（v1.05，2023-08）可以从 <https://theory.cm.utexas.edu/henkelman/code/bader/> 下载 Linux/macOS 可执行文件或 F90 源码。程序自动识别 CHGCAR 和 cube 格式，基本用法是 `bader chargefile`，常用选项：

| 选项 | 说明 |
|---|---|
| `-ref reference_charge` | 用参考密度确定 Bader 盆的边界，再对 chargefile 积分 |
| `-b neargrid / ongrid / weight` | 划分算法；`weight` 为 Yu–Trinkle 权重法（默认值以 `bader -h` 为准） |
| `-vac off / auto / 数值` | 是否把低密度点归入"真空"（1.03 起默认 off） |
| `-p all_atom / sel_atom ...` | 输出每个原子盆的体数据 |

主要输出文件是 `ACF.dat`（每个原子的坐标、`CHARGE`、到零通量面的最小距离 `MIN DIST`、`ATOMIC VOL`），以及 `BCF.dat`、`AtomVolumes.dat`。

### 7.2 QE 侧：价电子密度 + PAW 全电子参考

```fortran
&INPUTPP
  prefix   = 'ab'
  outdir   = './tmp'
  filplot  = 'ab.ae.pp'
  plot_num = 21
/
&PLOT
  nfile = 1, filepp(1) = 'ab.ae.pp'
  iflag = 3, output_format = 6, fileout = 'ab_ae.cube'
/
```

```bash
mpirun -np 16 pp.x -in rho_ab.pp.in > rho_ab.pp.out    # ab.rho.pp (plot_num=0)
mpirun -np 16 pp.x -in ab_rho.pp.in > ab_rho.pp.out    # ab.rho.pp -> ab_rho.cube
mpirun -np 16 pp.x -in ab_ae.pp.in  > ab_ae.pp.out     # plot_num=21 -> ab_ae.cube
bader ab_rho.cube -ref ab_ae.cube
```

`ab_rho.cube`（价电子）是被积分的密度，`ab_ae.cube`（价 + 芯）只用来确定零通量面。这与 VASP 中 `bader CHGCAR -ref CHGCAR_sum` 的思路相同。

### 7.3 读 ACF.dat 并换算净电荷

```text
    #         X           Y           Z       CHARGE      MIN DIST   ATOMIC VOL
 --------------------------------------------------------------------------------
    1         ...         ...         ...     ...          ...         ...
   ...
   13         ...         ...         ...     ...          ...         ...
   14         ...         ...         ...     ...          ...         ...
 --------------------------------------------------------------------------------
    VACUUM CHARGE:               0.0000
    VACUUM VOLUME:               0.0000
    NUMBER OF ELECTRONS:       ...
```

（示意输出，数值仅供参考；第 13、14 行对应 C 和 O。）

净电荷 q = Z_val − CHARGE。Z_val 从 UPF 头中读取：

```bash
grep -i 'z_valence' ../pseudo/Cu.paw.z_11.ld1.psl.v1.0.0-low.upf \
                    ../pseudo/C.pbe-n-kjpaw_psl.1.0.0.UPF ../pseudo/O.pbe-n-kjpaw_psl.0.1.UPF
grep 'number of electrons' ab.scf.out
```

`NUMBER OF ELECTRONS` 应与 `ab.scf.out` 中的 `number of electrons` 一致。如果差得多，通常是 cube 网格太粗或文件被截断。

### 7.4 网格收敛

INPUT_PP 明确说 `plot_num=17/21` "需要非常密的实空间网格"。实际操作：

1. 用默认 `ecutrho` 得到一组 Bader 电荷；
2. 在 **SCF 输入**中把 `nr1, nr2, nr3` 显式加大（例如各放大 1.5 倍），重新 SCF + pp.x + bader；
3. 关心的原子电荷变化小于 0.01–0.02 e 时认为收敛。

同时检查 `ACF.dat` 的 `MIN DIST`：如果它比赝势的截断半径还小，说明零通量面切进了芯区，结果需要谨慎对待。

> [!TIP]
> 只有 NC 赝势时（比如 §6 的 H2O），没有 `plot_num=21` 可用，可以直接 `bader h2o_rho.cube`。这时盆边界由价电子赝密度决定，电荷的绝对值偏差较大，但同一组设置下的**相对变化**（比如吸附前后）仍有参考价值。

---

## 8. 实战五：Tersoff–Hamann STM 图像

### 8.1 输入

```fortran
&INPUTPP
  prefix      = 'ab'
  outdir      = './tmp'
  filplot     = 'ab.stm.pp'
  plot_num    = 5
  sample_bias = -0.0735
/
&PLOT
  nfile = 1, filepp(1) = 'ab.stm.pp'
  iflag = 3, output_format = 6, fileout = 'ab_stm.cube'
/
```

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `plot_num` | 5 | STM（Tersoff–Hamann） |
| `sample_bias` | −0.0735 Ry（约 −1.0 eV） | 1 Ry = 13.6057 eV；负值对应占据态，正值对应空态 |

使用前的检查清单：

- **表面法向沿第三个晶格矢量**（源码要求）；本例 slab 法向沿 z，满足；
- **两侧对称的 slab** 可能处理不正确（源码注释），STM 示例最好用单面吸附的不对称 slab；
- **绝缘体/半导体**（没有用 smearing）：pp.x 用最低空带和最高占据带确定带隙中点作为 E_F，源码要求 `nbnd` > 占据带数 + 1（即至少两条空带），否则报错 `not enough bands`，需要在 SCF 中加大 `nbnd`；
- k 点要足够密，能量窗口小于 0.5 eV 时尤其如此；
- 真空层要足够厚，针尖高度处（表面上方 3–6 Å）的 LDOS 才不受周期像的影响。

### 8.2 从 3D 数据到 STM 图像

pp.x 给出的是 3D 的积分 LDOS。

- **恒高模式**：在 VESTA 中取表面上方某一高度的水平截面；
- **恒流模式**：画一个很小的等值面（例如 10⁻⁵–10⁻⁶ e/bohr³ 量级，按体系调整），从上方俯视并按高度着色；
- 换不同的 `sample_bias`（±0.5、±1.0、±2.0 eV）可以比较占据态和空态图像的差别。

对 CO/Cu(111)，顶位 CO 的 STM 衬度与偏压和针尖状态密切相关，实验中常见"凹陷"衬度。Tersoff–Hamann 近似假设 s 波针尖，不一定能复现实验衬度，这是方法本身的局限，不是输入写错了。

### 8.3 相关的 `plot_num=7` 与 `plot_num=23`

- `plot_num=7` + `kpoint`、`kband`：画某个 k 点上某条能带的 |ψ|²，例如 Γ 点的 HOMO/LUMO（分子）或表面态；`lsign=.true.` 且为 Γ 点时保留波函数的符号，便于看节面；
- `plot_num=23` + `emin`、`emax`（eV）：能量窗口内所有态的电荷密度，相当于 VASP 的 `LPARD` + `EINT`。

---

## 9. Löwdin 电荷与 COHP

### 9.1 projwfc.x 的 Löwdin 电荷

第6章用 projwfc.x 算 PDOS 时，输出末尾会附带 Löwdin 布居分析。本例单独给出一个最简输入：

```fortran
&PROJWFC
  prefix  = 'ab'
  outdir  = './tmp'
  filpdos = 'ab'
/
```

```bash
mpirun -np 16 projwfc.x -in ab.projwfc.in > ab.projwfc.out
grep -A20 'Lowdin Charges' ab.projwfc.out
grep 'Spilling Parameter' ab.projwfc.out
```

```text
Lowdin Charges:

     Atom #   1: total charge =  10.9...., s =  0.6..., p =  0.3..., d =  9.8...,
     ...
     Atom #  13: total charge =   3.7...., s =  1.2..., p =  2.4...,
     Atom #  14: total charge =   6.3...., s =  1.7..., p =  4.5...,
     Spilling Parameter:   0.0...
```

（示意输出，数值仅供参考。）

- `total charge` 是投影到正交化原子轨道上的电子数；净电荷 = Z_val − total charge；
- **Spilling Parameter** 衡量原子轨道基组能多好地展开 KS 态：越小越可信，较大（例如百分之几以上）时 Löwdin 电荷只能定性参考；
- Löwdin 电荷依赖赝势中自带的原子波函数（`PP_PSWFC`），不同赝势之间不可比。

### 9.2 各种"原子电荷"的比较

| 方法 | QE 中怎么得到 | 优点 | 局限 |
|---|---|---|---|
| Bader | pp.x `plot_num=0/21` + `bader` | 只依赖密度，与基组无关；可以和 VASP 结果比较趋势 | 需要密网格；赝势密度需全电子参考 |
| Löwdin | projwfc.x | 随 PDOS 一起得到，几乎不额外耗时 | 依赖赝势的原子轨道；spilling 大时不可靠 |
| 差分电荷 ΔQ(z) | §5 + 平面平均 | 直观显示"哪里多了、哪里少了" | 分界面的选择带有主观性 |
| 变形密度 | `plot_num=9` | 一次计算即可 | 参考态是孤立原子，不是分子 + 表面 |

> [!WARNING]
> 不同方法给出的数值可以差得很远，**不要**把 Bader 电荷和 Löwdin 电荷放在一起比较，也不要把它们等同于形式氧化态。论文中应写明方法、赝势、网格，以及做比较的参考体系。

### 9.3 LOBSTER COHP

晶体轨道哈密顿布居（COHP）把能带能量分解到原子对上，用来区分成键态（−COHP > 0）和反键态。[LOBSTER](http://www.cohp.de/) 是常用的实现，它把平面波波函数投影到局域原子基组上，支持包括 QE 在内的若干平面波程序。QE 侧对赝势类型、k 点对称性、输出文件的具体要求随 LOBSTER 版本变化，请以 LOBSTER 手册为准；社区教程 [Waranchit-TOTO/LOBSTER-Quantum-ESPRESSO-Tutorial](https://github.com/Waranchit-TOTO/LOBSTER-Quantum-ESPRESSO-Tutorial) 给出了一套 QE + LOBSTER 的示例流程，可以作为起点。VASP 侧的 COHP 流程见 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md)。

---

## 10. 可视化工具小结

| 工具 | 读取格式 | 适合做 |
|---|---|---|
| [VESTA](https://jp-minerals.org/vesta/en/) | cube、xsf | 等值面、截面、正负双色差分电荷，出版图 |
| [XCrySDen](http://www.xcrysden.org/) | xsf（`output_format=5`）、QE 输入/输出 | 直接打开 pw.x 输入检查结构，查看 3D 网格数据 |
| VMD | cube | 分子体系、轨迹 + 体数据 |
| Python（numpy/matplotlib） | cube | 平面平均、积分、批量出图（`qe_cube_diff.py`） |
| gnuplot | `iflag=1`、`output_format=0/7` | 1D/2D 快速查看 |

通用的可视化与后处理建议见 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md)；用 ASE 读写 cube 文件见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)。

---

## 参数选择与收敛建议

| 项目 | 建议 |
|---|---|
| SCF 精度 | 做密度分析时 `conv_thr` 至少 1e-8（小体系 1e-10），噪声会直接出现在差分电荷中 |
| `ecutrho` 与网格 | 决定 cube 的分辨率；Bader 需要做网格收敛（§7.4）；差分电荷三次计算必须相同 |
| 赝势 | ELF 用 NC；Bader 最好用 PAW + `plot_num=21`；同一组比较中不要混用不同赝势 |
| k 点 | 金属 slab 的密度对 k 点比较敏感；STM 的小能量窗口需要更密的 k 点 |
| 真空层 | STM 和平面平均需要足够的真空（≥ 12–15 Å），并检查偶极修正（第9章） |
| 等值面 | 报告等值面的数值和单位（e/bohr³ 或 e/Å³）；差分电荷用正负对称的等值面 |
| 文件大小 | 大体系的 cube 可达数百 MB；先用 `iflag=2` 的截面检查，再输出 3D |
| 自旋 | `nspin=2` 时可以用 `spin_component=1/2` 分别输出 ρ↑、ρ↓，或用 `plot_num=6` 直接得到自旋密度 |

---

## 常见问题与报错

**1. pp.x 一启动就因读取 `prefix.save` 中的文件失败而停止。**
原因：`prefix`、`outdir` 与 pw.x 不一致，或 pw.x 使用了 `disk_io='nowf'`/`'minimal'`，没有保存 pp.x 需要的数据。解决：统一 `prefix`/`outdir`；需要波函数的量（ELF、STM、`plot_num=7`）要求 SCF 保存了波函数，不要用 `disk_io='nowf'`。

**2. 差分电荷图满屏噪声，或者远离吸附位点的地方也有大块正负区域。**
原因：三次计算的 FFT 网格或原子位置不一致；或者 slab 被单独弛豫过。解决：`grep 'Dense  grid'` 比较三个输出；坐标从 AB 中原样拷贝；必要时显式写相同的 `nr1, nr2, nr3`。

**3. pp.x 组合文件时原子位置不对，或者报原子数相关的错误。**
原因：`&PLOT` 中第一个 `filepp` 不是原子数最多的文件（INPUT_PP：原子坐标从第一个文件读取）。解决：把 AB 的文件放在 `filepp(1)`。

**4. `elf/kin + US not fully implemented`。**
原因：ELF 的动能密度对 US/PAW 只包含平滑部分（这条提示每次都会打印）。解决：使用 NC 赝势时可以忽略；使用 US/PAW 时换成 PseudoDojo 或 SG15 ONCV 赝势重新计算。

**5. `elf+noncolin not yet implemented`。**
原因：非共线（`noncolin=.true.`）计算不支持 `plot_num=8`。解决：只能用共线计算（`nspin=1` 或 2）做 ELF。

**6. STM：`not enough bands`。**
原因：体系没有 smearing（绝缘体/半导体）时，`PP/src/stm.f90` 要求 `nbnd` 大于占据带数 + 1（至少两条空带）。解决：在 SCF 中把 `nbnd` 设得比占据带多（例如多 20%），重新 SCF。

**7. STM 图像看起来是 slab 底面，或者两面叠加在一起。**
原因：表面法向不沿第三个晶格矢量，或 slab 两面对称。解决：调整晶胞让法向沿 c 轴；使用单面吸附的不对称 slab，取 slab 上方区域的截面。

**8. `bader` 给出的电子总数与 pw.x 的 `number of electrons` 相差很多。**
原因：cube 网格太粗、文件被截断，或者把 `plot_num=21` 的全电子密度当成了被积分的 chargefile。解决：被积分的文件用 `plot_num=0`，全电子密度只放在 `-ref` 后面；加密网格后重做。

**9. Bader 电荷随网格明显变化。**
原因：网格没有收敛，PAW 全电子密度在核附近变化很剧烈。解决：按 §7.4 做网格收敛，报告收敛后的数值和使用的网格。

**10. Löwdin 的 `Spilling Parameter` 很大。**
原因：赝势自带的原子轨道不足以描述 KS 态（缺少高角动量轨道或半芯态处理不当）。解决：Löwdin 电荷只做定性参考，定量分析改用 Bader；PDOS 的相对形状通常仍然可用。

**11. cube 文件在 VESTA 中显示的晶胞或原子位置不对。**
原因：使用了 `output_format=3` 加自定义区域，或给定的 `e1/e2/e3` 不正交。解决：3D 数据优先用 `output_format=6` 且不指定 `e1/e2/e3`，输出整个 FFT 网格。

**12. 对 NC/US 赝势设置了 `plot_num=17` 或 `21`。**
原因：INPUT_PP 规定这两个量只适用于 PAW（需要 PAW 的重建信息），NC/US 计算得到的结果没有意义。解决：改用 PAW 赝势，或者只做价电子密度分析（NC 情形见 §7.4 的 TIP）。

---

## 练习

1. 对 Si 分别输出 `plot_num=0` 和 `plot_num=9`，在同一个截面上比较，说明"变形密度"展示了什么。
2. 把 `si_charge_line.pp.in` 的 `nx` 改成 51 和 401，比较 1D 曲线；再把 `ecutrho` 从 240 提高到 360 Ry，观察键中点密度的变化。
3. 用 `qe_cube_diff.py`（在 `co_on_cu_diff/` 目录下用 `python ../../../scripts/qe_cube_diff.py -h` 查看用法）对 `drho.cube` 做 z 方向平面平均，画出 ΔQ(z)，在 Cu 顶层与 C 原子之间取分界面，估计电荷转移的方向。
4. 对 CO/Cu(111) 分别取 `sample_bias` = −2.0、−0.5、+0.5、+2.0 eV（注意换算成 Ry），比较 STM 图像中 CO 处的衬度。
5. 用 PseudoDojo 或 SG15 赝势重算 H2O，并用 US 赝势再算一次 ELF，比较 O 孤对电子区域的差别，体会 §2.3 的局限。
6. 对 CO/Cu(111) 用两组网格（默认、各方向 ×1.5）做 Bader，判断 C、O 的 Bader 电荷是否收敛；再与 Löwdin 电荷对比，写出两者差异的原因。
7. 把 CO 换成 O 原子吸附（自旋极化体系，`nspin=2`），用 `plot_num=6` 画自旋密度，并检查 slab 和孤立 O 原子的参考态是否一致。

---

## 参考资料

**官方文档**

- pp.x 输入说明 INPUT_PP：<https://www.quantum-espresso.org/Doc/INPUT_PP.html>
- projwfc.x 输入说明 INPUT_PROJWFC：<https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html>
- pw.x 输入说明 INPUT_PW：<https://www.quantum-espresso.org/Doc/INPUT_PW.html>
- PWscf 用户手册：<https://www.quantum-espresso.org/Doc/pw_user_guide/>
- QE 源码中的 ELF 与 STM 实现（GitLab 主仓库）：`PP/src/elf.f90`、`PP/src/stm.f90`、`PP/src/punch_plot.f90`，见 <https://gitlab.com/QEF/q-e>
- 官方 PP 示例：<https://gitlab.com/QEF/q-e/-/tree/develop/PP/examples>（每个 example 目录的 README 说明其内容，`WorkFct_example` 对应第9章的功函数）
- Henkelman Bader 程序：<https://theory.cm.utexas.edu/henkelman/code/bader/>
- 赝势：[SSSP 1.3 Efficiency](https://www.materialscloud.org/discover/sssp/table/efficiency)、[PseudoDojo](https://www.pseudo-dojo.org/)、[SG15 ONCV](http://www.quantum-simulation.org/potentials/sg15_oncv/)

**课程与教程（只链接，不复制）**

- Cornell PARADIM 暑期学校 2016（Giustino）：<https://www.paradim.org/cu_ss_2016>。上机部分包含能带与 KS 波函数可视化，基于 QE 5.x 和 LDA 赝势；下载、编译和赝势部分已经过时，原理部分仍然适用。
- Imperial College MSE404-MM：<https://imperialcollegelondon.github.io/MSE404-MM/>，9 个 QE 实验，持续维护，适合配合本章练习 pp.x 的基本用法。
- pranabdas/espresso：<https://pranabdas.github.io/espresso/>，从安装到能带、DOS、声子、Wannier 的全流程（站点目前没有电荷密度专题，pp.x 部分以本章和 INPUT_PP 为准）。
- QE-SSP（《Quantum ESPRESSO Course for Solid-State Physics》配套）：<https://github.com/nguyen-group/QE-SSP>，其中 `gr/rho/` 是石墨烯电荷密度的 pp.x 算例，`bi-gr/elec-field/` 含外加电场下的 pp.x 后处理。
- LOBSTER + QE 社区示例：<https://github.com/Waranchit-TOTO/LOBSTER-Quantum-ESPRESSO-Tutorial>
- 官方 Python 后处理 postqe：<https://github.com/QEF/postqe>

**本仓库**

- 示例：[10-charge-elf-bader](../examples/10-charge-elf-bader/README.md)
- 脚本：`qe_cube_diff.py`（cube 差分、平面平均、积分）、`qe_workfunction.py`（第9章功函数）、`job_qe_slurm.sh`；用 `python ../scripts/<name> -h` 查看用法。
- 对照阅读：[VASP 第16章 电荷与成键分析](../../VASP/tutorials/16-charge-bonding.md)、[CP2K 第11章 cube 文件与电荷分析](../../CP2K/tutorials/11-cube-charge-workfunction.md)

**关键文献**

- A. D. Becke, K. E. Edgecombe, *J. Chem. Phys.* **92**, 5397 (1990)，ELF：<https://doi.org/10.1063/1.458517>
- J. Tersoff, D. R. Hamann, *Phys. Rev. B* **31**, 805 (1985)，STM 理论：<https://doi.org/10.1103/PhysRevB.31.805>
- G. Henkelman, A. Arnaldsson, H. Jónsson, *Comput. Mater. Sci.* **36**, 354 (2006)，网格 Bader 算法：<https://doi.org/10.1016/j.commatsci.2005.04.010>
- E. Sanville, S. D. Kenny, R. Smith, G. Henkelman, *J. Comput. Chem.* **28**, 899 (2007)：<https://doi.org/10.1002/jcc.20575>
- W. Tang, E. Sanville, G. Henkelman, *J. Phys.: Condens. Matter* **21**, 084204 (2009)，无晶格偏置的 Bader 算法：<https://doi.org/10.1088/0953-8984/21/8/084204>
- M. Yu, D. R. Trinkle, *J. Chem. Phys.* **134**, 064111 (2011)，Bader 积分的权重法：<https://doi.org/10.1063/1.3553716>
- S. Maintz, V. L. Deringer, A. L. Tchougréeff, R. Dronskowski, *J. Comput. Chem.* **37**, 1030 (2016)，LOBSTER：<https://doi.org/10.1002/jcc.24300>
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017)，QE 高级功能：<https://doi.org/10.1088/1361-648X/aa8f79>

---

[← 上一章](09-vdw-surface-workfunction.md) | [返回目录](../README.md) | [下一章 →](11-neb.md)
