# 第12章 声子：ph.x、q2r.x 与 matdyn.x（Phonons with DFPT: ph.x, q2r.x, matdyn.x）

本章讲解用 Quantum ESPRESSO 的密度泛函微扰理论（DFPT）计算声子：先用 pw.x 得到高精度基态（`conv_thr` ≤ 1e-12），再用 ph.x 在单个 q 点或均匀 q 网格（`ldisp=.true.`、`nq1/nq2/nq3`）上求动力学矩阵；用 q2r.x 傅里叶变换得到实空间原子间力常数（IFC），用 matdyn.x 插值出任意 q 路径上的声子色散（`q_in_band_form`）和声子态密度（`dos=.true.`），并施加声学求和规则（`asr='crystal'`）。示例覆盖半导体 Si 和金属 Al，并讨论虚频诊断、`recover`/`start_q`/`last_q` 续算与拆分、image 并行，以及与 phonopy 有限位移方法的对比。

> [!NOTE]
> 本章全部输入在示例目录 [12-phonon-dfpt](../examples/12-phonon-dfpt/README.md)。ph.x、q2r.x、matdyn.x 的关键词与默认值对照 QE 7.5/7.6 的 [INPUT_PH](https://www.quantum-espresso.org/Doc/INPUT_PH.html)、[INPUT_Q2R](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html)、[INPUT_MATDYN](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html) 核实。仓库不提交 UPF 文件；文中输出片段均为**示意输出，数值仅供参考**。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 DFPT 与有限位移（冻结声子）方法的区别，以及 QE 中"q 网格 → 力常数 → 插值"三步流程的物理含义；
- 为 ph.x 准备合格的基态：结构充分弛豫、`conv_thr` 足够小、保留对称性、金属用 smearing；
- 写出 Γ 点单 q 计算和 `ldisp` 网格计算的 ph.x 输入，读懂输出中的不可约表示、自洽收敛和频率；
- 用 q2r.x 和 matdyn.x 得到色散和 DOS，正确设置 `asr`、`q_in_band_form` 和 q 点单位；
- 用 `qe_phonon_plot.py` 或 gnuplot 绘图，并诊断虚频的来源；
- 用 `recover`、`start_q`/`last_q`、`-ni` 拆分和续算大型声子计算；
- 知道何时改用 phonopy 有限位移（杂化泛函、meta-GGA、TS-vdW，以及超出 ph.x 支持范围的 DFT+U 等情况）。

| 想做的事 | 方法 | 本章小节 | 示例文件 |
|---|---|---|---|
| Γ 点频率、拉曼/红外活性模的初步检查 | ph.x 单 q 点 | §4 | `si.ph_gamma.in` |
| 完整声子色散与 DOS | `ldisp` + q2r.x + matdyn.x | §5–§6 | `si.ph.in`、`si.q2r.in`、`si.matdyn.in` |
| 金属声子 | SCF 中 smearing，更密 k 点 | §7 | `al.*.in` |
| 判断动力学稳定性、找虚频原因 | 收敛检查 + 对称性检查 | §8 | — |
| 大体系拆分、续算、并行 | `recover`、`start_q`/`last_q`、`start_irr`/`last_irr`、`-ni` | §9 | `job_ph_slurm.sh` |
| 极性材料 LO–TO、介电常数、Born 电荷、IR/Raman | `epsil`、`lraman`、dynmat.x | [第13章](13-dielectric-ir-raman.md) | — |
| 电声耦合、超导 Tc | `electron_phonon`、EPW | [第16章](16-elph-epw.md) | — |

**前置知识**：结构优化（[第5章](05-relax-eos.md)）、收敛测试（[第4章](04-convergence.md)）、能带的 k 路径（[第6章](06-band-dos-pdos.md)）。

### 与 VASP 的对应关系

| VASP + phonopy | QE | 说明 |
|---|---|---|
| 超胞 + `IBRION=8`（DFPT，仅 q=0） | ph.x 在**原胞**中直接算任意 q | QE 的 DFPT 不需要超胞，q 网格对应于等效的超胞大小 |
| 超胞 + `IBRION=6`/phonopy 有限位移 | phonopy `--qe`（pw.x 超胞力） | 两种方法都能在 QE 中使用 |
| phonopy 的 `FORCE_CONSTANTS` | q2r.x 输出的 `flfrc` | 实空间力常数 |
| `band.conf` + `phonopy -p` | matdyn.x + `q_in_band_form` | |
| `mesh.conf`（DOS） | matdyn.x `dos=.true.` | |
| `LEPSILON=.TRUE.` → BORN | ph.x `epsil=.true.`（第13章） | 非解析项修正 |
| phonopy `FC_SYMMETRY` / 对称化 | `asr='crystal'`（q2r.x 的 `zasr`、matdyn.x 的 `asr`） | 声学求和规则 |

VASP 的声子流程见 [VASP 第12章](../../VASP/tutorials/12-phonon.md)；CP2K 的振动分析见 [CP2K 第14章](../../CP2K/tutorials/14-vibration-thermo.md)。

---

## 2. 原理速览

### 2.1 简谐近似与动力学矩阵

把总能量对原子位移 u 展开到二阶，原子间力常数为

$$
C_{\alpha\beta}(\kappa l,\kappa' l') = \frac{\partial^2 E}{\partial u_\alpha(\kappa l)\,\partial u_\beta(\kappa' l')}
$$

其中 κ 标记原胞中的原子，l 标记原胞。对波矢 q 的动力学矩阵和本征方程为

$$
D_{\alpha\beta}(\kappa\kappa',\mathbf q) = \frac{1}{\sqrt{M_\kappa M_{\kappa'}}}\sum_{l'} C_{\alpha\beta}(\kappa 0,\kappa' l')\,e^{i\mathbf q\cdot(\mathbf R_{l'})},\qquad
D(\mathbf q)\,\mathbf e = \omega^2\,\mathbf e
$$

ω² < 0 对应虚频。QE 的约定是把这种频率打印成**负数**（PHonon 用户手册：负频率实际上是虚频，表示动力学矩阵本征值 ω² < 0）。

### 2.2 DFPT：直接在 q 处求二阶导数

DFPT 把原子按波矢 q 的调制位移当作微扰，自洽求解电子密度的线性响应 Δρ（Sternheimer 方程），再由 Δρ 得到二阶导数。它的优点是：

- 每个 q 点的计算都在**原胞**中完成，与 q 的大小无关；
- 不需要选择位移幅度，没有有限差分误差；
- 对称性被充分利用：每个 q 点只需计算若干**不可约表示**（irreducible representation, irrep），ph.x 输出中的 `Representation #` 就是它们。

代价是每个 irrep 要做一次线性响应的自洽循环，收敛阈值由 `tr2_ph` 控制。

### 2.3 从 q 网格到任意 q：q2r.x 与 matdyn.x

在 nq1×nq2×nq3 的均匀 q 网格上算出 D(q) 后：

1. **q2r.x**：傅里叶变换得到实空间力常数 C(R)，R 只覆盖 nq1×nq2×nq3 的"等效超胞"；
2. **matdyn.x**：由 C(R) 反变换到任意 q，得到色散或 DOS。

所以 q 网格的作用相当于有限位移法中的超胞大小：力常数必须在这个范围内衰减到可忽略，插值才准确。半导体和简单金属 4×4×4 常常是一个合理的起点，长程力常数明显的体系（极性材料、Kohn 异常的金属）需要更密的网格。极性材料的长程偶极项由非解析修正单独处理（第13章）。

### 2.4 声学求和规则（ASR）

平移不变性要求 Γ 点的三个声学模频率严格为零。数值计算中，交换关联能在离散的实空间 FFT 网格上计算，平移不变性只对网格步长的平移成立，所以 Γ 点声学模会有小的非零频率。PHonon 用户手册指出它通常小于 10 cm⁻¹，有时可达 100 cm⁻¹；GGA、含真空的体系更明显。解决办法是在后处理中施加 ASR：q2r.x 的 `zasr`（作用于 Born 有效电荷）和 matdyn.x 的 `asr`（作用于力常数），推荐 `'crystal'`。

---

## 3. 准备基态

### 3.1 对 SCF 的要求

| 要求 | 原因 | 做法 |
|---|---|---|
| 结构充分弛豫 | 有残余力时，二阶导数不对应平衡结构，容易出现虚频 | 第5章 `relax`/`vc-relax`，`forc_conv_thr` 收紧到 1e-4 Ry/Bohr 量级 |
| `conv_thr` ≤ 1e-12（小体系） | DFPT 的线性响应以基态为起点，基态噪声会传递到频率 | 大体系按原子数适当放宽，但一般不宽于 1e-10 |
| 不能直接接在 `vc-relax` 后 | ph.x 源码检查：`The phonon code is not working after vc-relax` | 用优化后的晶格参数重新跑一次 `scf` |
| 保留对称性 | 对称性决定 irrep 的划分；接近但不精确的对称位置会导致对称性相关报错 | 尽量用正确的 `ibrav` 和 Wyckoff 位置；`ibrav=0` 时坐标要对称化 |
| 不能用 `K_POINTS gamma` | ph.x 报错 `cannot start from pw.x data file using Gamma-point tricks` | 分子也用 `K_POINTS automatic 1 1 1 0 0 0` |
| 金属用 smearing | 否则 ph.x 提示占据数可能错误 | `occupations='smearing'`，见 §7 |
| 泛函 | ph.x 不支持杂化泛函、meta-GGA、TS-vdW（源码 `phq_readin.f90` 中直接报错） | 这些情况改用 phonopy 有限位移（§10） |

### 3.2 Si 的 SCF 输入

`si.scf.in`：

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
  conv_thr = 1.0d-12
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

`celldm(1)=10.26` Bohr 只是示例值，请换成你自己用同一套赝势和截断能优化得到的晶格常数。赝势 `Si.pbe-n-rrkjus_psl.1.0.0.UPF` 来自 SSSP 1.3 PBE Efficiency（提示截断 30/240 Ry）。SSSP 的验证项中包含声子频率的收敛，但这是针对其测试体系的，你的体系仍应自行检查。

```bash
cd QE/examples/12-phonon-dfpt
mpirun -np 8 pw.x -nk 4 -in si.scf.in > si.scf.out
grep -e '^!' -e 'convergence has been achieved' -e 'JOB DONE' si.scf.out
```

---

## 4. 单个 q 点：Γ 点声子

### 4.1 输入

ph.x 的输入以**一行标题**开始，然后是 `&INPUTPH`，非 `ldisp` 计算时再跟一行 q 点坐标（2π/a 单位）：

```text
Si phonons at Gamma only (single q point)
&INPUTPH
  prefix = 'si'
  outdir = './tmp'
  tr2_ph = 1.0d-16
  fildyn = 'si.dynG'
/
0.0 0.0 0.0
```

```bash
mpirun -np 8 ph.x -nk 4 -in si.ph_gamma.in > si.ph_gamma.out
```

### 4.2 `&INPUTPH` 常用关键词

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `prefix`、`outdir` | `'pwscf'`、`ESPRESSO_TMPDIR` 或 `./` | 与 pw.x 相同 | 读取基态数据 |
| `tr2_ph` | 1e-12 | 1e-14 ~ 1e-16 | 线性响应自洽阈值；越小越准，代价增加有限 |
| `fildyn` | `'matdyn'` | 有意义的名字，如 `'si.dyn'` | 动力学矩阵文件 |
| `ldisp` | `.false.` | 色散时 `.true.` | 在 `nq1×nq2×nq3` 网格上计算 |
| `nq1, nq2, nq3` | 0, 0, 0 | 如 4, 4, 4 | 不带偏移的 Monkhorst–Pack q 网格 |
| `amass(i)` | 从数据文件读取 | 同位素效应时设置 | 原子质量（amu） |
| `epsil` | `.false.` | 绝缘体 Γ 点时 `.true.` | 介电常数与 Born 电荷（第13章）；**金属或 q ≠ 0 时不能设** |
| `recover` | `.false.` | 续算时 `.true.` | 从中断处继续 |
| `alpha_mix(1)` | 0.7 | 难收敛时 0.3 | 线性响应的混合因子 |
| `nmix_ph` | 4 | 8–20 | 混合使用的历史步数，加大可加速收敛，占用更多内存 |
| `max_seconds` | 1e7 | 略小于作业时限 | 到时平稳停止，便于 `recover` |
| `start_q`、`last_q` | 1、q 点总数 | 拆分任务时设置 | 只算一部分 q 点（仅 `ldisp`） |
| `start_irr`、`last_irr` | 1、3·nat | 拆分任务时设置 | 只算一部分 irrep |
| `search_sym` | `.true.` | 默认 | 设为 `.false.` 可关闭模式的对称性分析 |

### 4.3 读输出

ph.x 依次输出：对称性与 irrep 划分、每个 irrep 的自洽迭代、动力学矩阵对角化结果。

```bash
grep -e 'Representation #' -e 'Convergence has been achieved' si.ph_gamma.out
grep -A10 'Diagonalizing the dynamical matrix' si.ph_gamma.out
```

```text
     Diagonalizing the dynamical matrix

     q = (    0.000000000   0.000000000   0.000000000 )

 **************************************************************************
     freq (    1) =      -0.0...  [THz] =      -0.....  [cm-1]
     freq (    2) =      -0.0...  [THz] =      -0.....  [cm-1]
     freq (    3) =      -0.0...  [THz] =      -0.....  [cm-1]
     freq (    4) =      15.....  [THz] =     5.......  [cm-1]
     freq (    5) =      15.....  [THz] =     5.......  [cm-1]
     freq (    6) =      15.....  [THz] =     5.......  [cm-1]
 **************************************************************************
```

（示意输出，数值仅供参考。）

- 前三个是声学模，应接近 0；ph.x 这里**没有**施加 ASR，出现几个 cm⁻¹ 的正值或负值是正常的；
- 后三个是光学模，Si 的光学模在 Γ 点三重简并（Γ₂₅' 表示），频率与实验约 15.5 THz 同量级（以你的计算为准）；
- 每个 irrep 对应一行 `Convergence has been achieved`；如果出现 `No convergence has been achieved`，见常见问题。

想对单 q 点结果施加 ASR，可以用 dynmat.x（`fildyn='si.dynG'`、`asr='crystal'`），第13章会详细使用它。

---

## 5. q 网格：`ldisp` 计算

### 5.1 输入

```text
Si phonons on a 4x4x4 q-point grid
&INPUTPH
  prefix = 'si'
  outdir = './tmp'
  tr2_ph = 1.0d-16
  fildyn = 'si.dyn'
  ldisp  = .true.
  nq1 = 4, nq2 = 4, nq3 = 4
/
```

`ldisp=.true.` 时不在 namelist 后写 q 点（INPUT_PH：单独的 `xq` 行只在 `ldisp` 和 `qplot` 都不为真时读取）。ph.x 会自动找出 4×4×4 网格中的不可约 q 点，逐个计算。INPUT_PH 说明，`ldisp=.true.` 时输入的 `epsil` 会被忽略：程序对绝缘体自动在 Γ 点计算介电常数和 Born 电荷，对金属则不计算。

```bash
mpirun -np 8 ph.x -nk 4 -in si.ph.in > si.ph.out
```

### 5.2 输出与文件

```bash
grep -A12 'uniform grid of q-points' si.ph.out
ls si.dyn*
```

```text
     Dynamical matrices for ( 4, 4, 4)  uniform grid of q-points
     (   8 q-points):
       N         xq(1)         xq(2)         xq(3)
       1   0.000000000   0.000000000   0.000000000
       2  -0.250000000   0.250000000  -0.250000000
       ...
```

（示意输出，数值仅供参考；不可约 q 点的个数由晶体对称性决定。）

| 文件 | 内容 |
|---|---|
| `si.dyn0` | q 网格大小和不可约 q 点列表（q2r.x 依靠它找到其余文件） |
| `si.dyn1` … `si.dynN` | 每个不可约 q 点（及其星）的动力学矩阵；Γ 点文件中还包含介电张量和 Born 有效电荷（绝缘体） |
| `outdir/_ph0/prefix.phsave/` | 每个 q、每个 irrep 的中间结果（`dynmat.#iq.#irr.xml`），用于 `recover` 和任务拆分 |

---

## 6. q2r.x 与 matdyn.x

### 6.1 q2r.x：实空间力常数

```fortran
&INPUT
  fildyn = 'si.dyn'
  zasr   = 'crystal'
  flfrc  = 'si444.fc'
/
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `fildyn` | 与 ph.x 相同 | q2r.x 读取 `fildyn0` 获取网格信息，再读 `fildyn1…N` |
| `flfrc` | 如 `'si444.fc'` | 输出的实空间力常数文件（必需） |
| `zasr` | `'crystal'` | 对 **Born 有效电荷**施加的 ASR；非极性材料无影响 |
| `loto_2d` | `.false.` | 二维材料的 LO–TO 处理 |
| `write_lr` | `.false.` | 写出长程力常数；只有 matdyn.x 用 `asr='all'` 处理红外活性固体时才需要 |

非共线/自旋轨道计算中 ph.x 输出 xml 格式的动力学矩阵，这时 `fildyn` 要加 `.xml` 扩展名。

```bash
q2r.x -in si.q2r.in > si.q2r.out
```

### 6.2 matdyn.x：色散

```fortran
&INPUT
  asr            = 'crystal'
  flfrc          = 'si444.fc'
  flfrq          = 'si.freq'
  q_in_band_form = .true.
/
6
  0.000  0.000  0.000  40
  1.000  0.000  0.000  20
  1.000  0.500  0.000  20
  1.000  1.000  0.000  60
  0.000  0.000  0.000  40
  0.500  0.500  0.500   1
```

q 点列表的格式：第一行是点数 `nq`，之后每行 `q_x q_y q_z nptq`。`q_in_band_form=.true.` 时，`nptq` 是**该点到下一点之间**生成的点数。坐标默认是**笛卡尔坐标、2π/a 单位**；设 `q_in_cryst_coord=.true.` 可改用晶体坐标。

对 fcc 晶格（2π/a 笛卡尔单位）：Γ = (0,0,0)、X = (1,0,0)、W = (1,1/2,0)、L = (1/2,1/2,1/2)、K = (3/4,3/4,0)。本例路径 Γ–X–W–X'–Γ–L 中 X' = (1,1,0) 与 X 等价，X'→Γ 这一段经过 K 点，所以常被标作 Γ–X–W–K–Γ–L。

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `flfrc` | q2r.x 的输出 | 必需 |
| `asr` | `'crystal'` | `'no'`（默认）、`'simple'`、`'crystal'`、`'all'`、`'one-dim'`、`'zero-dim'` |
| `flfrq` | 如 `'si.freq'` | 频率文件，默认 `'matdyn.freq'`；同时写出 `flfrq.gp`（便于作图） |
| `flvec` / `fleig` | 默认 `'matdyn.modes'` / `'matdyn.eig'` | 归一化位移 / 动力学矩阵本征矢 |
| `q_in_band_form` | `.true.` | 按路径生成 q 点 |
| `q_in_cryst_coord` | `.false.` | q 点是否用晶体坐标 |
| `dos` | 色散时 `.false.` | `.true.` 时计算 DOS，忽略 q 点列表 |
| `nk1, nk2, nk3` | 如 20, 20, 20 | DOS 的均匀 q 网格（`dos=.true.` 时必需） |
| `deltaE` | 1.0 cm⁻¹ | DOS 能量步长 |
| `degauss` | 0（四面体） | DOS 展宽（cm⁻¹）；0 表示用四面体方法 |
| `fldos` | 默认 `'matdyn.dos'` | DOS 输出文件 |
| `amass` | 默认从 `flfrc` 读取 | 可替换原子质量做同位素研究 |
| `loto_disable` | `.false.` | 不在 Γ 点施加 LO–TO 劈裂 |
| `read_lr` | `.false.` | 与 q2r.x 的 `write_lr` 配合 |

`asr='all'` 同时施加 3 个平移、3 个转动 ASR 和 15 个 Huang 条件（应力为零）。INPUT_MATDYN 提醒：对红外活性固体使用 `asr='all'` 时，q2r.x 中要设 `write_lr=.true.`，matdyn.x 中要设 `read_lr=.true.`。

```bash
matdyn.x -in si.matdyn.in > si.matdyn.out
head -3 si.freq.gp
```

`si.freq.gp` 第一列是沿路径累计的 q 距离（2π/a 单位），其后每列是一支声子频率（cm⁻¹）；负值表示虚频。源码 `matdyn.f90` 中，ω² < 0 时频率取 −√|ω²|。

### 6.3 matdyn.x：声子 DOS

```fortran
&INPUT
  asr   = 'crystal'
  flfrc = 'si444.fc'
  flfrq = 'si.dos.freq'
  dos   = .true.
  fldos = 'si.phdos'
  nk1 = 20, nk2 = 20, nk3 = 20
  deltaE = 1.0
/
```

`fldos` 第一行是注释 `# Frequency[cm^-1] DOS PDOS`，第一列为频率（cm⁻¹），第二列为总 DOS（states/cm⁻¹），后面是各原子的投影 DOS。INPUT_MATDYN 说明 DOS **归一化到 3·nat**（声子支数），并提醒四面体方法在非立方晶系中可能工作不正常，这时可以设 `degauss`（cm⁻¹）改用展宽。

### 6.4 绘图

本仓库脚本 `qe_phonon_plot.py` 读取 `flfrq.gp` 和 `fldos`，画色散和 DOS 并标出虚频：

```bash
python ../../scripts/qe_phonon_plot.py -h        # 查看用法（教程目录下为 python ../scripts/qe_phonon_plot.py -h）
```

也可以直接用 gnuplot：

```bash
gnuplot -e "set terminal png; set output 'si_phonon.png'; set ylabel 'Frequency (cm^{-1})'; \
            plot for [i=2:7] 'si.freq.gp' u 1:i w l lc rgb 'black' notitle"
```

QE 自带的 plotband.x 也能处理 `flfrq` 文件（交互式输入），这里不展开。

### 6.5 成功标准

1. ph.x 输出中每个 irrep 都 `Convergence has been achieved`，没有 `No convergence`；
2. 施加 ASR 后，Γ 点三支声学模为 0，靠近 Γ 时线性下降；
3. 整个路径上没有虚频（`si.freq.gp` 中没有明显负值）；
4. 网格上的 q 点处，插值频率与 ph.x 直接计算的频率一致（插值的自洽性检查）；
5. 把 q 网格从 4×4×4 加密到 6×6×6 时，色散变化小于你设定的容差；
6. DOS 积分等于 3·nat。

---

## 7. 金属：fcc Al

### 7.1 要点

金属的 DFPT 有几点不同：

- SCF 必须用 `occupations='smearing'`。PHonon 用户手册 Troubleshooting：如果 ph.x 提示 `occupation numbers probably wrong`，就是金属或自旋极化体系没有设置 smearing；
- 不要设 `epsil=.true.`（报错 `no elec. field with metals`）；`ldisp` 时程序会自动对金属关闭 `epsil`；
- 声子对 k 点和 `degauss` 收敛**很慢**，用户手册还特别提到半芯态体系、以及 q 与 k 网格不公度时收敛更慢；
- 有 Kohn 异常的金属需要更密的 q 网格。

### 7.2 输入

`al.scf.in` 的关键部分：

```fortran
&SYSTEM
  ibrav       = 2
  celldm(1)   = 7.65
  nat         = 1
  ntyp        = 1
  ecutwfc     = 30.0
  ecutrho     = 240.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr = 1.0d-12
/
K_POINTS automatic
  16 16 16 0 0 0
```

`al.ph.in` 与 Si 的形式相同（`fildyn='al.dyn'`，4×4×4 q 网格）；`al.q2r.in`、`al.matdyn.in` 只改文件名。Al 每个原胞只有一个原子，只有 3 支声学支，`asr='simple'` 已经足够。

```bash
mpirun -np 8 pw.x -nk 4 -in al.scf.in > al.scf.out
mpirun -np 8 ph.x -nk 4 -in al.ph.in > al.ph.out
q2r.x -in al.q2r.in > al.q2r.out
matdyn.x -in al.matdyn.in > al.matdyn.out
```

### 7.3 k 点与 smearing 的联合收敛

建议用一个高对称 q 点（例如 X 点）上的频率做联合测试：

| k 网格 | `degauss`（Ry） | ω(X)（cm⁻¹） |
|---|---|---|
| 12×12×12 | 0.04 | 你的结果 |
| 16×16×16 | 0.02 | 你的结果 |
| 20×20×20 | 0.02 | 你的结果 |
| 24×24×24 | 0.01 | 你的结果 |

（表格只是测试设计，没有给出数值。）当 k 网格加密、`degauss` 减小后，频率变化小于 1–2% 时可以认为收敛。单 q 点的测试很便宜：ph.x 输入中把 `ldisp` 去掉，在 namelist 后写 X 点坐标 `1.0 0.0 0.0` 即可。

---

## 8. 虚频诊断

色散中出现虚频（负频率）时，先区分是**数值问题**还是**真实的结构不稳定**。

| 现象 | 可能原因 | 检查 / 解决 |
|---|---|---|
| 只有 Γ 点附近的声学支出现小负值 | ASR 未施加或施加不充分 | matdyn.x 设 `asr='crystal'`；PHonon 手册建议用 dynmat.x 施加 ASR，若声学模降到 1 cm⁻¹ 以下而其他模基本不变，结果可信 |
| 声学支在 Γ 附近的小 q 处向下弯曲 | 力常数截断（q 网格太粗）、`conv_thr`/`tr2_ph` 不够 | 收紧 `conv_thr`（≤1e-12）和 `tr2_ph`（1e-16）；加密 q 网格 |
| 网格 q 点上的直接计算没有虚频，插值后出现 | 插值误差，q 网格太粗 | 加密 q 网格；检查非解析项（极性材料，第13章） |
| 网格 q 点上 ph.x 直接得到大的负频率 | 结构没有充分弛豫，或确实不稳定 | 收紧 `forc_conv_thr` 重新弛豫；检查截断能、k 点（金属尤其） |
| 对称性错误、频率简并被破坏 | 原子位置接近但不精确在对称位置 | 用正确的 `ibrav` 和 Wyckoff 位置重做 SCF（PHonon 手册建议） |
| 低维材料（二维）的 ZA 支在 Γ 附近出现小虚频 | 转动 ASR、真空中 GGA 的数值噪声 | `asr='all'` 或 `'one-dim'`/`'zero-dim'` 的变体；`loto_2d=.true.`；提高 `ecutrho` |
| 确认是真实的软模 | 高温相在 0 K 不稳定（如钙钛矿的八面体旋转） | 沿软模本征矢扰动结构、在超胞中重新弛豫；有限温度稳定性需要 SSCHA 等非谐方法 |

PHonon 用户手册对"糟糕的声子"给出的检查清单：读错数据文件、原子质量设错、`conv_thr` 或 `tr2_ph` 太大、`ecutwfc`（以及 US/PAW 的 `ecutrho`）不够、金属的 k 点不够。

---

## 9. 续算、拆分与并行

### 9.1 `recover`

ph.x 把每个 q、每个 irrep 的结果写入 `outdir/_ph0/prefix.phsave/`。作业中断后，在 `&INPUTPH` 中加入 `recover = .true.`，用同一个输入重新运行即可：已完成的 irrep 不会重算。配合 `max_seconds`（略小于作业时限）可以让 ph.x 平稳停止，避免写坏重启文件。

如果报错无法 recover，用户手册建议删除 `outdir` 中的 `recover*` 文件；`prefix.phsave` 中的结果仍然可以复用。

### 9.2 `start_q`/`last_q` 与 `start_irr`/`last_irr`

大体系可以把 q 点（`start_q`、`last_q`，仅 `ldisp` 时有效）和 irrep（`start_irr`、`last_irr`）分给不同的作业：

```text
Si phonons, q points 1-4 only
&INPUTPH
  prefix  = 'si'
  outdir  = './tmp'
  tr2_ph  = 1.0d-16
  fildyn  = 'si.dyn'
  ldisp   = .true.
  nq1 = 4, nq2 = 4, nq3 = 4
  start_q = 1
  last_q  = 4
/
```

INPUT_PH 的附加说明描述了完整的 GRID 流程：先用 `start_irr=0, last_irr=0` 预先生成位移模式文件（`patterns.#iq.xml`），拷贝到各个任务的 `prefix.phsave` 中；各任务计算互补的 q/irrep 子集；最后把所有 `dynmat.#iq.#irr.xml`（包括 `dynmat.#iq.0.xml`）收集到同一个 `prefix.phsave`，不带 `start_q`/`last_q`/`start_irr`/`last_irr` 再运行一次 ph.x，收集并对角化。注意**不同任务生成的位移模式一般不同，不能混用**，所以必须共享同一套 patterns 文件。官方示例在 PHonon/examples 的 `GRID_example` 目录。

### 9.3 image 并行（`-ni`）

PHonon 用户手册（Parallelism）：ph.x 支持 image、平面波和 k 点（pool）并行。image 并行把不同的 q 点或 irrep 分配给不同的进程组：

```bash
# pw.x: no -ni; the same -nk
mpirun -np 8  pw.x -nk 2 -in si.scf.in > si.scf.out
# ph.x: 8 images x 8 ranks, 2 pools inside each image
mpirun -np 64 ph.x -ni 8 -nk 2 -in si.ph.in > si.ph.out
# collect: run again without -ni and with recover = .true.
mpirun -np 8  ph.x -nk 2 -in si.ph_collect.in > si.ph_collect.out
```

手册说明：image 并行结束后，各 image 计算的 q 点的动力学矩阵不在工作目录中，需要**不带 `-ni`**、加 `recover=.true.` 再运行一次 ph.x 来收集。示例目录的 `job_ph_slurm.sh` 在 `NI>1` 时自动完成这一步。image 之间负载不均衡时，所有 image 都要等最慢的那个完成。

通用的作业脚本见 `job_qe_slurm.sh`（用 `bash ../scripts/job_qe_slurm.sh -h` 查看用法），作业调度基础见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)。

---

## 10. 与 phonopy 有限位移对比

phonopy 通过 `--qe` 接口读取 pw.x 输入、生成带位移的超胞，再从 pw.x 输出中读取力：

```bash
phonopy-init --qe -d --dim="2 2 2" -c si.unitcell.in        # phonopy >= 4 (phonopy-init)
# ... build full pw.x inputs for supercell-001.in etc. and run pw.x ...
phonopy-init -f si-001.out                                   # creates FORCE_SETS
phonopy -p --config band.conf
```

（命令格式以 phonopy 官方 QE 页面为准；phonopy 4.x 起准备步骤使用 `phonopy-init`，旧版本为 `phonopy --qe -d ...`。phonopy 只读取 `ibrav=0` 形式的 pw.x 结构（`CELL_PARAMETERS` + `ATOMIC_POSITIONS`），所以 `si.unitcell.in` 要从本章的 `ibrav=2` 输入改写成 `ibrav=0`；生成的 `supercell-XXX.in` 只含结构部分，需要自己拼上 `&CONTROL`/`&SYSTEM`/`&ELECTRONS`、`ATOMIC_SPECIES` 和 `K_POINTS`。超胞计算中带位移的结构**不能弛豫**。）

| 方面 | DFPT（ph.x） | 有限位移（phonopy + pw.x） |
|---|---|---|
| 计算单元 | 原胞，任意 q | 超胞，等价于 q 网格 |
| 数值参数 | `tr2_ph` | 位移幅度、超胞大小 |
| 支持的泛函 | LDA/GGA（含 DFT+U 的部分情形）；**不支持**杂化泛函、meta-GGA、TS-vdW | 任何能算力的设置 |
| 极性材料 | `epsil` 直接给出 ε∞ 和 Z* | 需要另外提供 BORN（可由 ph.x 计算） |
| 代价 | 每个 irrep 一次线性响应 | 每个不等价位移一次超胞 SCF |
| 后处理 | q2r.x、matdyn.x | phonopy（热力学、群速度、QHA、phono3py） |

两种方法在相同泛函、赝势和收敛参数下应给出一致的频率（差异在百分之几以内），这是一个很好的交叉检查。DFT+U 声子在 ph.x 中有限制（源码中只支持 `atomic` 投影、单个 Hubbard 通道，不支持与 Raman、电声同时计算），超出这些限制时 phonopy 是更稳妥的选择。

---

## 参数选择与收敛建议

| 参数 | 建议 |
|---|---|
| `conv_thr`（pw.x） | ≤ 1e-12（小体系）；大体系按原子数放宽，但一般不宽于 1e-10 |
| `tr2_ph` | 1e-14 ~ 1e-16；默认 1e-12 对声学支往往不够 |
| `ecutwfc`/`ecutrho` | 用于声子时比总能收敛更严；以某个高对称点的光学模频率变化 < 1–2% 为准 |
| k 网格 | 半导体可与 SCF 相同；金属需要显著更密，并与 `degauss` 联合测试（§7.3） |
| q 网格 | 半导体从 4×4×4 起步；检查 6×6×6 时色散的变化；q 网格最好与 k 网格公度 |
| `asr` | matdyn.x 用 `'crystal'`；二维材料考虑 `'all'` 并配合 `write_lr`/`read_lr` |
| 结构 | 声子计算前重新做严格的 `relax`/`vc-relax`，再单独跑 `scf` |
| 对称性 | 保留正确的 `ibrav` 和 Wyckoff 位置；不要让 `nosym` 出现在声子的基态计算中 |
| 并行 | 先用 `-nk` 做 pool 并行；q 点多时加 `-ni`，记得最后用 `recover` 收集 |
| 原子质量 | ph.x 默认从数据文件读取；同位素效应可在 ph.x 或 matdyn.x 中用 `amass` 替换 |

---

## 常见问题与报错

**1. `cannot start from pw.x data file using Gamma-point tricks`**
原因：SCF 用了 `K_POINTS gamma`。解决：改为 `K_POINTS automatic 1 1 1 0 0 0`（分子）或正常的 k 网格，重新 SCF。

**2. `The phonon code is not working after vc-relax`**
原因：ph.x 直接读取了 `vc-relax` 的数据。解决：用优化后的晶格常数和坐标重新运行 `calculation='scf'`。

**3. `no elec. field with metals`**
原因：金属（smearing）体系中设了 `epsil=.true.`。解决：去掉 `epsil`；金属没有宏观介电常数。

**4. `gamma is needed for elec.field`**
原因：在 q ≠ 0 的单 q 计算中设了 `epsil`、`zue`、`lraman` 或 `elop`。解决：这些量只能在 Γ 点计算，单独做一个 q = 0 的计算。

**5. ph.x 提示 `occupation numbers probably wrong` 后继续运行。**
原因：金属或自旋极化体系的 SCF 没有设置 smearing（PHonon 用户手册）。解决：SCF 中设 `occupations='smearing'` 和合适的 `degauss`，重新计算。

**6. 输出中出现 `No convergence has been achieved`。**
原因：线性响应自洽在 `niter_ph`（默认 100）步内未达到 `tr2_ph`。解决：减小 `alpha_mix(1)`（如 0.3）、增大 `nmix_ph`（8–20）；检查基态 `conv_thr` 是否足够；金属检查 smearing 设置。

**7. `The phonon code with hybrid functionals is not yet available` / `... meta-GGA functionals ...` / `... TS-VdW is not yet available`**
原因：ph.x 不支持这些泛函或 vdW 方法。解决：用 phonopy 有限位移（§10），或换成 ph.x 支持的泛函（LDA/GGA；DFT-D3 等其他 vdW 处理请以你所用版本的源码检查为准）。

**8. `The phonon code with US-PP and raman or elop not yet available` / `The phonon code with paw and raman or elop is not yet available` / `third order derivatives not implemented with GGA`**
原因：`lraman`/`elop` 只支持 NC 赝势（`phq_readin.f90`），且不支持 GGA 类泛函（`phq_setup.f90`），即需要 NC + LDA。解决：第13章中改用 LDA 的 NC 赝势单独计算 Raman 张量。

**9. `cannot recover` 或读取 recover 文件出错。**
原因：之前的作业异常终止，重启文件损坏。解决：删除 `outdir` 中的 `recover*` 文件（PHonon 用户手册），保留 `_ph0/prefix.phsave`，再以 `recover=.true.` 运行。

**10. `Wrong degeneracy`、`Wrong representation`、`symmetry operation is non orthogonal` 等对称性错误。**
原因：原子位置接近但不精确等于对称位置，或 q 点与高对称点相差约 1e-5（`star_q` 的容差）。解决：用正确的 `ibrav` 和 Wyckoff 位置重做 SCF；把 q 点写成精确的高对称坐标。

**11. q2r.x 报找不到 `fildyn0` 或读入的 q 点数不对。**
原因：`fildyn` 与 ph.x 不同，或 ph.x 用了 image 并行但没有执行收集步骤，导致部分 `dynN` 文件缺失。解决：统一文件名；`-ni` 运行后再不带 `-ni`、加 `recover=.true.` 运行一次 ph.x。

**12. ph.x 用 `lshift_q=.true.` 后 q2r.x 无法使用。**
原因：INPUT_PH 明确说明半步偏移的 q 网格不能用 q2r.x。解决：做色散时不要设 `lshift_q`。

**13. matdyn.x 色散的 q 路径看起来错位，高对称点不对。**
原因：把晶体坐标的 q 点当作笛卡尔坐标（2π/a）输入，或反之。解决：默认是笛卡尔 2π/a；用晶体坐标时设 `q_in_cryst_coord=.true.`。

---

## 练习

1. 只改 `tr2_ph`（1e-12、1e-14、1e-16），比较 Si 在 Γ 点的声学模频率（不施加 ASR），理解阈值的影响。
2. 把 Si 的 q 网格从 2×2×2、4×4×4 改到 6×6×6，比较 L 点和 X 点的 TA 支频率，判断 4×4×4 是否足够。
3. 在 matdyn.x 中分别用 `asr='no'`、`'simple'`、`'crystal'`，比较 Γ 附近声学支的形状。
4. 用 `amass` 把 Si 的质量改成 ²⁸Si 和 ³⁰Si，比较光学模频率比是否接近 √(30/28)。
5. 对 Al 按 §7.3 的表做 k 点与 `degauss` 联合测试，画出 ω(X) 的收敛曲线。
6. 用 `start_q`/`last_q` 把 Si 的 4×4×4 计算拆成两个作业，最后合并并与一次性计算的结果比较。
7. 用 phonopy（2×2×2 常规超胞，或与 4×4×4 q 网格等价的超胞）计算 Si 的色散，与 DFPT 结果叠加作图。

---

## 参考资料

**官方文档**

- ph.x 输入说明 INPUT_PH：<https://www.quantum-espresso.org/Doc/INPUT_PH.html>
- q2r.x 输入说明 INPUT_Q2R：<https://www.quantum-espresso.org/Doc/INPUT_Q2R.html>
- matdyn.x 输入说明 INPUT_MATDYN：<https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html>
- dynmat.x 输入说明 INPUT_DYNMAT：<https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html>
- PHonon 用户手册（使用方法、并行、Troubleshooting）：<https://www.quantum-espresso.org/Doc/ph_user_guide/>
- 官方 PHonon 示例：<https://gitlab.com/QEF/q-e/-/tree/develop/PHonon/examples>（`example01` 起；`GRID_example`、`Image_example`、`Recover_example` 分别演示拆分、image 并行和续算）
- phonopy 的 QE 接口：<https://phonopy.github.io/phonopy/qe.html>

**课程与教程（只链接，不复制）**

- Cornell PARADIM 暑期学校 2018（Giustino、Pamuk）：<https://www.paradim.org/cu_ss_2018>。上机包括金刚石声子和 GaAs 的 LO–TO，每讲有视频，原理讲解清晰；基于 QE 6.x 和 LDA 赝势，下载/编译部分已过时。
- Cornell PARADIM 2016：<https://www.paradim.org/cu_ss_2016>，包含声子色散一讲。
- Imperial College MSE404-MM：<https://imperialcollegelondon.github.io/MSE404-MM/>，其中的声子实验（DFPT、q2r、matdyn）与本章流程一致，持续维护，强烈推荐配合练习。
- Advanced QE School 2023（EPFL THEOS / Materials Cloud）：<https://github.com/materialscloud-org/QuantumESPRESSO-school-2023>，DFPT 进阶上机。
- pranabdas/espresso：<https://pranabdas.github.io/espresso/>，含声子色散与 Jupyter 绘图。
- AlexanderKvashnin/QE_tutorial（Skoltech）：<https://github.com/AlexanderKvashnin/QE_tutorial>，含声子 Lab。
- 非谐与热输运延伸：[phono3py](https://github.com/phonopy/phono3py)（晶格热导率）、[thermo_pw](https://github.com/dalcorso/thermo_pw)（准谐近似、热膨胀）、[python-sscha](https://github.com/SSCHAcode/python-sscha)（随机自洽谐波近似）。

**本仓库**

- 示例：[12-phonon-dfpt](../examples/12-phonon-dfpt/README.md)
- 脚本：`qe_phonon_plot.py`（色散与 DOS 作图、虚频标记）、`job_qe_slurm.sh`；用 `python ../scripts/<name> -h` 查看用法。
- 对照阅读：[VASP 第12章 声子](../../VASP/tutorials/12-phonon.md)、[CP2K 第14章 振动分析与热力学](../../CP2K/tutorials/14-vibration-thermo.md)

**关键文献**

- S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, *Rev. Mod. Phys.* **73**, 515 (2001)，DFPT 综述：<https://doi.org/10.1103/RevModPhys.73.515>
- X. Gonze, *Phys. Rev. B* **55**, 10337 (1997)：<https://doi.org/10.1103/PhysRevB.55.10337>
- X. Gonze, C. Lee, *Phys. Rev. B* **55**, 10355 (1997)，动力学矩阵、Born 电荷与力常数：<https://doi.org/10.1103/PhysRevB.55.10355>
- C. Lin, S. Poncé, N. Marzari, *npj Comput. Mater.* **8**, 236 (2022)，晶格动力学的不变性与平衡条件（`asr='all'` 的依据）：<https://doi.org/10.1038/s41524-022-00920-6>
- A. Togo, I. Tanaka, *Scr. Mater.* **108**, 1 (2015)，phonopy：<https://doi.org/10.1016/j.scriptamat.2015.07.021>
- A. Togo, L. Chaput, T. Tadano, I. Tanaka, *J. Phys.: Condens. Matter* **35**, 353001 (2023)：<https://doi.org/10.1088/1361-648X/acd831>
- G. Prandini et al., *npj Comput. Mater.* **4**, 72 (2018)，SSSP：<https://doi.org/10.1038/s41524-018-0127-2>
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017)：<https://doi.org/10.1088/1361-648X/aa8f79>

---

[← 上一章](11-neb.md) | [返回目录](../README.md) | [下一章 →](13-dielectric-ir-raman.md)
