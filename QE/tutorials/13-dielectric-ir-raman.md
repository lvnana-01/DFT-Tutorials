# 第13章 介电常数、Born 有效电荷、红外与拉曼（Dielectric Tensor, Born Charges, IR & Raman）

本章讲 Quantum ESPRESSO 中与"电场响应"有关的性质：用 `ph.x` 的密度泛函微扰理论（DFPT）计算高频介电张量 ε∞（`epsil`）、Born 有效电荷 Z*（`zeu`/`zue`）和非共振拉曼张量（`lraman`）；用 `dynmat.x` 施加声学求和规则（ASR）、加入 LO-TO 非解析项，得到红外（IR）强度、拉曼活性和静态介电常数 ε0；用 `q2r.x`/`matdyn.x` 画出带 LO-TO 劈裂的声子色散；最后用 `epsilon.x` 计算独立粒子近似（IPA）下的光学介电函数 ε(ω)，并简单介绍 turboTDDFT 与 Yambo（GW/BSE）。示例体系为闪锌矿 GaAs 和 H2O 分子。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 GaAs 的 scf + `ph.x`（`epsil = .true.`、`trans = .true.`、`lraman = .true.`）输入，从输出中读出 ε∞ 与 Z*；
- 用 `dynmat.x` 得到 Γ 点 TO/LO 频率、IR 强度、拉曼活性，并理解 `asr`、`q`、`lperm` 的作用；
- 让 `q2r.x` + `matdyn.x` 自动把 Z* 和 ε∞ 用于非解析修正，画出极性材料的声子色散；
- 对孤立分子（H2O）在大盒子中做 Γ 点振动分析，得到红外和拉曼谱线；
- 知道拉曼计算的限制（**模守恒赝势 + LDA + 非自旋极化**），知道哪些情况下 ph.x 会直接报错；
- 用 `epsilon.x` 做 IPA 光学谱，并知道什么时候需要 TDDFT 或 GW/BSE。

| 目标 | 推荐程序 / 关键词 | 说明 |
|---|---|---|
| 电子部分介电张量 ε∞ | `ph.x`：`epsil = .true.` | 只对绝缘体、只在 q = 0 |
| Born 有效电荷 Z* | `ph.x`：`epsil` + `trans`（默认 `zeu = epsil`） | 极性材料 LO-TO 劈裂、IR 强度都需要它 |
| Γ 点频率、IR 强度、LO 模 | `dynmat.x` | 读 `fildyn`，施加 `asr`，用 `q` 指定 LO 方向 |
| 拉曼张量与拉曼活性 | `ph.x`：`lraman = .true.` + `dynmat.x` | NC + LDA，不能自旋极化、不能 US/PAW |
| 极性材料声子色散 | `ph.x`（`ldisp`）+ `q2r.x` + `matdyn.x` | 绝缘体 `ldisp` 时 epsil 自动打开 |
| 静态介电常数 ε0 | `dynmat.x`：`lperm = .true.` | ε0 = ε∞ + 各极性模贡献 |
| 光学介电函数 ε(ω)（IPA） | `pw.x` nscf + `epsilon.x` | 无局域场、无激子，代价低 |
| 分子吸收谱 / 局域场效应 | turboTDDFT（`turbo_lanczos.x` 等） | 线性响应 TDDFT |
| 准粒子带隙、激子 | Yambo（`p2y` 转换 QE 数据） | GW/BSE，超出本章 |

前置知识：[第12章 声子：ph.x、q2r.x 与 matdyn.x](12-phonon-dfpt.md) 中的 DFPT 流程、`tr2_ph`、`asr`、q 网格；[第6章 能带结构、态密度与投影态密度](06-band-dos-pdos.md) 中的 nscf 与 `nbnd`。

---

## 2. 原理速览

### 2.1 三种二阶导数

把总能 E 对原子位移 u 和宏观电场 **E** 做展开：

$$
\varepsilon^{\infty}_{\alpha\beta} = \delta_{\alpha\beta} - \frac{4\pi}{\Omega}\frac{\partial^2 E}{\partial \mathcal{E}_\alpha \partial \mathcal{E}_\beta},\qquad
Z^{*}_{\kappa,\alpha\beta} = -\frac{\partial^2 E}{\partial \mathcal{E}_\alpha\,\partial u_{\kappa\beta}}
= \frac{\Omega}{e}\frac{\partial P_\alpha}{\partial u_{\kappa\beta}}
= \frac{1}{e}\frac{\partial F_{\kappa\beta}}{\partial \mathcal{E}_\alpha}
$$

- ε∞（"clamped-ion"，离子不动）：只来自电子的极化。
- Z*：一个原子位移时产生多少宏观极化，或者一个电场作用时原子受多大的力。`ph.x` 默认用"d Force / dE"的形式计算（`zeu`），也可以用"d P / du"形式（`zue`）交叉验证，两者在数值误差内应一致。
- 力常数（dynamical matrix）：来自原子位移的二阶导数，见第 12 章。

DFPT 中电场扰动需要波函数对 k 的导数，因此 `ph.x` 在 epsil 计算中会先通过线性方程求解 $\partial\psi/\partial k$ 相关的量（输出中可以看到先做电场扰动、再做原子位移）。`dek`、`eth_rps`、`eth_ns` 是**拉曼（二阶响应）专用**的数值参数（INPUT_PH：分别为对 k 求导的步长、Pc R|ψ⟩ 与非自洽波函数的收敛阈值），一般保持默认值。

### 2.2 LO-TO 劈裂与非解析项

在极性绝缘体中，长波纵光学（LO）声子伴随宏观电场，Γ 点动力学矩阵多出一个与 q 方向有关的非解析项：

$$
\widetilde{C}^{\mathrm{NA}}_{\kappa\alpha,\kappa'\beta}(\mathbf q\to 0) =
\frac{4\pi}{\Omega}\frac{(\mathbf q\cdot Z^{*}_{\kappa})_\alpha\,(\mathbf q\cdot Z^{*}_{\kappa'})_\beta}{\mathbf q\cdot\varepsilon^{\infty}\cdot\mathbf q}
$$

`ph.x` 在 q = 0 处输出的动力学矩阵**不包含**这一项，也**不施加** ASR（PHonon 用户手册明确说明）。所以 Γ 点的正确 TO/LO 频率要靠 `dynmat.x`（单个 Γ 点）或 `matdyn.x`（整条色散）补上。

Lyddane–Sachs–Teller（LST）关系可以用来检查结果是否自洽：

$$
\frac{\varepsilon^{0}}{\varepsilon^{\infty}} = \prod_m \frac{\omega_{\mathrm{LO},m}^2}{\omega_{\mathrm{TO},m}^2}
$$

### 2.3 红外强度与拉曼活性

- **IR 强度**：模 m 的振子强度正比于 $\left|\sum_{\kappa\beta} Z^{*}_{\kappa,\alpha\beta}\,e_{m}(\kappa\beta)/\sqrt{M_\kappa}\right|^2$，因此只要有 Z*，`dynmat.x` 就能给出每个模的 IR 强度，单位为 (D/Å)²/amu。
- **拉曼活性**：非共振拉曼强度由极化率对简正坐标的导数 ∂χ/∂u 决定，它是能量对两个电场和一个位移的三阶导数。QE 用 Lazzeri–Mauri 的二阶响应方法（2n+1 定理）计算拉曼张量，`dynmat.x` 再把它投影到简正模上，给出拉曼活性（Å⁴/amu）和退偏比（depolarization factor）。

> [!IMPORTANT]
> 拉曼张量只对**模守恒（NC）赝势 + LDA** 实现。QE 7.6 源码 `PHonon/PH/phonon.f90` 头部的功能表写明 "Raman tensor NC [1]"，其中 [1] 表示 LDA；交换关联三阶导数例程 `d2mxc` 只实现了 Perdew–Zunger LDA。US/PAW、非共线、自旋极化（LSDA）时 `ph.x` 直接报错退出；GGA 赝势不会报错，但拉曼张量中的三阶交换关联项仍按 LDA 计算，结果只是近似，不应作为正式结果。

### 2.4 ph.x 的功能矩阵（QE 7.6）

下表根据 `PHonon/PH/phonon.f90` 头部注释与 `phq_readin.f90` 中的检查整理：

| 物理量 | NC | US | PAW | 泛函要求 |
|---|---|---|---|---|
| 介电常数 ε∞ | 可以 | 可以 | 可以 | LDA/GGA/LSDA；非共线非磁（含 SOC）可以；非共线**磁性**体系报错 |
| Born 有效电荷 Z* | 可以 | 可以 | 可以 | 同上 |
| 拉曼张量（`lraman`） | 可以 | 报错 | 报错 | 仅 LDA，非自旋极化 |
| 电光张量（`elop`） | 可以 | 报错 | 报错 | 仅 LDA |
| 任何 ph.x 计算 | — | — | — | 杂化泛函、meta-GGA、TS-vdW、外加电场、约束磁化、非周期边界条件都不支持；DFT-D2 与非局域 vdW-DF 支持；DFT-D3 需要另行提供 D3 Hessian（INPUT_PH 的 `dftd3_hess`，实验性） |

> [!WARNING]
> QE 7.6 的 release notes 写明："Calculation of Raman coefficients with second-order response was broken since v.7.4.1"，已在 7.6 修复。如果你用 7.4.1–7.5 算过拉曼，请用 7.6 重算。

---

## 3. 准备工作

示例目录：[../examples/13-dielectric-ir-raman/README.md](../examples/13-dielectric-ir-raman/README.md)。

### 3.1 赝势

因为要算拉曼，本章统一使用 PseudoDojo v0.4 的 **LDA（PW）模守恒标准精度**赝势包：

```bash
mkdir -p pseudo && cd pseudo
wget https://www.pseudo-dojo.org/pseudos/nc-sr-04_pw_standard_upf.tgz
tar -xzf nc-sr-04_pw_standard_upf.tgz Ga.upf As.upf O.upf H.upf
cd ..
grep -m1 'functional=' pseudo/*.upf      # every file: "SLA  PW   NOGX NOGC"  -> LDA
grep -m1 'pseudo_type=' pseudo/*.upf     # every file: "NC"
grep -m1 'z_valence=' pseudo/Ga.upf      # Ga includes 3d semicore: 13 valence electrons
```

| 元素 | 文件名 | 价电子数 | 说明 |
|---|---|---|---|
| Ga | `Ga.upf` | 13（含 3d） | PseudoDojo 截断提示 normal ≈ 40 Ha（80 Ry） |
| As | `As.upf` | 15（含 3d） | |
| O | `O.upf` | 6 | |
| H | `H.upf` | 1 | |

PseudoDojo 的截断提示单位是 Hartree，QE 的 `ecutwfc` 单位是 Ry（乘 2）。模守恒赝势 `ecutrho = 4 × ecutwfc` 即可。赝势选择的一般原则见 [第3章 赝势选择](03-pseudopotentials.md)。

> [!NOTE]
> 这套 LDA 包中的 Ga、As、O 含非线性芯修正（`core_correction="T"`，H 不含）。QE 7.6 的 `PHonon/PH/dvpsi_e2.f90` 在计算拉曼所需的三阶交换关联项时只用价电荷密度（加上芯电荷的那一行被注释掉），因此对含 NLCC 的赝势，拉曼张量中这一项是近似的。需要严格对比时，可以改用无 NLCC 的 LDA 模守恒赝势（PseudoDojo GitHub 上的 `ONCVPSP-LDA-SR-noCCv0.4` 仓库）并比较两者的结果。ε∞、Z* 和频率不受这一点影响。

> [!TIP]
> 如果只需要 ε∞、Z* 和 LO-TO（不算拉曼），可以用 PBE 赝势，包括 SSSP 中的 US/PAW 赝势（例如 SSSP 1.3 Efficiency 的 `Ga.pbe-dn-kjpaw_psl.1.0.0.UPF`、`As.pbe-n-rrkjus_psl.0.2.UPF`）。但同一个计算中所有赝势必须使用同一种泛函。

### 3.2 结构与收敛要求

DFPT 对基态的精度要求比普通 scf 高：

| 项目 | 推荐值 | 原因 |
|---|---|---|
| 结构 | 在**同一泛函、同一赝势、同一截断**下充分弛豫 | 残余应力/力会导致 Γ 点声学模不为零、出现虚频 |
| `conv_thr` | 1e-12 或更小（小体系 1e-14） | DFPT 线性方程的精度受限于基态 |
| `tr2_ph` | 1e-14 ~ 1e-16 | 默认 1e-12，Z* 的 ASR 偏差会偏大 |
| k 网格 | 比总能收敛值更密 | ε∞ 对 k 点比总能敏感，需单独测试 |
| `occupations` | 绝缘体用默认 `'fixed'` | 有 smearing 时 epsil 会报错 |

本章 GaAs 示例用晶格常数 a = 5.65 Å（实验值附近）只作演示；正式计算请先按 [第5章](05-relax-eos.md) 用 LDA 优化晶格常数，否则 ε∞ 与频率都会带上应变误差。

---

## 4. GaAs：Γ 点 ε∞、Z* 与拉曼张量

### 4.1 第 1 步：scf

```fortran
! GaAs (zincblende), LDA, norm-conserving PseudoDojo pseudopotentials.
! Ground state for DFPT: dielectric tensor, Born charges, Raman tensor.
! NC + LDA is REQUIRED if lraman=.true. in ph.x.
&CONTROL
  calculation = 'scf'
  prefix      = 'gaas'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  verbosity   = 'high'
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 2
  ecutwfc   = 80.0      ! start value from PseudoDojo hints (Ga: 40 Ha = 80 Ry); converge it
  ecutrho   = 320.0     ! 4 x ecutwfc for norm-conserving pseudopotentials
  ! occupations default = 'fixed' (insulator). Do NOT use smearing: epsil fails for metals.
/
&ELECTRONS
  conv_thr    = 1.0d-12  ! tight SCF threshold before DFPT
  mixing_beta = 0.5
/
ATOMIC_SPECIES
  Ga  69.723  Ga.upf
  As  74.922  As.upf
CELL_PARAMETERS angstrom
  -2.8250   0.0000   2.8250
   0.0000   2.8250   2.8250
  -2.8250   2.8250   0.0000
ATOMIC_POSITIONS crystal
  Ga  0.00  0.00  0.00
  As  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

```bash
mpirun -np 8 pw.x -nk 4 -in gaas.scf.in > gaas.scf.out
grep -e '^!' -e 'highest occupied' -e 'convergence has been achieved' gaas.scf.out
```

检查：输出中应出现 `highest occupied, lowest unoccupied level`（绝缘体、无 smearing 时的格式），而不是 `the Fermi energy is`。如果出现后者，说明设置了 smearing，下一步 epsil 会失败。

### 4.2 第 2 步：ph.x（epsil + trans + lraman）

```fortran
GaAs at Gamma: epsilon_inf, Born effective charges, Raman tensor
&INPUTPH
  prefix  = 'gaas'
  outdir  = './tmp'
  fildyn  = 'gaas.dynG'
  tr2_ph  = 1.0d-16     ! tighter than the 1e-12 default for clean Z* and ASR
  epsil   = .true.      ! macroscopic dielectric tensor (insulators only, q = 0 only)
  trans   = .true.      ! phonons at q = 0; with epsil -> Born effective charges (zeu)
  lraman  = .true.      ! non-resonant Raman tensor (NC + LDA, no spin, CPU build)
/
0.0 0.0 0.0
```

`&INPUTPH` 之前的第一行是任务标题（必须有）；namelist 之后的一行是 q 点（笛卡尔坐标，单位 2π/a），Γ 点写 `0.0 0.0 0.0`。

| 关键词 | 推荐值 | 说明（据 INPUT_PH） |
|---|---|---|
| `prefix` / `outdir` | 与 scf 相同 | 必须一致，否则找不到 `gaas.save` |
| `fildyn` | `'gaas.dynG'` | 动力学矩阵输出文件（默认 `'matdyn'`）；ε∞、Z*、拉曼张量也写在这里 |
| `tr2_ph` | 1e-14 ~ 1e-16 | 自洽阈值，默认 1e-12 |
| `epsil` | `.true.` | q = 0 的非金属计算宏观介电常数；金属或 q ≠ 0 会报错 |
| `trans` | `.true.`（默认） | 计算声子；与 `epsil` 同时为真时计算有效电荷 |
| `zeu` | 默认等于 `epsil` | 用 dF/dE 计算 Z*（默认算法）；`epsil=.true.` + `zeu=.false.` 时只算介电张量 |
| `zue` | `.false.`（默认） | 用 dP/du 计算 Z*，可与 zeu 交叉验证 |
| `lraman` | `.true.` | 用二阶响应（Lazzeri–Mauri）计算非共振拉曼系数 |
| `elop` | `.false.`（默认） | 电光张量，限制与 `lraman` 相同 |
| `lrpa` / `lnoloc` | `.false.` | 设为真时分别给出 RPA（ΔV_xc = 0）和无局域场（ΔV_Hxc = 0）的介电常数，用于分析。只能单独算介电常数：必须同时设 `trans = .false.`（否则报 `only dielectric constant with lrpa or lnoloc`），且会自动关闭 `zeu` 与 `lraman` |
| `asr` | 晶体不用设 | ph.x 内的 ASR 只与 Γ-Γ 技巧（只有 k = 0）配合使用，适用于分子；晶体在 `dynmat.x`/`matdyn.x` 中施加 |
| `amass(i)` | 不设 | 不设时读取 pw.x 中的质量 |
| `max_seconds` / `recover` | 按需 | 长时间作业的断点续算，见第 12 章 |

```bash
mpirun -np 8 ph.x -nk 4 -in gaas.ph_epsil.in > gaas.ph_epsil.out
```

### 4.3 读 ph.x 输出

```bash
grep -A4 "Dielectric constant in cartesian axis" gaas.ph_epsil.out
grep -A8 "Effective charges (d Force / dE) in cartesian axis without acoustic sum rule" gaas.ph_epsil.out
grep -A8 "Effective charges (d Force / dE) in cartesian axis with asr applied" gaas.ph_epsil.out
grep -A1 "Effective charges Sum" gaas.ph_epsil.out
grep "Raman tensor" gaas.dynG
```

输出格式如下（示意输出，数值仅供参考）：

```text
          Dielectric constant in cartesian axis 

          (      1x.xxxxxxxxx       0.000000000       0.000000000 )
          (       0.000000000      1x.xxxxxxxxx       0.000000000 )
          (       0.000000000       0.000000000      1x.xxxxxxxxx )

          Effective charges (d Force / dE) in cartesian axis without acoustic sum rule applied (asr)

           atom      1 Ga     Mean Z*:        2.xxxxx
      Ex  (        2.xxxxx        0.00000        0.00000 )
      Ey  (        0.00000        2.xxxxx        0.00000 )
      Ez  (        0.00000        0.00000        2.xxxxx )
           atom      2 As     Mean Z*:       -2.xxxxx
      ...
          Effective charges Sum: Mean:        0.0xxxx
```

成功判据：

- 立方晶体的 ε∞ 应为对角阵且三个对角元相等，非对角元为 0（或 1e-6 量级的数值噪声）。
- `Effective charges Sum` 应接近 0（电荷中性 / ASR）。偏离较大（例如 > 0.05）说明 k 点、截断或 `tr2_ph` 不够。
- `with asr applied` 部分中 Ga 和 As 的 Z* 大小相等、符号相反。
- 如果用 `lrpa`/`lnoloc`（需另开一次 `trans = .false.` 的 ph.x 计算），标题会分别变成 `RPA dielectric constant in cartesian axis (DV_xc=0)` 和 `Dielectric constant in cartesian axis (DV_Hxc=0)`。

`gaas.dynG` 中依次写有动力学矩阵、`Dielectric Tensor:`、`Effective Charges E-U: Z_{alpha}{s,beta}`，打开了 `lraman` 时还有 `Raman tensor (A^2)`。这些块正是后面 `dynmat.x` 读取的内容。

> [!NOTE]
> 物理背景：LDA/GGA 低估带隙，通常会让 ε∞ 偏大，这是泛函误差而不是输入错误。比较实验时要说明所用泛函和晶格常数。

### 4.4 第 3 步：dynmat.x —— ASR、LO 模、IR 与拉曼

```fortran
&INPUT
  fildyn = 'gaas.dynG'
  asr    = 'crystal'          ! acoustic sum rule for a 3D crystal
  q(1)   = 1.0                ! direction of q -> 0 used for the LO non-analytic term
  q(2)   = 0.0
  q(3)   = 0.0
  filout = 'gaas.dynmat.modes'
  filmol = 'gaas.dynmat.mold'
  filxsf = 'gaas.dynmat.axsf'
/
```

| 关键词 | 推荐值 | 说明（据 INPUT_DYNMAT） |
|---|---|---|
| `fildyn` | ph.x 的 `fildyn` | 输入动力学矩阵文件，默认 `'matdyn'` |
| `asr` | 晶体 `'crystal'`；分子 `'zero-dim'`；一维 `'one-dim'` | 默认 `'no'`；`'crystal'` 用投影法施加 3 个平移 ASR，`'zero-dim'` 再加 3 个转动 ASR；除 `'no'` 以外都会对 Z* 做简单修正 |
| `q(1:3)` | 想看的 LO 方向，例如 `(1,0,0)` | 默认 (0,0,0) 时不加非解析项，输出中没有 LO-TO 劈裂 |
| `lperm` | 单独一次、不设 `q` | 计算各 Γ 模对介电张量的贡献（得到 ε0）；若同时给了非零 `q`，dynmat.x 打印 `BEWARE: phonon contribution to permittivity computed with TO-LO splitting`，结果混入了 LO 模，不应采用 |
| `lplasma` | 按需 | 输出各模的有效等离子体频率，会自动打开 `lperm` |
| `amass(i)` | 不设 | 默认从 `fildyn` 读入；做同位素替换时修改 |
| `filout` / `fileig` | 按需 | 归一化位移 / 正交本征矢 |
| `filmol` / `filxsf` | 按需 | Molden 格式、XCrySDen 的 axsf 格式，用于可视化振动模 |
| `loto_2d` | 二维材料设 `.true.` | 二维 LO-TO 处理 |

```bash
dynmat.x -in gaas.dynmat.in > gaas.dynmat.out
grep -B2 -A8 "# mode" gaas.dynmat.out
```

输出格式（示意输出，数值仅供参考）：

```text
     IR activities are in (D/A)^2/amu units
     Raman activities are in A^4/amu units
     multiply Raman by ... for Clausius-Mossotti correction

# mode   [cm-1]    [THz]      IR          Raman   depol.fact
    1      0.00    0.0000    0.0000         0.0000    0.0000
    2      0.00    0.0000    0.0000         0.0000    0.0000
    3      0.00    0.0000    0.0000         0.0000    0.0000
    4    2xx.xx    x.xxxx    x.xxxx         x.xxxx    0.xxxx
    5    2xx.xx    x.xxxx    x.xxxx         x.xxxx    0.xxxx
    6    2xx.xx    x.xxxx    x.xxxx         x.xxxx    0.xxxx
```

读法：

- 模 1–3：声学模，施加 ASR 后应为 0（或 < 1 cm⁻¹）。没有 ASR 时可能出现几十 cm⁻¹ 甚至虚频，这是数值问题而非物理。
- 模 4–6：闪锌矿只有一组 T2 光学模，加了 `q = (1,0,0)` 方向的非解析项后，劈裂成 2 个 TO + 1 个 LO（LO 更高）。
- 闪锌矿 T2 模同时具有 IR 和拉曼活性；含反演中心的晶体（例如岩盐结构 MgO、金刚石结构 Si）满足互斥规则：IR 活性模没有拉曼活性，反之亦然。可以用这一点快速检查结果。
- ε0 要用**不带 `q`** 的另一次 `dynmat.x` 计算（示例目录中的 `gaas.dynmat_eps0.in`，`lperm = .true.`），它输出各模对介电张量的贡献，求和加上 ε∞ 即 ε0，可以用 LST 关系交叉检查：

```fortran
&INPUT
  fildyn = 'gaas.dynG'
  asr    = 'crystal'
  lperm  = .true.             ! q is left at (0,0,0): TO modes only, as required for eps0
  filout = 'gaas.dynmat_eps0.modes'
/
```

```bash
dynmat.x -in gaas.dynmat_eps0.in > gaas.dynmat_eps0.out
grep -A4 "Electronic dielectric permittivity tensor" gaas.dynmat_eps0.out          # eps_inf
grep -A4 "with zone-center polar mode contributions" gaas.dynmat_eps0.out         # eps0
```

> [!TIP]
> 把 `q` 改成 (0,0,1)、(1,1,0)、(1,1,1) 分别运行一次，立方晶体的 LO 频率应保持不变；非立方晶体（如纤锌矿 ZnO、GaN）的 LO 频率会依赖方向，这正是需要指定 `q` 的原因。

### 4.5 画红外/拉曼谱

`dynmat.x` 给出的是"线谱"（频率 + 强度）。实验谱有线宽，常用洛伦兹展宽：

```python
# broaden_ir_raman.py -- illustrative only; parse the "# mode" table of dynmat.x output
import numpy as np
import matplotlib.pyplot as plt

rows = []
with open("gaas.dynmat.out") as f:
    in_table = False
    for line in f:
        if line.startswith("# mode"):
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) < 4 or not parts[0].isdigit():
                break
            rows.append([float(x) for x in parts[1:]])
data = np.array(rows)              # columns: cm-1, THz, IR, (Raman, depol.)
freq, ir = data[:, 0], data[:, 2]
w = np.linspace(0, 400, 2001)
gamma = 5.0                        # Lorentzian HWHM in cm-1 (choose to match experiment)
spec = sum(a * gamma / np.pi / ((w - f) ** 2 + gamma ** 2) for f, a in zip(freq, ir) if f > 10)
plt.plot(w, spec)
plt.xlabel("Frequency (cm$^{-1}$)")
plt.ylabel("IR intensity (arb. units)")
plt.savefig("gaas_ir.png", dpi=150)
```

拉曼列（第 5 列）同理。实验拉曼强度还与激发光频率、温度（Bose 因子）和散射几何有关，`dynmat.x` 给出的是非共振、各向平均的活性。

---

## 5. LO-TO 劈裂进入声子色散：ph.x + q2r.x + matdyn.x

单个 Γ 点只能看 TO/LO 频率；整条色散需要在 q 网格上做 DFPT，然后傅里叶插值。对极性绝缘体，长程偶极-偶极相互作用会让力常数衰减很慢，QE 的做法是：在 `q2r.x` 中把非解析部分（由 Z* 和 ε∞ 决定）从力常数中扣除，`matdyn.x` 插值时再加回来。只要 Γ 点的动力学矩阵文件里有 Z* 和 ε∞，这一过程是自动的。

### 5.1 ph.x 在 q 网格上

```fortran
GaAs phonon dispersion on a 4x4x4 q grid (epsil/zeu are set automatically for insulators)
&INPUTPH
  prefix  = 'gaas'
  outdir  = './tmp'
  fildyn  = 'gaas.dyn'
  tr2_ph  = 1.0d-14
  ldisp   = .true.
  nq1 = 4, nq2 = 4, nq3 = 4
/
```

INPUT_PH 对 `epsil` 的说明：`ldisp = .true.` 时输入的 `epsil` 会被忽略，代码对绝缘体自动设为 `.true.`、对金属设为 `.false.`。所以这里不写 `epsil`。

```bash
mpirun -np 16 ph.x -nk 4 -in gaas.ph_disp.in > gaas.ph_disp.out
ls gaas.dyn*                                 # gaas.dyn0 (q list) + gaas.dyn1 ... gaas.dynN
grep -c "Dielectric Tensor" gaas.dyn1         # should be 1: Z* and eps are stored at Gamma
```

### 5.2 q2r.x 与 matdyn.x

```fortran
&INPUT
  fildyn = 'gaas.dyn'
  zasr   = 'crystal'     ! ASR on the Born effective charges
  flfrc  = 'gaas.fc'
/
```

```fortran
&INPUT
  asr            = 'crystal'
  flfrc          = 'gaas.fc'
  flfrq          = 'gaas.freq'
  q_in_band_form = .true.
  q_in_cryst_coord = .true.
/
5
  0.000  0.500  0.500  40   ! X  (crystal coordinates of the fcc primitive cell above)
  0.000  0.000  0.000  40   ! Gamma
  0.000  0.500  0.000  40   ! L
  0.000  0.625  0.375  10   ! U
  0.000  0.500  0.500   1   ! X
```

| 程序 | 关键词 | 说明 |
|---|---|---|
| q2r.x | `fildyn`、`flfrc` | 必需；`fildyn` 与 ph.x 相同（不带编号） |
| q2r.x | `zasr` | 对 Born 有效电荷施加 ASR，可选 `'simple'`、`'crystal'`、`'one-dim'`、`'zero-dim'` |
| q2r.x | `loto_2d` | 二维材料的 LO-TO 处理 |
| matdyn.x | `asr` | 对插值后的力常数施加 ASR（默认 `'no'`） |
| matdyn.x | `q_in_band_form` / `q_in_cryst_coord` | 只给每段端点 + 点数；端点用晶体坐标 |
| matdyn.x | `flfrq` | 频率输出文件（还会写 `flfrq.gp` 供 gnuplot 使用） |

> [!NOTE]
> 上面的 X、L、U 坐标是**相对于本例 `CELL_PARAMETERS` 那组原胞矢量**的倒格子晶体坐标（与 `ibrav = 2` 的约定相同）。如果你换了原胞矢量，必须重新换算；用 SeeK-path 等工具生成路径时，要使用工具给出的标准原胞。

```bash
q2r.x    -in gaas.q2r.in    > gaas.q2r.out
matdyn.x -in gaas.matdyn.in > gaas.matdyn.out
python ../../scripts/qe_phonon_plot.py -h    # run from the example dir; check usage, then plot gaas.freq.gp
```

检查：色散在 Γ 点两侧光学支应该"跳"一下——沿 X–Γ 与 Γ–L 趋近 Γ 时都能看到 LO 高于 TO。如果 Γ 点光学支完全简并，说明 Z* 没有被读入（检查 `gaas.dyn1` 是否含 `Dielectric Tensor`、scf 是否误加了 smearing）。色散的其余收敛问题（q 网格、虚频）见 [第12章](12-phonon-dfpt.md)。

---

## 6. 分子：H2O 的 Γ 点振动、IR 与拉曼

孤立分子用"大盒子 + 单个 k 点"处理。有两个容易踩的坑：

1. `ph.x` 不能读取用 `K_POINTS gamma`（Γ 点技巧）生成的数据，会报 `cannot start from pw.x data file using Gamma-point tricks`。要写成 `K_POINTS automatic` + `1 1 1 0 0 0`。
2. `ph.x` 不支持非周期边界条件，所以不要设置 `assume_isolated`。盒子要足够大（10–15 Å），并检查频率随盒子大小的变化。

### 6.1 弛豫

```fortran
! H2O molecule in a 12 A cubic box, LDA, NC pseudopotentials.
! Relax tightly, then run ph.x at Gamma directly on the relaxed structure.
! K_POINTS must NOT be 'gamma': ph.x cannot read data written with Gamma tricks.
&CONTROL
  calculation   = 'relax'
  prefix        = 'h2o'
  outdir        = './tmp'
  pseudo_dir    = './pseudo'
  forc_conv_thr = 1.0d-5
  etot_conv_thr = 1.0d-7
/
&SYSTEM
  ibrav     = 1
  A         = 12.0
  nat       = 3
  ntyp      = 2
  ecutwfc   = 90.0      ! start value (H/O PseudoDojo hints are ~37-42 Ha); converge it
  ecutrho   = 360.0
  ! no assume_isolated: ph.x does not support non-periodic boundary conditions
/
&ELECTRONS
  conv_thr    = 1.0d-12
  mixing_beta = 0.5
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  O  15.999  O.upf
  H   1.008  H.upf
ATOMIC_POSITIONS angstrom
  O   6.0000   6.0000   6.0000
  H   6.7572   6.5865   6.0000
  H   5.2428   6.5865   6.0000
K_POINTS automatic
  1 1 1 0 0 0
```

`relax` 结束时 pw.x 已把最终结构写入 `tmp/h2o.save`，ph.x 可以直接在其上计算（与第 12 章"先 relax 再 scf"的做法等价，但更省事）。若你想保险，也可以把最终坐标复制进一个新的 scf 输入再跑一次。

### 6.2 ph.x 与 dynmat.x

```fortran
H2O normal modes, IR and Raman at Gamma
&INPUTPH
  prefix  = 'h2o'
  outdir  = './tmp'
  fildyn  = 'h2o.dynG'
  tr2_ph  = 1.0d-16
  epsil   = .true.     ! needed for IR intensities (Born charges)
  trans   = .true.
  lraman  = .true.     ! needed for Raman activities (NC + LDA)
  ! no asr here: lraman=.true. switches off the gamma_gamma trick, and ph.x then
  ! resets asr to .false.; the ASR is imposed afterwards by dynmat.x (asr='zero-dim')
/
0.0 0.0 0.0
```

```fortran
&INPUT
  fildyn = 'h2o.dynG'
  asr    = 'zero-dim'     ! 3 translations + 3 rotations for a molecule
  filout = 'h2o.dynmat.modes'
  filmol = 'h2o.dynmat.mold'
  filxsf = 'h2o.dynmat.axsf'
/
```

```bash
mpirun -np 8 pw.x -in h2o.relax.in    > h2o.relax.out
mpirun -np 8 ph.x -in h2o.ph_gamma.in > h2o.ph_gamma.out
dynmat.x -in h2o.dynmat.in > h2o.dynmat.out
grep -A10 "# mode" h2o.dynmat.out
```

当 BZ 求和只有 k = 0 时，ph.x 默认使用 "gamma_gamma" 技巧（INPUT_PH 中 `nogg` 的说明），利用对称性和 ASR 把需要求解的线性响应数目降到最少；ph.x 的 `asr = .true.` 只在这种情况下有效。但 QE 7.6 的 `phq_readin.f90` 中，`lraman = .true.`（或 `modenum > 0`）会关闭 gamma_gamma 技巧，随后 `asr` 被重置为 `.false.`。所以本例的 ph.x 输入不写 `asr`，声学求和规则统一交给 `dynmat.x` 的 `asr = 'zero-dim'`。只算 IR（不开 `lraman`）时，可以在 ph.x 中加 `asr = .true.` 以利用 gamma_gamma 技巧加速，QE 官方 example15 的 CO2 输入就是这样写的。

检查要点：

- 9 个模中前 6 个（3 平动 + 3 转动）应接近 0；转动模受盒子大小、截断和 egg-box 效应影响，常残留几十 cm⁻¹，`asr = 'zero-dim'` 会把它们投影掉。
- 剩下 3 个内振动模依次是弯曲、对称伸缩、反对称伸缩。用 `h2o.dynmat.mold` 在 Molden/Jmol 中可以看到振动动画。
- 分子中 `dynmat.x` 不需要 `q`（没有 LO-TO）。
- 谐振频率与实验基频比较时，要记住非谐效应（O–H 伸缩约差 100–200 cm⁻¹ 量级）和泛函误差。

> [!TIP]
> QE 官方示例 `PHonon/examples/example15` 演示了 CO2 分子和 ZnO 晶体的 pw.x + ph.x + dynmat.x 红外/拉曼流程（其中 `lraman` 被注释掉，因为用的是 PBE 赝势），`example05` 演示了 AlAs 的拉曼张量与电光张量（LDA）。两者都值得对照阅读：[PHonon/examples](https://github.com/QEF/q-e/tree/develop/PHonon/examples)。

---

## 7. 光学介电函数：epsilon.x（独立粒子近似）

`ph.x` 给出的是静态（ω → 0）电子介电常数。要看吸收谱、折射率随频率的变化，最便宜的做法是 PP 包中的 `epsilon.x`：在 RPA 层面计算宏观介电函数，**不含局域场效应**，也不含激子效应。源码 `PP/src/epsilon.f90` 的头注释写明 "Compute the complex macroscopic dielectric function, at the RPA level, neglecting local field effects"。

### 7.1 第 1 步：带空带的 nscf

```fortran
! Non-SCF run on a dense full grid with empty bands, for epsilon.x (independent-particle optics).
! Run AFTER gaas.scf.in (same prefix/outdir). It overwrites the wavefunctions in gaas.save.
&CONTROL
  calculation = 'nscf'
  prefix      = 'gaas'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 2
  ecutwfc   = 80.0
  ecutrho   = 320.0
  nbnd      = 40        ! 14 occupied bands (28 valence electrons) + empty bands
  nosym     = .true.    ! full grid, no symmetry reduction (safe for every QE version)
  noinv     = .true.
/
&ELECTRONS
  conv_thr        = 1.0d-10
  diago_full_acc  = .true.   ! converge empty states as tightly as occupied ones
/
ATOMIC_SPECIES
  Ga  69.723  Ga.upf
  As  74.922  As.upf
CELL_PARAMETERS angstrom
  -2.8250   0.0000   2.8250
   0.0000   2.8250   2.8250
  -2.8250   2.8250   0.0000
ATOMIC_POSITIONS crystal
  Ga  0.00  0.00  0.00
  As  0.25  0.25  0.25
K_POINTS automatic
  12 12 12 0 0 0
```

为什么用 `nosym`：7.6 之前的 `epsilon.x` 没有对偶极矩阵元做对称化，必须在全布里渊区网格上计算；7.6 加入了对称化（见 §7.3）。为了让同一输入在各版本上都正确，本例保守地使用全网格。7.6 用户可以尝试去掉 `nosym`/`noinv` 并与全网格结果对比。

### 7.2 第 2 步：epsilon.x

```fortran
&INPUTPP
  prefix      = 'gaas'
  outdir      = './tmp'
  calculation = 'eps'
/
&ENERGY_GRID
  intersmear = 0.10     ! interband broadening (eV)
  wmin       = 0.0      ! eV
  wmax       = 15.0     ! eV
  nw         = 1500
/
```

`epsilon.x` 没有 INPUT_*.html 页面，说明文档是源码树中的 `PP/Doc/eps_man.tex`；下表按 QE 7.6 源码 `PP/src/epsilon.f90` 中的 namelist 与默认值整理：

| 关键词（namelist） | 默认值 | 说明 |
|---|---|---|
| `prefix`、`outdir`（INPUTPP） | `'pwscf'`、`$ESPRESSO_TMPDIR` 或 `./` | 与 nscf 一致 |
| `calculation`（INPUTPP） | `'eps'` | 还可选 `'jdos'`（联合态密度）、`'offdiag'`（非对角分量） |
| `intersmear`（ENERGY_GRID） | 0.136 | 带间展宽，eV |
| `intrasmear`（ENERGY_GRID） | 0.0 | 金属带内（Drude）项展宽，eV |
| `wmin` / `wmax` / `nw` | 0.0 / 30.0 / 600 | 频率网格（eV） |
| `nbndmin` / `nbndmax` | 1 / 全部能带 | 参与跃迁的能带范围 |
| `shift` | 0.0 | 剪刀算符（eV），把导带整体上移 |
| `smeartype` | `'gauss'` | jdos 使用的展宽类型（`'gauss'` 或 `'lorentz'`） |

```bash
mpirun -np 8 pw.x -nk 4 -in gaas.nscf_eps.in > gaas.nscf_eps.out
mpirun -np 8 epsilon.x -in gaas.eps.in > gaas.eps.out
ls epsr_gaas.dat epsi_gaas.dat eels_gaas.dat ieps_gaas.dat
head -3 epsi_gaas.dat
```

输出文件名为 `<量>_<prefix>.dat`：`epsr`（实部）、`epsi`（虚部）、`eels`（能量损失谱）、`ieps`（虚频轴上的 ε），每个文件的列依次为能量（eV）和 x、y、z 三个对角分量。如果 `nbnd` 不足（没有空带），程序报 `bad band number`。

成功判据：

- ε2（`epsi`）在 PBE/LDA 带隙以下接近 0，带隙以上出现吸收峰；
- 立方晶体三列相同；
- ε1（`epsr`）在 ω → 0 处的值与 `ph.x` 中 `lnoloc = .true.`（无局域场，ΔV_Hxc = 0）的介电常数接近，而不是与默认 `epsil` 结果比较——后者包含局域场与交换关联响应；
- 增加 `nbnd` 和 k 点后谱形不再明显变化。

### 7.3 QE 7.6 的变化

QE 7.6 release notes："Symmetry and the nonlocal pseudopotential contribution added to PP/epsilon for eps_calc, now allowing USPP/PAW pseudos. Calculations that previously used nonlocal pseudopotentials (ex: ONCVs) may yield different (more accurate) numerical results"。也就是说：

- 7.6 起 `epsilon.x` 在偶极矩阵元中加入了非局域赝势的对易子项，并支持对称化；
- 7.6 起可以用 US/PAW 赝势；
- 对同一模守恒（如 ONCV/PseudoDojo）输入，7.6 与旧版本的结果会有差别，新结果更准确。比较文献数据时要写清版本。

### 7.4 超越 IPA：turboTDDFT 与 Yambo

| 方法 | QE 生态中的工具 | 适用 | 入口 |
|---|---|---|---|
| 线性响应 TDDFT（Liouville–Lanczos） | turboTDDFT：`turbo_lanczos.x` + `turbo_spectrum.x`（另有 `turbo_davidson.x`、`turbo_eels.x`） | 分子吸收谱、含局域场的固体 EELS | [INPUT_Lanczos](https://www.quantum-espresso.org/Doc/INPUT_Lanczos.html)、[TDDFPT examples](https://github.com/QEF/q-e/tree/develop/TDDFPT/examples) |
| GW 准粒子能带 | Yambo（`p2y` 把 `prefix.save` 转成 Yambo 数据库） | 带隙、能带排列 | [Yambo wiki](https://www.yambo-code.eu/wiki/index.php/Main_Page) |
| BSE 激子光谱 | Yambo | 激子峰、二维材料光学 | 同上 |

turboTDDFT 的典型流程：pw.x scf（分子常用大盒子 + Γ 点）→ `turbo_lanczos.x` 迭代 Lanczos 系数 → `turbo_spectrum.x` 生成吸收谱。Yambo 流程：pw.x scf → nscf（大量空带）→ `p2y` → `yambo -i`（初始化）→ GW/BSE 输入。Yambo 官网的 [h-BN 教程](https://www.yambo-code.eu/wiki/index.php/Bulk_material:_h-BN) 给出了从 QE 开始的完整步骤，这里不再重复。

---

## 8. 与 VASP 的对应关系

从 VASP 切换过来的读者可以对照下表（VASP 一侧的详细讲解见 [VASP 第14章 光学性质与介电函数](../../VASP/tutorials/14-optics-dielectric.md) 和 [VASP 第12章 声子](../../VASP/tutorials/12-phonon.md)）：

| 任务 | VASP | Quantum ESPRESSO |
|---|---|---|
| ε∞ 与 Z*（DFPT） | `LEPSILON = .TRUE.` | `ph.x`：`epsil = .true.`（+ `trans` 得到 Z*） |
| ε∞ 与 Z*（有限电场） | `LCALCEPS = .TRUE.` | `pw.x`：`lelfield = .true.` + `efield_cart`（Berry 相位有限电场） |
| Γ 点力常数 | `IBRION = 7/8` | `ph.x` q = 0 |
| LO-TO 与 IR 强度 | 后处理（phonopy `--nac`、自写脚本） | `dynmat.x`（`q` 方向、`asr`） |
| 色散中的非解析修正 | phonopy BORN 文件 | `q2r.x` + `matdyn.x` 自动处理 |
| 拉曼 | 有限位移 + 有限电场（外部脚本） | `ph.x`：`lraman = .true.`（NC + LDA） |
| IPA 介电函数 | `LOPTICS = .TRUE.` + 大 `NBANDS` | nscf（大 `nbnd`）+ `epsilon.x` |
| GW/BSE | `ALGO = EVGW0`、`BSE` | Yambo（`p2y`）、BerkeleyGW |

---

## 参数选择与收敛建议

| 参数 | 建议 | 判据 |
|---|---|---|
| `ecutwfc` | 从赝势提示值起，每次 +10 Ry | ε∞ 变化 < 1%，光学模频率变化 < 2 cm⁻¹ |
| scf k 网格 | 8³ → 10³ → 12³（GaAs 类半导体） | ε∞ 收敛（比总能慢），Z* 求和接近 0 |
| `conv_thr` | 1e-12（晶体），1e-14（分子） | ph.x 无 "convergence not achieved" |
| `tr2_ph` | 1e-14 ~ 1e-16 | Z* 求和 < 0.01~0.05，声学模 < 1 cm⁻¹（施加 ASR 前也不应过大） |
| 晶格常数 | 同一泛函下优化 | 用实验晶格常数时要注明 |
| ph.x q 网格（色散） | 4³ 起，检查 6³ | 色散曲线不再变化，无虚频 |
| 分子盒子 | 10 → 12 → 15 Å | 内振动频率变化 < 5 cm⁻¹ |
| `epsilon.x` 的 `nbnd` | 占据带 2–3 倍起 | ε1(0) 与高能谱形收敛 |
| `epsilon.x` 的 k 网格 | 远密于 scf（12³、16³、20³） | 谱形不再出现锯齿 |
| `intersmear` | 0.05–0.2 eV | 只影响展宽，报告时注明 |

并行提示：`ph.x` 支持 `-nk`（pools）和 `-ni`（images，把 q 点和不可约表示分到不同 image）。Γ 点单个计算一般用 pools；q 网格计算可以用 images 或手动用 `start_q`/`last_q` 拆分（见第 12 章、第 18 章）。

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 处理 |
|---|---|---|
| `no elec. field with metals` | scf 使用了 smearing 或四面体占据，系统被当作金属 | 绝缘体去掉 `occupations='smearing'`；真正的金属不能算 ε∞ 和 Z* |
| `The phonon code with US-PP and raman or elop not yet available` | `lraman`/`elop` 用了超软赝势 | 换成模守恒 LDA 赝势（如 PseudoDojo PW 版） |
| `The phonon code with paw and raman or elop is not yet available` | `lraman`/`elop` 用了 PAW | 同上 |
| `lraman, elop, and noncolin not programmed` | 非共线/SOC 计算中开了 `lraman` | 不能在非共线下算拉曼张量 |
| `Raman for GPU not present in this version` | GPU 编译版本开了 `lraman` | 用 CPU 版本 ph.x |
| `raman: spin-polarized case not implemented` | `nspin = 2` | 拉曼只能非自旋极化 |
| `The phonon code with hybrid functionals is not yet available` | scf 用了 `input_dft='hse'` 等 | ph.x 不支持杂化泛函；pw.x 的有限电场（`lelfield`）与杂化泛函同用也会报 `hybrid XC and electric fields untested`，声子只能用有限位移（phonopy） |
| `The phonon code with noncollinear magnetism and electric field is not implemented` | 非共线磁性体系上设了 `epsil`/`zue`/`lraman`/`elop` | 不支持；非磁 SOC 体系可以算 ε∞ 和 Z* |
| `The phonon code with Raman and Hubbard U is not implemented` | DFT+U 基态 + `lraman` | 拉曼不支持 DFT+U |
| `The phonon code with meta-GGA functionals is not yet available` | scf 用了 SCAN/r2SCAN | 同上 |
| `The phonon code with TS-VdW is not yet available` | `vdw_corr='ts'` | 换 DFT-D2 或非局域 vdW-DF（ph.x 支持） |
| `The phonon code is not working after vc-relax` | 直接在 vc-relax 的 save 目录上跑 ph.x | 把优化后的晶胞写进新的 scf 输入，重跑 scf |
| `cannot start from pw.x data file using Gamma-point tricks` | 分子用了 `K_POINTS gamma` | 改为 `K_POINTS automatic` + `1 1 1 0 0 0` |
| `dynmat.x`：`A direction for q was not specified:TO-LO splitting will be absent` | `q` 未设置 | 要看 LO 模就在 `&INPUT` 中设置 `q(1)`、`q(2)`、`q(3)`；算 ε0（`lperm`）时则正是要不设 `q` |
| `dynmat.x`：`BEWARE: phonon contribution to permittivity computed with TO-LO splitting` | `lperm = .true.` 的同时设了非零 `q` | ε0 用不带 `q` 的单独一次 dynmat.x 计算 |
| `only dielectric constant with lrpa or lnoloc` | `lrpa`/`lnoloc` 与 `trans = .true.`（默认）同时使用 | 设 `trans = .false.`，单独计算介电常数 |
| `gamma is needed for elec.field` | 非 Γ 点（q ≠ 0）上设了 `epsil`/`zue`/`lraman`/`elop` | 电场相关量只能在 q = 0 计算 |
| `epsilon`：`bad band number` | nscf 中没有空带 | 增大 `nbnd` |
| Z* 求和明显不为零 | k 点、截断或 `tr2_ph` 不足 | 加密 k 点、收紧 `tr2_ph`；`dynmat.x`/`q2r.x` 施加 ASR 只是修正，不能替代收敛 |
| 色散在 Γ 附近光学支没有 LO-TO 劈裂 | Γ 点 dyn 文件里没有 Z*（体系被当成金属，或用了 `epsil=.false.`） | 检查 `gaas.dyn1` 中是否有 `Dielectric Tensor`，检查 scf 占据方式 |
| 分子转动模频率几十 cm⁻¹ | 盒子小、截断低、egg-box 效应 | `asr='zero-dim'`；增大盒子与截断 |
| 7.4.1–7.5 算得的拉曼强度异常 | 已知 bug（7.6 修复） | 升级到 QE 7.6 |
| 新旧版本 `epsilon.x` 结果不同 | 7.6 加入非局域项与对称化 | 正常现象，以 7.6 为准并注明版本 |

---

## 练习

1. 对 GaAs 把 scf k 网格从 6³ 增加到 12³，记录 ε∞ 和 Ga 的 Z* 的变化，画收敛曲线。
2. 在 `gaas.ph_epsil.in` 中同时打开 `zeu` 与 `zue`，比较两种算法得到的 Z*。
3. 另写两个 ph.x 输入（`trans = .false.`、去掉 `lraman`），分别用 `lrpa = .true.` 和 `lnoloc = .true.` 计算 ε∞，与默认结果比较，讨论局域场和交换关联响应的贡献；再把 `lnoloc` 的结果与 `epsilon.x` 的 ε1(0) 对比。
4. 用 `dynmat.x` 的 `lperm = .true.` 得到 ε0，用 LST 关系由 TO/LO 频率计算 ε0/ε∞ 并与之比较。
5. 构建岩盐结构 MgO（LDA NC 赝势），计算 ε∞、Z*、Γ 点频率，验证 IR/拉曼互斥规则。
6. 把 H2O 的盒子从 10 Å 增加到 15 Å，观察 3 个内振动频率和 IR 强度的变化。
7. 对 GaAs 的 `epsilon.x` 计算分别取 `nbnd = 20, 40, 60`，比较 ε2 谱形；再用 `shift` 加一个剪刀算符，观察吸收边的移动。

---

## 参考资料

**官方文档**

- [INPUT_PH（ph.x）](https://www.quantum-espresso.org/Doc/INPUT_PH.html) — `epsil`、`trans`、`zeu`/`zue`、`lraman`、`elop`、`lrpa`、`lnoloc`、`asr`、`nogg`
- [INPUT_DYNMAT（dynmat.x）](https://www.quantum-espresso.org/Doc/INPUT_DYNMAT.html) — `asr`、`q`、`lperm`、`lplasma`、`filmol`、`loto_2d`
- [INPUT_Q2R](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html) / [INPUT_MATDYN](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html) — `zasr`、`asr`、`q_in_band_form`
- [PHonon User's Guide](https://www.quantum-espresso.org/Doc/ph_user_guide/) — "Single-q calculation" 一节说明 Γ 点不含非解析项和 ASR，需要 `dynmat.x`
- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — `lelfield`、`efield_cart`、`nbnd`、`nosym`
- [INPUT_Lanczos（turbo_lanczos.x）](https://www.quantum-espresso.org/Doc/INPUT_Lanczos.html)
- [QE 7.6 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.6) — 拉曼 bug 修复与 `epsilon.x` 变化
- [PHonon/examples](https://github.com/QEF/q-e/tree/develop/PHonon/examples) — example05（AlAs 拉曼/电光张量）、example15（CO2、ZnO 的 IR/拉曼）、example10（非共线 Z*）
- [Yambo wiki](https://www.yambo-code.eu/wiki/index.php/Main_Page) — 从 QE 出发的 GW/BSE 教程

**课程与教程（只链接）**

- [PARADIM Summer School 2018（Cornell）](https://www.paradim.org/cu_ss_2018) — Giustino、Pamuk 主讲，含 Cl2 振动、金刚石声子、GaAs 的 LO-TO/IR/介电常数上机与视频；基于旧版 QE 与 LDA 赝势，原理部分非常适合本章
- [MSE404-MM（Imperial College London）](https://imperialcollegelondon.github.io/MSE404-MM/) — 声子实验课，讲 ph.x 基本流程
- [QuantumESPRESSO-school-2023](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023) — 高级学校上机材料，含 DFPT 部分
- [pranabdas/espresso](https://pranabdas.github.io/espresso/) — 安装到 EPW 的系统教程，含声子与光学
- [nguyen-group/QE-SSP](https://github.com/nguyen-group/QE-SSP) — 《QE Course for Solid-State Physics》配套算例，含声子、拉曼

**相关 GitHub 仓库（只链接）**

- [QEF/q-e](https://github.com/QEF/q-e) — 官方源码镜像（开发在 GitLab）
- [phonopy/phonopy](https://github.com/phonopy/phonopy) — 有限位移声子，可读 QE 输出，与 DFPT 对比；杂化泛函/meta-GGA 下的替代方案
- [yambo-code/yambo](https://github.com/yambo-code/yambo)、[yambo-code/yambopy](https://github.com/yambo-code/yambopy) — 基于 QE 的 GW/BSE 与后处理

**关键文献**

- S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, *Rev. Mod. Phys.* 73, 515 (2001)（DFPT 综述）. DOI: [10.1103/RevModPhys.73.515](https://doi.org/10.1103/RevModPhys.73.515)
- S. Baroni, R. Resta, *Phys. Rev. B* 33, 7017 (1986)（DFPT 宏观介电常数）. DOI: [10.1103/PhysRevB.33.7017](https://doi.org/10.1103/PhysRevB.33.7017)
- P. Giannozzi, S. de Gironcoli, P. Pavone, S. Baroni, *Phys. Rev. B* 43, 7231 (1991)（半导体声子色散、非解析项与插值）. DOI: [10.1103/PhysRevB.43.7231](https://doi.org/10.1103/PhysRevB.43.7231)
- X. Gonze, C. Lee, *Phys. Rev. B* 55, 10355 (1997)（Born 电荷、介电张量与力常数）. DOI: [10.1103/PhysRevB.55.10355](https://doi.org/10.1103/PhysRevB.55.10355)
- M. Lazzeri, F. Mauri, *Phys. Rev. Lett.* 90, 036401 (2003)（二阶响应拉曼张量，ph.x `lraman` 的方法）. DOI: [10.1103/PhysRevLett.90.036401](https://doi.org/10.1103/PhysRevLett.90.036401)
- O. B. Malcıoğlu, R. Gebauer, D. Rocca, S. Baroni, *Comput. Phys. Commun.* 182, 1744 (2011)（turboTDDFT）. DOI: [10.1016/j.cpc.2011.04.020](https://doi.org/10.1016/j.cpc.2011.04.020)
- D. Sangalli et al., *J. Phys.: Condens. Matter* 31, 325902 (2019)（Yambo）. DOI: [10.1088/1361-648X/ab15d0](https://doi.org/10.1088/1361-648X/ab15d0)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* 29, 465901 (2017)（QE 进展）. DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)

---

[← 上一章](12-phonon-dfpt.md) | [返回目录](../README.md) | [下一章 →](14-hybrid-metagga.md)
