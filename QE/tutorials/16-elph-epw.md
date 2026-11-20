# 第16章 电声耦合与超导：ph.x 与 EPW 入门（Electron-Phonon Coupling & Superconductivity: ph.x and EPW）

本章介绍 Quantum ESPRESSO 中计算电声耦合（electron-phonon coupling, EPC）的两条路线。第一条只用 QE 自带程序：`pw.x`（`la2F=.true.`）→ `ph.x`（`electron_phonon='interpolated'`、`el_ph_sigma`、`el_ph_nsigma`）→ `q2r.x`/`matdyn.x`（`la2F=.true.`）→ `lambda.x`，得到 Eliashberg 谱函数 α²F(ω)、耦合常数 λ、ω_log 和 Allen-Dynes/McMillan 公式给出的 T_c。第二条用 EPW（`epw.x`）做 Wannier 插值：粗网格 DFPT → `pp.py` → nscf → Wannier 化（`wannierize`、`nbndsub`、`proj`）→ 在细网格（`nkf`/`nqf`）上插值，求解 Migdal-Eliashberg 方程（`eliashberg`、`liso`、`limag`、`muc`），也可以继续算迁移率等输运性质。示例体系为 fcc Al（ph.x 路线）和 fcc Pb（EPW 路线）。

> [!NOTE]
> 本章所有输入在 [examples/16-elph-epw](../examples/16-elph-epw/README.md)。仓库**没有运行**这些计算，文中输出片段都标注为"（示意输出，数值仅供参考）"，只用来说明格式和检查方法。关键词核对的是官方 [INPUT_PH](https://www.quantum-espresso.org/Doc/INPUT_PH.html)、[INPUT_MATDYN](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html)、`PHonon/PH/lambda.f90` 源码注释和 [EPW v6.1 输入文档](https://docs.epw-code.org/Inputs/Inputs.html)。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 原理 | 电声矩阵元 g、声子线宽 γ_qν、λ_qν、α²F(ω)、λ、ω_log，以及 McMillan/Allen-Dynes 公式 |
| ph.x 路线 | 两次 SCF（稠密 k 网格写 `a2Fsave` + 声子网格），`ph.x` 同时算声子和 λ_qν，`q2r.x`/`matdyn.x` 插值线宽，`lambda.x` 汇总 |
| 读结果 | 在 `ph.x` 输出中找 `lambda(` 行，看 `matdyn` 写出的 `lambda` 文件和 `lambda.x` 的 T_c 表；学会用 smearing 平台判断收敛 |
| EPW 路线 | 粗网格 DFPT → `pp.py` 收集 dV_SCF → nscf（完整 k 列表）→ `epw.x` Wannier 化 → 细网格插值 → 各向同性 Eliashberg 方程 |
| 版本 | EPW 与 QE 联动发布（QE 7.5 ↔ EPW 6.0，QE 7.6 ↔ EPW 6.1），旧教程的关键字可能已经变化 |
| 成本 | 估计粗/细网格的计算量，知道哪些步骤可以并行、哪些不行 |

适用场景：

- **常规金属超导体**（Al、Pb、Nb、MgB₂、氢化物等）的 λ 与 T_c 估算；
- 声子线宽（声子寿命中的电声部分），对比非弹性中子/X 射线散射测得的线宽；
- 半导体载流子迁移率、电子自能、带隙重整化等（EPW 路线；本章只做概览）。

先修内容：DFPT 声子（[第12章](12-phonon-dfpt.md)）是本章的直接基础；EPW 路线还要用到 Wannier 函数（[第15章](15-wannier90.md)）。

与其它代码/工具的对应关系：

| 任务 | QE 生态 | 其它常见选择 | 说明 |
|---|---|---|---|
| 声子 | `ph.x`（DFPT） | VASP + phonopy（有限位移，见 [VASP 第12章](../../VASP/tutorials/12-phonon.md)） | 电声耦合需要 dV_SCF，DFPT 天然提供 |
| λ、α²F、T_c（粗糙但便宜） | `ph.x` + `lambda.x` | — | 依赖稠密 k 网格上的双 δ 函数求和，收敛较慢 |
| λ、Eliashberg 方程、各向异性能隙 | EPW | — | Wannier 插值到很密的 k/q 网格 |
| 迁移率（Boltzmann 方程） | EPW | Perturbo 等 | EPW 支持 SERTA 与迭代 BTE |
| 读取 QE/EPW 数据做模型 | [elphmod](https://github.com/janberges/elphmod) | — | Python 包，适合做电声模型和教学 |

---

## 2. 原理速览

### 2.1 电声矩阵元与声子线宽

在 Born-Oppenheimer 和谐近似下，原子按声子模 (q, ν) 振动时，自洽势的变化为 ∂_{qν}V。电声矩阵元为

$$
g_{mn\nu}(\mathbf k,\mathbf q)=\left(\frac{\hbar}{2M\omega_{\mathbf q\nu}}\right)^{1/2}\langle\psi_{m\mathbf k+\mathbf q}|\partial_{\mathbf q\nu}V|\psi_{n\mathbf k}\rangle .
$$

DFPT（[第12章](12-phonon-dfpt.md)）求解声子时已经得到 ∂_{qν}V（自洽势的线性响应 dV_SCF），所以只要把它与 Kohn-Sham 波函数做内积就能得到 g。这也是为什么 `ph.x` 做电声计算时必须设置 `fildvscf` 保存 dV_SCF（不设会报 `El-ph needs a DeltaVscf file`）。

金属中，声子可以把能量交给费米面附近的电子，由此产生的声子线宽（Fermi 黄金规则，双 δ 函数近似）为

$$
\gamma_{\mathbf q\nu}=2\pi\omega_{\mathbf q\nu}\sum_{mn}\int\frac{d\mathbf k}{\Omega_{BZ}}|g_{mn\nu}(\mathbf k,\mathbf q)|^2\,\delta(\epsilon_{n\mathbf k}-\epsilon_F)\,\delta(\epsilon_{m\mathbf k+\mathbf q}-\epsilon_F).
$$

每个模式的耦合强度为

$$
\lambda_{\mathbf q\nu}=\frac{\gamma_{\mathbf q\nu}}{\pi N(\epsilon_F)\,\omega_{\mathbf q\nu}^2},
$$

其中 N(ε_F) 是费米能处每自旋的电子态密度（与 `lambda.f90` 和 `elphon.f90` 源码注释中的定义一致）。

### 2.2 α²F(ω)、λ 与 ω_log

对布里渊区的 q 和所有声子支求和，得到 Eliashberg 谱函数

$$
\alpha^2F(\omega)=\frac{1}{2}\sum_{\nu}\int\frac{d\mathbf q}{\Omega_{BZ}}\,\omega_{\mathbf q\nu}\lambda_{\mathbf q\nu}\,\delta(\omega-\omega_{\mathbf q\nu}),
$$

以及

$$
\lambda=2\int_0^\infty\frac{\alpha^2F(\omega)}{\omega}d\omega=\sum_\nu\int\frac{d\mathbf q}{\Omega_{BZ}}\lambda_{\mathbf q\nu},\qquad
\omega_{\log}=\exp\left[\frac{2}{\lambda}\int_0^\infty\frac{d\omega}{\omega}\alpha^2F(\omega)\ln\omega\right].
$$

### 2.3 McMillan / Allen-Dynes 公式

`lambda.x` 用的是 Allen-Dynes 修正的 McMillan 公式（不含强耦合修正因子 f₁f₂）：

$$
T_c=\frac{\omega_{\log}}{1.2}\exp\left[-\frac{1.04(1+\lambda)}{\lambda-\mu^*(1+0.62\lambda)}\right].
$$

μ\* 是库仑赝势，通常取经验值 0.1–0.15，**不是第一性原理算出来的**。λ 较大（例如 λ > 1.5）时这个公式明显低估 T_c，应改为求解完整的 Migdal-Eliashberg 方程（EPW 的 `eliashberg`）。T_c 对 μ\* 的指数依赖意味着：不同 μ\* 给出的 T_c 可能相差一倍，报告结果时必须写明 μ\*。

### 2.4 为什么电声计算收敛这么慢

γ_qν 中有两个 δ 函数，而且被积函数只在费米面附近有贡献。粗 k 网格上用宽 smearing 代替 δ 函数会带来系统误差，用窄 smearing 又会产生很大的噪声。两条路线的应对办法不同：

- **ph.x 路线**（`electron_phonon='interpolated'`）：额外做一次稠密 k 网格的 SCF（`la2F=.true.`），把本征值写进 `a2Fsave` 文件；`ph.x` 在稠密网格上插值来计算双 δ 求和，并同时给出一组 smearing（`el_ph_sigma` 的整数倍，共 `el_ph_nsigma` 个）对应的结果，用户从中找平台。
- **EPW 路线**：先把电子哈密顿量、动力学矩阵和 g 都变换到最局域的 Wannier 表象（实空间），再 Fourier 插值到任意细的 k、q 网格，代价很低。细网格可以达到 50³–100³ 甚至更高，这是 EPW 能做定量 T_c 和迁移率计算的原因。

> [!NOTE]
> 若要深入理论（Migdal-Eliashberg 方程、各向异性能隙、迁移率），推荐 F. Giustino 的综述 *Rev. Mod. Phys.* **89**, 015003 (2017)，以及 EPW 官网的 [Theory](https://docs.epw-code.org/Theory.html) 页面。

---

## 3. 两条路线怎么选

| 比较项 | ph.x + `lambda.x` | EPW |
|---|---|---|
| 需要的程序 | `pw.x`、`ph.x`、`q2r.x`、`matdyn.x`、`lambda.x`（随 PHonon 编译） | 另外编译 `epw.x`（`make epw`，会同时编译 wannier90） |
| 费米面积分 | 稠密 k 网格（`a2Fsave`）上插值 | Wannier 插值，细网格可以非常密 |
| q 网格 | 每个 q 都要做一次完整 DFPT；α²F 的 q 积分靠 `matdyn.x` 插值线宽 | 粗 q 网格 DFPT，然后插值到细 q 网格 |
| 结果 | λ、α²F、ω_log、Allen-Dynes T_c、声子线宽 | 以上全部，外加各向同性/各向异性 Eliashberg 方程、能隙、迁移率、自能、极化子等 |
| 适合 | 简单金属、快速估算、教学 | 定量研究、多带/各向异性超导体（如 MgB₂）、半导体输运 |
| 主要难点 | smearing 平台不明显；q 网格收敛贵 | Wannier 化质量（冻结窗口、投影）；内存和磁盘 |

建议的学习顺序：先用 Al 跑通 ph.x 路线，弄清 λ、α²F 和 smearing 平台的含义，再学 EPW。

---

## 4. 实战一：fcc Al 的 ph.x 电声耦合

### 4.1 工作流

```text
al.scf.fit.in  (k 24x24x24, la2F=.true.)  --> tmp/al.save + tmp/al.a2Fsave
       |
al.scf.in      (k 12x12x12, same prefix/outdir)
       |
al.elph.in     ph.x: ldisp 4x4x4, fildvscf, electron_phonon='interpolated'
       |        --> al.dyn0..N, elph_dir/elph.inp_lambda.*, elph_dir/a2Fq2r.*
       |
al.q2r.in      q2r.x (la2F=.true.)        --> al444.fc + elph_dir/a2Fmatdyn.*
       |
al.matdyn.*.in matdyn.x (la2F=.true.)     --> dispersion + linewidths; a2F.dos*, lambda
       |
al.lambda.in   lambda.x                    --> lambda, omega_log, Tc (Allen-Dynes)
```

这个流程与官方 `PHonon/examples/example03`（`run_example` 中的 Al 电声算例）一致。官方算例使用旧的 LDA `Al.pz-vbc.UPF` 和很低的截断能，只为演示流程；本仓库的输入改用 PseudoDojo PBE 模守恒赝势，并加密了 k 网格。

### 4.2 第一步：稠密 k 网格 SCF（`la2F=.true.`）

```fortran
! Step 1 of the ph.x electron-phonon route for fcc Al.
! Pseudopotential: PseudoDojo nc-sr-05 PBE standard, Al.upf (norm-conserving).
&CONTROL
  calculation = 'scf'
  prefix      = 'al'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav       = 2
  celldm(1)   = 7.63          ! Bohr, about 4.04 Angstrom; use your relaxed value
  nat         = 1
  ntyp        = 1
  ecutwfc     = 40.0          ! Ry; check against the PseudoDojo hints and chapter 4
  occupations = 'smearing'
  smearing    = 'mv'          ! Marzari-Vanderbilt cold smearing
  degauss     = 0.02          ! Ry
  la2F        = .true.        ! write a2Fsave for the el-ph Fermi-surface sum
/
&ELECTRONS
  conv_thr    = 1.0d-12
  mixing_beta = 0.7
/
ATOMIC_SPECIES
  Al  26.9815  Al.upf
ATOMIC_POSITIONS crystal
  Al  0.0  0.0  0.0
K_POINTS automatic
  24 24 24  0 0 0
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `la2F` | `.true.`（仅这一步） | 在 `outdir` 写出 `al.a2Fsave`（文件名为 `prefix.a2Fsave`；内容为能带数、k 点数、本征值、k 点、权重、网格与对称操作），供 `ph.x` 做插值。这个关键词 pw.x 能识别（见 `Modules/read_namelists.f90`），但 7.5 版 INPUT_PW 页面没有列出 |
| `occupations` | `'smearing'` | 电声计算只对金属有意义；`ph.x` 对非金属会报 `Electron-phonon only for metals` |
| k 网格 | 远密于声子用的网格 | 必须与后面的 q 网格相容（24 是 4 的倍数） |
| `conv_thr` | 1e-12 或更小 | DFPT 需要非常精确的基态 |

### 4.3 第二步：声子用的 SCF

`al.scf.in` 与上一步只差 `K_POINTS 12 12 12 0 0 0`，并去掉 `la2F`。两步的 `prefix='al'`、`outdir='./tmp'` 必须相同：第二步覆盖 `al.save` 中的波函数和电荷密度，`tmp/al.a2Fsave` 保留下来给 `ph.x` 用。**顺序不能颠倒**。另外两步的能带数 `nbnd` 必须相同（本例都用默认值）：`ph.x` 读 `a2Fsave` 时会检查，不一致时报 `Error in routine elphsum: wrong file read`。

### 4.4 第三步：ph.x 同时计算声子与 λ_qν

```text
Al: phonons + electron-phonon coefficients on a 4x4x4 q grid
&INPUTPH
  prefix          = 'al'
  outdir          = './tmp'
  fildyn          = 'al.dyn'
  fildvscf        = 'aldv'           ! dV_SCF is required for el-ph
  tr2_ph          = 1.0d-15
  ldisp           = .true.
  nq1 = 4, nq2 = 4, nq3 = 4
  trans           = .true.           ! phonons and lambdas in the same run
  electron_phonon = 'interpolated'   ! needs a2Fsave from al.scf.fit.in
  el_ph_sigma     = 0.005            ! Ry, spacing of double-delta smearings
  el_ph_nsigma    = 10               ! 10 smearings: 0.005 ... 0.050 Ry
/
```

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `electron_phonon` | `' '` | `'interpolated'` | 其它选项：`'simple'`（在 `nk1..3` 网格上直接求和）、`'lambda_tetra'`/`'gamma_tetra'`（优化四面体法，scf 需用 `occupations='tetrahedra_opt'`，否则报 `"electron_phonon" and "occupation" are inconsistent.`）、`'epa'`、`'ahc'`（电子自能，供 `postahc.x`，要求 `trans=.false.`）。`'simple'`/`'interpolated'` 只用于金属并要求 smearing |
| `el_ph_sigma` | 0.02 | 0.005 | 双 δ 函数展宽的**间隔**（Ry），实际使用 σ = `el_ph_sigma` × 1, 2, …, `el_ph_nsigma` |
| `el_ph_nsigma` | 10 | 10 | 展宽个数 |
| `fildvscf` | `' '` | 任意文件名 | 保存 dV_SCF；电声计算必须设置 |
| `trans` | `.true.` | `.true.` | `.true.`：同一次运行中算声子和 λ（k 网格相同）；`.false.`：读之前存下的 dV_SCF 和动力学矩阵，可以换更大的 k 网格 |
| `tr2_ph` | 1e-12 | 1e-14–1e-16 | 声子自洽阈值；电声计算建议比默认值更严 |
| `ldisp`、`nq1..3` | — | 4×4×4 起步 | q 网格要与 k 网格相容 |

运行（电声计算**不能用 `-ni` image 并行**，源码 `phq_readin.f90` 中会报 `el-ph with images not implemented`）：

```bash
mpirun -np 16 ph.x -nk 4 -in al.elph.in > al.elph.out
```

每个 q 点结束后，输出中会有每个 smearing 的态密度和每个模式的 λ、γ（输出格式见 `PHonon/PH/elphon.f90`；示意输出，数值仅供参考）：

```text
     Gaussian Broadening:   0.005 Ry, ngauss=   0
     DOS =  x.xxxxxx states/spin/Ry/Unit Cell at Ef=  x.xxxxxx eV
     lambda(    1)=  0.xxxx   gamma=    x.xx GHz
     lambda(    2)=  0.xxxx   gamma=    x.xx GHz
     lambda(    3)=  0.xxxx   gamma=   xx.xx GHz
     Gaussian Broadening:   0.010 Ry, ngauss=   0
     ...
```

同时在工作目录生成：

| 文件 | 内容 | 谁来读 |
|---|---|---|
| `al.dyn0` | q 网格与不可约 q 点列表 | `q2r.x` |
| `al.dyn1` … `al.dynN` | 每个不可约 q 的动力学矩阵（含星中所有 q） | `q2r.x` |
| `elph_dir/elph.inp_lambda.N` | 第 N 个不可约 q 的 ω²、各 σ 下的 N(ε_F)、λ_qν、γ_qν | `lambda.x` |
| `elph_dir/a2Fq2r.<50+i>.N` | 第 i 个 σ 的线宽（i = 1…`el_ph_nsigma`，即 `a2Fq2r.51.N` … `a2Fq2r.60.N`），供 `q2r.x` 做实空间变换 | `q2r.x`（`la2F=.true.`） |

检查：

```bash
grep -c "lambda(" al.elph.out                  # = N_q(irr) x el_ph_nsigma x 3 modes
ls al.dyn[0-9]* | wc -l                       # = N_q(irr) + 1 (al.dyn0 is the q list)
ls elph_dir/
grep "JOB DONE" al.elph.out
```

> [!TIP]
> Γ 点三个声学模的 λ 和 γ 应为 0 或很小的数。若 Γ 点声学支频率明显不为 0（几十 cm⁻¹ 以上），先回到 [第12章](12-phonon-dfpt.md) 检查 `tr2_ph`、`conv_thr` 和截断能。

### 4.5 第四步：q2r.x 与 matdyn.x（`la2F=.true.`）

```fortran
&INPUT
  fildyn = 'al.dyn'
  flfrc  = 'al444.fc'
  zasr   = 'simple'
  la2F   = .true.        ! also Fourier-transform the el-ph linewidths
/
```

`q2r.x` 的 namelist 中有 `la2F`（见 `PHonon/PH/q2r.f90`），设为 `.true.` 时会把每个 σ 的线宽也变换到实空间，写入 `elph_dir/a2Fmatdyn.*`。

声子 DOS 和 α²F：

```fortran
&INPUT
  asr      = 'simple'
  flfrc    = 'al444.fc'
  flfrq    = 'al.dos.freq'
  la2F     = .true.       ! writes a2F.dos1 ... a2F.dos10 and the file "lambda"
  dos      = .true.
  fldos    = 'al.phdos'
  nk1 = 24, nk2 = 24, nk3 = 24
  ndos     = 200
/
```

| 关键词（matdyn.x） | 推荐值 | 说明 |
|---|---|---|
| `la2F` | `.true.` | 同时插值电声系数（INPUT_MATDYN：*if .true. interpolates also the el-ph coefficients*） |
| `dos` | `.true.`（DOS/α²F）或 `.false.`（色散 + 线宽） | `dos=.true.` 时写 `a2F.dosN` 和 `lambda`；`dos=.false.` 时写 `gam.lines`（沿路径的线宽，单位 THz） |
| `nk1..3` | ≥ 20 | 插值 α²F 用的 q 网格；越密 α²F 越平滑 |
| `asr` | `'simple'` 或 `'crystal'` | 声学求和规则 |

`matdyn.x` 写出的文件（以 `matdyn.f90` 源码为准）：

- `a2F.dos1` … `a2F.dos10`：每个 σ 一个文件，第 1 列为频率（Ry），第 2 列为总 α²F，后面是各模式贡献，文件末尾一行给出积分得到的 λ；
- `lambda`：每个 σ 一行，`Broadening … lambda … dos(Ef) … omega_ln [K] …`。

```bash
cat lambda
```

（示意输出，数值仅供参考）

```text
 Electron-phonon coupling constant, lambda

 Broadening   0.0050 lambda       0.xxxx dos(Ef)  x.xxxx omega_ln [K]    xxx.xxxx
 Broadening   0.0100 lambda       0.xxxx dos(Ef)  x.xxxx omega_ln [K]    xxx.xxxx
 ...
```

### 4.6 第五步：lambda.x 计算 T_c

`lambda.x` 从标准输入读入（格式来自 `PHonon/PH/lambda.f90` 的注释）：

```text
emax  degaussq  ngaussq        ! emax (THz) above the highest phonon; q-sum smearing (THz); smearing type
nks                            ! number of irreducible q points
q(1) q(2) q(3) wk              ! one line per q: coordinates and weight (number of q in the star)
...
filelph(1)                     ! elph_dir/elph.inp_lambda.1
...
mustar                         ! mu* for the Allen-Dynes formula
```

本例的完整输入见 [al.lambda.in](../examples/16-elph-epw/phx-al/al.lambda.in)。准备时要注意：

1. q 点个数等于不可约 q 点数，即 `al.dyn1` … `al.dynN` 的个数（`head -2 al.dyn0` 的第二行）；
2. 权重是每个 q 的"星"中 q 点数，可以数对应 `al.dynN` 中 `Dynamical  Matrix in cartesian axes` 出现的次数（fcc 有反演对称，−q 总在星内）；`lambda.x` 会自动归一化，所以不必除以总数；
3. 顺序必须与 `elph.inp_lambda.1 … N` 一致；当前版本的 `lambda.x` 不检查 q 坐标是否一致（这段检查在源码里被注释掉了），写错了也不会报错，只会算错；
4. `emax` 要大于最高声子频率（THz）。Al 的最高声子频率在 10 THz 左右，本例取 12 THz。

```bash
lambda.x < al.lambda.in > al.lambda.out
tail -15 al.lambda.out
```

（示意输出，数值仅供参考）

```text
     lambda = 0.xxxxxx (   0.xxxxxx )  <log w>=  xxx.xxx K  N(Ef)=  x.xxxxxx at degauss= 0.005
     lambda = 0.xxxxxx (   0.xxxxxx )  <log w>=  xxx.xxx K  N(Ef)=  x.xxxxxx at degauss= 0.010
     ...
lambda        omega_log          T_c
   0.xxxxx       xxx.xxx            x.xxx
   0.xxxxx       xxx.xxx            x.xxx
   ...
```

同时生成 `lambda.dat`（每个 σ 的 λ、ω_log、N(ε_F)）和 `alpha2F.dat`（第 1 列频率 THz，其余各列是不同 σ 下的 α²F）。

### 4.7 如何判断结果可信

- **smearing 平台**：把 λ 对 σ 作图。σ 很小时噪声大，σ 很大时 λ 被平滑得偏离；中间若有一段 λ 基本不变，就取这段的值。如果没有平台，说明 `a2Fsave` 的 k 网格不够密，需要加密 `al.scf.fit.in` 的 k 网格（例如 24³ → 32³），或者加密 q 网格。
- **两种 λ 的一致性**：`lambda.x` 括号内的值是由 α²F 积分得到的 λ，应与括号外直接对 q 求和得到的 λ 接近。
- **与 `matdyn` 的 `lambda` 文件对比**：两者的 q 积分方式不同（直接求和 vs. 插值到稠密 q 网格），差别反映 q 网格的收敛程度。
- **与文献/实验比较**：Al 的 λ 文献值约 0.4，实验 T_c 约 1.2 K。T_c 对 λ 和 μ\* 都呈指数依赖，量级对就可以了，不要要求精确吻合。

---

## 5. 实战二：fcc Pb 的 EPW 入门

### 5.1 EPW 的版本与安装

- EPW 随 QE 源码一起发布，版本联动：QE 7.5 对应 EPW 6.0，QE 7.6 对应 EPW 6.1（GPU、image + pool 并行、DFPT+U、LSDA 等新功能）。EPW 官网也提供单独打包的源码（如 `q-e-EPW-6.1.tar.gz`）。
- 编译：在 QE 源码目录执行 `./configure` 后运行 `make epw`，会同时编译 `pw.x`、`ph.x` 和 wannier90；可执行文件在 `bin/epw.x`（或 `EPW/bin/epw.x`）。测试：`cd test-suite && make run-tests-epw`。
- 输入关键字以 [EPW 输入文档](https://docs.epw-code.org/Inputs/Inputs.html) 为准。网上很多 EPW 4.x/5.x 时代的教程已经过时，照抄可能报未知关键字，或者含义已经改变（例如 `vme` 在 v5.4 之前是逻辑型，现在是字符串 `'wannier'`/`'dipole'`）。

### 5.2 工作流

```text
pb.scf.in   (k 8x8x8)            pw.x
pb.ph.in    (q 6x6x6, fildvscf)  ph.x       --> pb.dyn*, _ph0/ (dvscf, patterns)
pp.py       (prefix "pb", run in the ph.x dir) --> save/ (dyn, dvscf, patterns)
pb.q2r.in                        q2r.x       (optional: check phonon dispersion)
pb.nscf.in  (full 6x6x6 list, nbnd=10)  pw.x
pb.epw1.in  wannierize, epwwrite     epw.x   --> pb.epmatwp, pb.ukk, *.fmt
pb.epw2.in  epwread, eliashberg      epw.x   --> pb.a2f, pb.imag_iso_*, pb.pade_iso_*
```

参数参照 EPW 官方教程 "Superconducting properties" 中的 Pb 练习，并做了缩减：细网格减小、温度点减少、不做迭代解析延拓（`lacon`），便于入门试跑。

### 5.3 粗网格 SCF 与 DFPT

SCF 输入见 [pb.scf.in](../examples/16-elph-epw/epw-pb/pb.scf.in)（Pb，PseudoDojo 模守恒赝势，`ecutwfc=60`，`smearing='mp'`，`degauss=0.04`，8×8×8 k）。DFPT 输入：

```text
Pb: DFPT phonons on a coarse 6x6x6 q grid, saving dV_SCF for EPW
&INPUTPH
  prefix    = 'pb'
  outdir    = './'           ! pp.py looks for ./_ph0 and ./pb.save
  fildyn    = 'pb.dyn'
  fildvscf  = 'dvscf'        ! EPW reads these via pp.py -> save/
  tr2_ph    = 1.0d-15
  ldisp     = .true.
  nq1 = 6, nq2 = 6, nq3 = 6  ! coarse q grid; must divide the coarse k grid of the nscf run
/
```

这里**不需要**设置 `electron_phonon`：EPW 自己根据 dV_SCF 和 nscf 波函数计算 g。

运行 `ph.x` 后，用 EPW 自带的脚本收集文件。`pp.py`（见 `EPW/bin/pp.py` 源码）在**当前目录**里同时寻找 `pb.dyn*`、`_ph0/pb.phsave/`、`_ph0/pb.dvscf*` 和 `pb.save/data-file-schema.xml`，所以 Pb 例子的所有输入都用 `outdir='./'`（与 EPW 官方教程相同）。如果把 `outdir` 设成别的目录，`pp.py` 会因找不到 `_ph0/` 或 `pb.save/` 而失败。

```bash
python3 /path/to/q-e/EPW/bin/pp.py     # asks for the prefix: type  pb
ls save/                               # pb.dyn_q*, pb.dvscf_q*, pb.phsave/ (patterns)
```

### 5.4 nscf：完整（不约化）的粗 k 网格

EPW 需要粗 k 网格上**全部** k 点的波函数（不能用对称性约化），通常用 Wannier90 自带的 `kmesh.pl` 生成 `K_POINTS crystal` 列表：

```bash
cat pb.nscf.in.head > pb.nscf.in
perl /path/to/wannier90/utility/kmesh.pl 6 6 6 >> pb.nscf.in    # 216 points with weights
grep -c "" pb.nscf.in
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `calculation` | `'nscf'` | 与 scf 相同的 `prefix`、`outdir` |
| `nbnd` | 10 | 要覆盖 Wannier 化需要的所有能带（5 条 5d + 4 条 sp + 余量） |
| 粗 k 网格 | 6×6×6 | 必须是粗 q 网格的整数倍；EPW 官方教程建议 k、q 粗网格取相同 |

### 5.5 EPW 第 1 步：Wannier 化

```fortran
&INPUTEPW
  prefix        = 'pb'
  outdir        = './'
  amass(1)      = 207.2
  dvscf_dir     = './save'          ! created by EPW/bin/pp.py

  elph          = .true.
  epwwrite      = .true.
  epwread       = .false.
  etf_mem       = 0
  asr_typ       = 'crystal'

  wannierize    = .true.
  nbndsub       = 4
  bands_skipped = 'exclude_bands = 1:5'   ! drop the five 5d bands
  num_iter      = 500
  proj(1)       = 'Pb:sp3'
  dis_froz_min  = -1.0             ! eV; frozen window, check against E_F in pb.scf.out
  dis_froz_max  = 16.5             ! eV; placeholder, adapt to your band structure
  wdata(1)      = 'dis_num_iter = 1000'
  wdata(2)      = 'use_ws_distance = .true.'

  nk1 = 6, nk2 = 6, nk3 = 6         ! coarse k grid (= nscf grid)
  nq1 = 6, nq2 = 6, nq3 = 6         ! coarse q grid (= ph.x grid)
  nkf1 = 1, nkf2 = 1, nkf3 = 1
  nqf1 = 1, nqf2 = 1, nqf3 = 1
/
```

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `elph` | `.false.` | 计算电声系数 |
| `epwwrite` / `epwread` | `.true.` / `.false.` | 写出 / 读入 Wannier 表象的电声矩阵元（`epwdata.fmt`、`*.epmatwp` 等）。`epwread=.true.` 需要先跑过一次 `epwwrite=.true.`，并同时设置 `kmaps=.true.` |
| `wannierize` | `.false.` | 调用 Wannier90 库计算 Wannier 函数并写出 `filukk`；为 `.false.` 时从磁盘读 |
| `nbndsub` | 0 | Wannier 函数个数 |
| `bands_skipped` | `''` | 不参与 Wannier 化的能带，例如 `'exclude_bands = 1:5'` |
| `proj(:)` | `''` | 初始投影，语法同 Wannier90 |
| `dis_win_min/max`、`dis_froz_min/max` | ±9999 | 解纠缠外窗口与冻结窗口（eV），含义同 Wannier90 |
| `wdata(:)` | `''` | 原样写入 `prefix.win` 的额外 Wannier90 关键字 |
| `etf_mem` | 1 | 内存/IO 的取舍；0 最快但最耗内存，`epw_memdist=.true.` 只能与 0 一起用 |
| `asr_typ` | `'simple'` | 实空间声学求和规则 |
| `nk1..3`、`nq1..3` | 0 | 粗 k、q 网格，必须与 nscf 和 ph.x 一致 |
| `nkf1..3`、`nqf1..3` | 0 | 细网格；Wannier 化这一步可以都设为 1 |

检查：

```bash
grep -A5 "Wannier Function centers" pb.epw1.out
grep -i "spread" pb.wout | tail -3
ls pb.epmatwp pb.ukk crystal.fmt epwdata.fmt wigner.fmt
```

（示意输出，数值仅供参考）

```text
     Wannier Function centers (cartesian, alat) and spreads (ang):
     (   x.xxxxx   x.xxxxx   x.xxxxx) :   x.xxxxx
     ...
```

4 个 sp3 型 Wannier 函数的 spread 应彼此相等（对称性要求），中心应围绕 Pb 原子对称分布。Wannier 插值能带必须在冻结窗口内与 DFT 能带重合；做法是在 `wdata` 里加 `bands_plot = .true.` 和 `begin kpoint_path … end kpoint_path`，与 [第15章](15-wannier90.md) 中 Wannier90 的写法相同。

### 5.6 EPW 第 2 步：细网格插值与各向同性 Eliashberg 方程

```fortran
&INPUTEPW
  prefix        = 'pb'
  outdir        = './'
  amass(1)      = 207.2
  dvscf_dir     = './save'

  elph          = .true.
  epwwrite      = .false.
  epwread       = .true.
  kmaps         = .true.            ! required together with epwread for a restart
  etf_mem       = 0
  asr_typ       = 'crystal'

  wannierize    = .false.
  nbndsub       = 4
  bands_skipped = 'exclude_bands = 1:5'

  fsthick       = 0.2               ! eV, Fermi-surface window
  degaussw      = 0.05              ! eV, electron smearing for the delta functions
  degaussq      = 0.15              ! meV, phonon smearing for the q sum

  eliashberg    = .true.
  a2f_iso       = .true.
  liso          = .true.            ! isotropic Migdal-Eliashberg
  limag         = .true.            ! imaginary (Matsubara) axis
  lpade         = .true.            ! Pade continuation to the real axis
  nsiter        = 500
  conv_thr_iaxis = 1.0d-4
  wscut         = 0.1               ! eV, Matsubara frequency cutoff
  muc           = 0.1               ! Coulomb pseudopotential mu*
  temps         = 0.5 1.0 2.0 3.0 4.0 5.0 6.0   ! K

  nk1 = 6, nk2 = 6, nk3 = 6
  nq1 = 6, nq2 = 6, nq3 = 6
  mp_mesh_k     = .true.
  nkf1 = 24, nkf2 = 24, nkf3 = 24   ! fine k grid: converge (e.g. 24 -> 36 -> 48)
  nqf1 = 12, nqf2 = 12, nqf3 = 12   ! fine q grid: nkf must be a multiple of nqf
/
```

| 关键词 | 默认值 | 推荐值/说明 |
|---|---|---|
| `fsthick` | 1e10 eV | 0.2–0.5 eV。只保留费米能附近窗口内的态，窗口越窄越省 |
| `degaussw` | 0.025 eV | 能量守恒 δ 函数的电子展宽；与细 k 网格一起收敛 |
| `degaussq` | 0.05 meV | q 求和中的声子展宽 |
| `eliashberg` | `.false.` | 求解 Eliashberg 方程和/或计算谱函数 |
| `a2f_iso` | `.false.` | 直接在线计算各向同性 α²F，不写电声矩阵元平方（g 的模平方）文件；只能与 `ephwrite=.false.`、`laniso=.false.`、`liso=.true.` 同用 |
| `liso` / `laniso` | `.false.` | 各向同性 / 各向异性 Eliashberg 方程。`laniso` 需要先用 `ephwrite=.true.` 写出 `prefix.ephmat` 目录 |
| `limag` | `.false.` | 在虚轴（Matsubara 频率）上求解 |
| `lpade`、`lacon` | `.false.` | Padé 近似 / 迭代解析延拓到实轴；`lacon` 需要 `limag` 和 `lpade` 都为 `.true.` |
| `muc` | 0.0 | μ\*；默认值为 0，一定要显式设置 |
| `wscut` | 1.0 eV | Matsubara 频率截断；`nswi > 0` 时被忽略 |
| `temps`、`nstemp` | 300 K | 温度列表（K）；最多 50 个 |
| `nsiter` | 40 | Eliashberg 自洽迭代上限 |
| `conv_thr_iaxis` | 1e-5 | 虚轴方程收敛阈值 |
| `mp_mesh_k` | `.false.` | 细 k 网格只取不可约部分 |
| `nkf1..3`、`nqf1..3` | 0 | 细网格；`ephwrite` 要求 `nkf` 是 `nqf` 的整数倍 |

运行：

```bash
mpirun -np 32 epw.x -nk 32 -in pb.epw2.in > pb.epw2.out
grep "Estimated Tc using McMillan" pb.epw2.out
ls pb.a2f pb.a2f_proj pb.imag_iso_* pb.pade_iso_*
```

（示意输出，数值仅供参考）

```text
     Estimated Tc using McMillan expression =   x.xxxx K for muc =   0.1000
```

输出文件（EPW 官方教程中的说明）：

| 文件 | 内容 |
|---|---|
| `pb.a2f` | 第 1 列频率（meV），第 2 列 α²F，第 3 列累积 λ(ω)；不同 smearing 依次排列 |
| `pb.a2f_proj` | 按原子投影的 α²F |
| `pb.imag_iso_XXX.XX` | 每个温度一个文件：Matsubara 频率、重整化函数 Z、能隙 Δ（eV） |
| `pb.pade_iso_XXX.XX` | Padé 延拓到实轴的 Z、Δ（实部、虚部） |

T_c 的读法：每个温度取 `pb.imag_iso_*` 中最低 Matsubara 频率处的 Δ，Δ(T) 降到 0 的温度就是 Eliashberg T_c。也可以设 `tc_linear=.true.` 求解线性化的 Eliashberg 方程直接得到 T_c（参见 EPW 官方教程）。

> [!WARNING]
> 示例中的细网格（24³/12³）只用于试跑。EPW 官方 Pb 练习用的是 48³/24³ 细网格，正式研究必须系统加密 `nkf`/`nqf`，同时测试 `degaussw`、`fsthick` 和粗网格，直到 λ 和 T_c 收敛。Pb 的自旋轨道耦合对声子和电声耦合也有可观影响，本例为了简单忽略了 SOC。

### 5.7 EPW 能做的其它事情（链接为主）

| 任务 | 主要关键字（EPW 6.1 文档） | 官方教程 |
|---|---|---|
| 各向异性超导能隙（MgB₂） | `ephwrite`、`laniso`、`limag`、`fbw` | [Tutorial 4](https://docs.epw-code.org/tutorials/tutorial_04/index.html) |
| 电子/声子能带插值与电声矩阵元检查 | `band_plot`、`filkf`、`filqf`、`prtgkk` | [Tutorial 1](https://docs.epw-code.org/tutorials/tutorial_01/index.html) |
| 迁移率（SERTA / 迭代 BTE） | `scattering`、`scattering_serta`、`iterative_bte`、`carrier`、`ncarrier`、`int_mob`、`temps` | [Tutorial 7](https://docs.epw-code.org/tutorials/tutorial_07/index.html) |
| 极化子 | 见文档第 7 节 | [Tutorial 5](https://docs.epw-code.org/tutorials/tutorial_05/index.html) |
| 声子辅助光吸收 | 见文档 | [Tutorial 6](https://docs.epw-code.org/tutorials/tutorial_06/index.html) |
| 特殊位移法（ZG，温度对能带的影响） | `ZG.x` | [Tutorial 2](https://docs.epw-code.org/tutorials/tutorial_02/index.html) |

极性材料（GaAs、GaN 等）需要考虑长程 Fröhlich 相互作用（`lpolar=.true.`），并在 ph.x 中先算 Born 有效电荷和介电常数（[第13章](13-dielectric-ir-raman.md)）。

---

## 6. 计算量与并行

| 步骤 | 主要耗时 | 并行建议 |
|---|---|---|
| 稠密 k 网格 SCF（`la2F`） | k 点多 | `pw.x -nk`（pool 并行），见 [第18章](18-automation-parallel-troubleshooting.md) |
| `ph.x` 电声 | 每个不可约 q 一次 DFPT | 只能 `-nk`/`-nd`，不能 `-ni`；可以用 `start_q`/`last_q` 把 q 点拆成多个作业（[第12章](12-phonon-dfpt.md)），但都要能访问同一个 `outdir` |
| `pp.py` | 文件复制 | 串行 |
| nscf（完整 k 列表） | 无对称约化 | `-nk` |
| `epw.x` | Wannier 化 + 细网格插值 | `-nk`（k 点 pool）和 `-ni`（q 点 image），nprocs = npool × nimage；粗网格阶段 image 数最多等于不可约 q 点数 |
| `epw.x` GPU（EPW 6.1） | — | 通常每块 GPU 一个 MPI 进程 |

磁盘与内存：

- `ph.x` 的 `_ph0/` 目录保存每个 q 的 dV_SCF，三维体系的大小与 FFT 网格和原子数成正比；
- EPW 的 `prefix.epmatwp` 大小与 Wannier 函数个数的平方、原子数和实空间格点数成正比；
- `ephwrite=.true.` 写出的 `prefix.ephmat/` 在细网格很密时可能达到几十 GB，只做各向同性 Eliashberg 时可以改用 `a2f_iso=.true.`。

---

## 参数选择与收敛建议

| 参数 | 起步值 | 收敛方法 |
|---|---|---|
| `ecutwfc` | 赝势 hint 的高值 | 声子频率变化 < 1–2% |
| `conv_thr`（scf） | ≤ 1e-12 | DFPT 的前提 |
| `tr2_ph` | 1e-15 | Γ 点声学支频率接近 0 |
| q 网格（ph.x） | 简单金属 4³–6³ | 依次加密，看 λ、ω_log 变化 |
| `a2Fsave` k 网格 | 声子 k 网格的 2–3 倍 | 出现 smearing 平台；λ 变化 < 5% |
| `el_ph_sigma` × `el_ph_nsigma` | 0.005 × 10 | 找平台，不要只取单个 σ |
| 粗 k/q（EPW） | 相同，例如 6³/6³ | Wannier 插值能带与 DFT 一致；插值声子与 `matdyn.x` 一致 |
| 细 k/q（EPW） | 先 24³/12³，再加密 | λ、T_c 收敛；`nkf` 是 `nqf` 的整数倍 |
| `fsthick` | 0.2–0.4 eV | 应覆盖若干倍最高声子能量和 `degaussw` |
| `degaussw` | 0.02–0.05 eV | 与细 k 网格一起测试 |
| μ\* | 0.1–0.13 | 报告时注明；最好给出 T_c 随 μ\* 变化的范围 |

> [!TIP]
> 把"Wannier 插值能带 vs. DFT 能带"和"EPW 插值声子 vs. `matdyn.x` 声子"这两张对比图作为每次 EPW 计算的必做检查。这两张图对不上，后面的 λ 和 T_c 都不可信。

---

## 常见问题与报错

| 报错/现象 | 原因 | 解决 |
|---|---|---|
| `Error in routine phq_readin: Electron-phonon only for metals` | scf 用了 `occupations='fixed'`，或体系是绝缘体 | 金属用 `occupations='smearing'`；半导体的电声问题（迁移率、带隙重整化）改用 EPW 或 `electron_phonon='ahc'` |
| `El-ph needs a DeltaVscf file` | `ph.x` 输入没有 `fildvscf` | 加上 `fildvscf='xxdv'` |
| `el-ph with images not implemented` | `ph.x` 电声计算用了 `-ni` | 去掉 `-ni`，只用 `-nk`；或手动用 `start_q`/`last_q` 拆分 |
| `Electron-phonon with Hubbard U is not supported` | DFT+U 体系的 `ph.x` 电声 | ph.x 路线不支持；EPW 6.1 起支持 DFPT+U，参见 EPW 文档 |
| `lgamma_gamma and electron-phonon are incompatible` | 只算 Γ 点且启用了 Γ 技巧 | 电声计算用一般 q 网格；必要时设 `nogg=.true.` |
| `ph.x` 读不到 `a2Fsave`（打开文件出错） | 没有先跑 `la2F=.true.` 的稠密 SCF，或 `prefix`/`outdir` 不一致 | 按 4.1 节的顺序重跑 |
| λ 随 σ 单调变化、没有平台 | `a2Fsave` 的 k 网格太粗或 q 网格太粗 | 加密 `al.scf.fit.in` 的 k 网格；加密 q 网格 |
| `lambda.x`：`wrong/too many gauss.broad.` 或 `inconsistent DOS(Ef) read` | `elph.inp_lambda.*` 来自不同的计算，或 σ 设置不一致 | 删除旧的 `elph_dir/` 后重跑 ph.x |
| `lambda.x` 给出 NaN | 有虚频（负 ω²） | 先解决虚频：结构是否充分弛豫、`tr2_ph`、截断能、smearing |
| T_c 为 0 或极小 | λ 很小，或 μ\* 太大 | 检查 λ；Allen-Dynes 公式在 λ ≲ μ\*(1+0.62λ) 时没有意义 |
| EPW：Wannier spread 很大或不收敛 | 投影、冻结窗口、`bands_skipped` 设置不当 | 先用独立的 Wannier90 流程（[第15章](15-wannier90.md)）调好参数，再搬到 `wdata`/`proj` 里 |
| EPW：插值能带在费米面附近与 DFT 不一致 | 冻结窗口没有覆盖费米面附近的能带 | 调整 `dis_froz_min/max`，确认其包含 E_F ± 若干 eV |
| EPW：读取 `save/` 失败 | 没有运行 `pp.py`，或 `dvscf_dir` 路径错 | 重新运行 `pp.py`，检查 `dvscf_dir` |
| EPW：`epwread` 重启报错 | 没设 `kmaps=.true.`，或缺少第一步写出的 `*.fmt` 文件 | 设置 `kmaps=.true.`，确认 `epwdata.fmt`、`crystal.fmt`、`*.epmatwp` 在运行目录/`outdir` 中 |
| `pp.py` 报找不到 `_ph0/pb.phsave/control_ph.xml` 或 `pb.save/data-file-schema.xml` | `ph.x`/`pw.x` 的 `outdir` 不是 `pp.py` 的运行目录 | 在 `outdir` 所在目录运行 `pp.py`，或像本例一样用 `outdir='./'` |
| `Error in routine elphsum: wrong file read` | `a2Fsave` 中的能带数与声子用的 SCF 不同（两次 SCF 的 `nbnd` 不一致） | 两次 SCF 用相同的 `nbnd`（或都用默认值） |
| EPW 跑到一半内存不足 | 细网格太密、`etf_mem=0` | 增加 `-nk`，设置 `epw_memdist=.true.`，或改用 `etf_mem=1`；只做各向同性时用 `a2f_iso` |
| 旧 EPW 教程的输入报 unknown variable | 关键字在新版本中已改名或删除 | 对照 [EPW 输入文档](https://docs.epw-code.org/Inputs/Inputs.html) 和 [Releases](https://docs.epw-code.org/Releases.html) |

---

## 练习

1. 在 Al 例子中，分别用 `a2Fsave` k 网格 16³、24³、32³ 重跑第 1–5 步，画出 λ 对 σ 的曲线，观察平台怎样随 k 网格出现。
2. 把 q 网格从 4×4×4 改为 6×6×6（注意 k 网格要与之相容），比较 λ、ω_log 和 T_c 的变化。
3. 在 `lambda.x` 中把 μ\* 从 0.08 改到 0.15，画出 T_c(μ\*)，体会 T_c 对 μ\* 的敏感程度。
4. 用 `matdyn.x`（`la2F=.true.`、`dos=.false.`）画出 Al 沿 Γ–X–W–L–Γ 的声子色散，并用 `gam.lines` 中的线宽给色散曲线加误差带，找出耦合最强的支和区域。
5. 在 Pb 的 EPW 例子中，比较细网格 24³/12³ 与 36³/18³ 得到的 λ 和 `Estimated Tc using McMillan` 值。
6. 选做：参照 EPW Tutorial 4 的 MgB₂ 练习，体会各向异性 Eliashberg 方程为什么会给出两个能隙。

---

## 参考资料

**官方文档**

- [INPUT_PH（`electron_phonon`、`el_ph_sigma`、`el_ph_nsigma`、`fildvscf`、`trans`）](https://www.quantum-espresso.org/Doc/INPUT_PH.html)
- [INPUT_MATDYN（`la2F`、`dos`）](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html) / [INPUT_Q2R](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html)
- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html)（`la2F` 未在页面列出，见源码 `Modules/read_namelists.f90`）
- QE 官方电声示例：[PHonon/examples/example03](https://github.com/QEF/q-e/tree/develop/PHonon/examples/example03)（fcc Al，`run_example` 中给出完整流程）
- `lambda.x` 的输入格式：源码 [PHonon/PH/lambda.f90](https://github.com/QEF/q-e/blob/develop/PHonon/PH/lambda.f90) 开头的注释
- [EPW 文档首页](https://docs.epw-code.org/) / [输入关键字（v6.1）](https://docs.epw-code.org/Inputs/Inputs.html) / [下载与安装](https://docs.epw-code.org/DownloadAndInstall.html) / [教程列表](https://docs.epw-code.org/Tutorials.html) / [EPW Schools](https://docs.epw-code.org/Schools.html) / [用户论坛](https://forum.epw-code.org/)

**课程与教程（只链接）**

- [EPW 官方教程 Tutorial 4：Superconducting properties](https://docs.epw-code.org/tutorials/tutorial_04/index.html)：Pb 各向同性和 MgB₂ 各向异性超导的完整上机流程，本章 EPW 部分的参数来源
- [EPW Tutorial 1：电声矩阵元插值](https://docs.epw-code.org/tutorials/tutorial_01/index.html)：入门必看，讲 Wannier 化质量检查
- [pranabdas/espresso](https://pranabdas.github.io/espresso/)：从安装到 EPW 的英文教程网站，有电声/超导部分的流程
- [nguyen-group/QE-SSP](https://github.com/nguyen-group/QE-SSP)：《Quantum ESPRESSO Course for Solid-State Physics》配套算例，覆盖从收敛测试到电声耦合
- [Imperial MSE404-MM](https://imperialcollegelondon.github.io/MSE404-MM/)：Lab 8（Vibrational Normal Modes and Phonons）是本章前置知识的好教材
- Cornell PARADIM 2016/2018 暑校（[2018](https://www.paradim.org/cu_ss_2018)）：Giustino 讲授的声子与电声原理讲座和视频，适合补理论；上机部分基于旧版 QE

**相关仓库**

- [QEF/q-e](https://github.com/QEF/q-e)：`PHonon/examples/`、`EPW/examples/`、`test-suite/epw_*` 是最权威的输入范例
- [janberges/elphmod](https://github.com/janberges/elphmod)：读取 QE/EPW 数据构建电声模型的 Python 包
- [phonopy/phonopy](https://github.com/phonopy/phonopy)：有限位移声子，可与 ph.x 的声子部分交叉检查

**本仓库脚本**：`qe_phonon_plot.py`（画声子色散/DOS，标记虚频）、`qe_parse_out.py`（解析 pw.x 输出）、`job_qe_slurm.sh`（作业模板）。用 `python ../scripts/<name> -h` 查看用法。

**关键文献**

- W. L. McMillan, *Phys. Rev.* **167**, 331 (1968). DOI: [10.1103/PhysRev.167.331](https://doi.org/10.1103/PhysRev.167.331)
- P. B. Allen, R. C. Dynes, *Phys. Rev. B* **12**, 905 (1975)，Allen-Dynes 公式. DOI: [10.1103/PhysRevB.12.905](https://doi.org/10.1103/PhysRevB.12.905)
- S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, *Rev. Mod. Phys.* **73**, 515 (2001)，DFPT 综述. DOI: [10.1103/RevModPhys.73.515](https://doi.org/10.1103/RevModPhys.73.515)
- F. Giustino, *Rev. Mod. Phys.* **89**, 015003 (2017)，第一性原理电声相互作用综述. DOI: [10.1103/RevModPhys.89.015003](https://doi.org/10.1103/RevModPhys.89.015003)
- F. Giustino, M. L. Cohen, S. G. Louie, *Phys. Rev. B* **76**, 165108 (2007)，Wannier 插值电声矩阵元. DOI: [10.1103/PhysRevB.76.165108](https://doi.org/10.1103/PhysRevB.76.165108)
- E. R. Margine, F. Giustino, *Phys. Rev. B* **87**, 024505 (2013)，基于 Wannier 函数的各向异性 Migdal-Eliashberg 理论. DOI: [10.1103/PhysRevB.87.024505](https://doi.org/10.1103/PhysRevB.87.024505)
- S. Poncé, E. R. Margine, C. Verdi, F. Giustino, *Comput. Phys. Commun.* **209**, 116 (2016)，EPW. DOI: [10.1016/j.cpc.2016.07.028](https://doi.org/10.1016/j.cpc.2016.07.028)
- H. Lee et al., *npj Comput. Mater.* **9**, 156 (2023)，EPW 代码综述. DOI: [10.1038/s41524-023-01107-3](https://doi.org/10.1038/s41524-023-01107-3)
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017)，QE 进展. DOI: [10.1088/1361-648X/aa8f79](https://doi.org/10.1088/1361-648X/aa8f79)

---

[← 上一章](15-wannier90.md) | [返回目录](../README.md) | [下一章 →](17-md-cp.md)
