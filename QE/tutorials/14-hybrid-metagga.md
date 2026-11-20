# 第14章 杂化泛函与 meta-GGA（Hybrid Functionals & meta-GGA: HSE, PBE0, ACE, SCAN）

PBE 等半局域泛函系统性低估带隙。本章讲 Quantum ESPRESSO 中两类常用的"超越 GGA"方法：杂化泛函（`input_dft = 'hse'`、`'pbe0'`、`'b3lyp'`，EXX q 网格 `nqx1/2/3`、`ecutfock`、`exxdiv_treatment`、`x_gamma_extrapolation`、默认开启的 ACE 自适应压缩交换算符）和 meta-GGA（`'scan'`、`'r2scan'`，需要 Libxc 与模守恒赝势）。重点讲清楚 pw.x 外层 EXX 循环的输出怎么读、为什么杂化泛函不能直接用 `calculation='bands'` 画能带（推荐 Wannier 插值），以及 US/PAW 赝势、应力、ph.x 等方面的限制。示例体系为 Si（PBE vs HSE 带隙）和 MgO（PBE0）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 HSE/PBE0 的 pw.x 输入，合理设置 `nqx1/2/3`、`ecutfock`、`nbnd`，读懂 `est. exchange err (dexx)` 和 `EXX self-consistency reached`；
- 知道 pw.x 杂化泛函的限制：US/PAW 下不能算力（relax/md 报错）、`ecutfock` 只对 NC 有效、应力在 7.6.1 之前有误；
- 用"显式全网格自洽 + Wannier90"得到杂化泛函能带，知道零权重 k 点法的前提；
- 在启用 Libxc 的 pw.x 中运行 SCAN/r2SCAN，知道 meta-GGA 对赝势和非共线的限制；
- 用 `-nk`、`-nb` 等并行选项控制杂化泛函的计算量。

| 目标 | 推荐 | 说明 |
|---|---|---|
| 半导体/绝缘体带隙、能带排列、缺陷能级 | HSE（`input_dft='hse'`） | 最常用的"可信带隙"方法，代价高 |
| 分子、大带隙绝缘体 | PBE0 / B3LYP | 无屏蔽，k 点收敛慢，分子用 Γ 点盒子较常见 |
| 结构与能量兼顾、代价适中 | r2SCAN（Libxc） | 比 PBE 均衡；只能用 NC 赝势 |
| 金属、磁性体系常规计算 | PBE / PBE+U | 杂化泛函对金属代价高、未必更好，见 [第8章](08-dftu-hp.md) |
| 声子、介电、拉曼 | 不能用杂化/meta-GGA 的 ph.x | 见 [第13章](13-dielectric-ir-raman.md) 的功能矩阵 |
| 准粒子能带、激子 | GW / BSE（Yambo） | 超出本章 |

---

## 2. 原理速览

### 2.1 Jacob 阶梯与 QE 中的名字

LDA → GGA（依赖 ∇ρ）→ meta-GGA（再依赖动能密度 τ）→ 杂化泛函（混入非局域 Hartree–Fock 交换）。QE 用 `input_dft` 覆盖赝势文件中的泛函；QE 7.6 源码 `XClib/qe_dft_list.f90` 中的部分短名称如下：

| `input_dft` | 类型 | 说明 |
|---|---|---|
| `'pbe0'` | 全局杂化 | 25% 精确交换（`exx_fraction` 默认 0.25） |
| `'hse'` | 屏蔽杂化 | 短程 25% 精确交换，`screening_parameter` 默认 0.106 bohr⁻¹ |
| `'b3lyp'` | 全局杂化 | `exx_fraction` 默认 0.20 |
| `'gaupbe'` | 屏蔽杂化（Gau-PBE） | 无 q → 0 发散，需 `exxdiv_treatment='none'`、`x_gamma_extrapolation=.false.` |
| `'scan'`、`'r2scan'`、`'rscan'`、`'scan0'`、`'tb09'` | meta-GGA（SCAN0 为杂化） | 短名称映射到 **Libxc**，没有 Libxc 时报 `libxc needed for this functional` |
| `'tpss'`、`'m06l'` | meta-GGA | QE 内置实现 |

> [!NOTE]
> 关于 HSE06 的参数：QE 的 `screening_parameter = 0.106` bohr⁻¹ ≈ 0.20 Å⁻¹，与 VASP 中 HSE06 的 `HFSCREEN = 0.2`（Å⁻¹）相当；QE 源码注释指出 Gaussian 中 HSE06 取 0.11 bohr⁻¹，并且梯度标度关系略有不同，两者差异很小（单个 He 原子约 1e-5 Ha）。

### 2.2 精确交换的代价

Fock 算符作用在一个 Bloch 态上需要对所有 q（= k − k'）和所有占据带求和，每一对都要做 FFT 求解泊松方程。QE 官方 `PW/examples/EXX_example/README` 给出的估计是代价约为 `NBND × NQS × (2 FFT + ...)`，即使 NQS = 1，也比普通 GGA 计算至少贵一个数量级。

两种降低代价的手段：

- **q 网格**（`nqx1/2/3`）：Fock 算符中 q 的采样可以比 k 网格粗，默认等于 k 网格。HSE 的短程屏蔽让它对 q 网格收敛更快。
- **ACE**（Adaptively Compressed Exchange，Lin Lin 2016）：把 Fock 算符压缩成低秩形式，内层 SCF 中反复使用，只在外层循环更新。`ace` 默认 `.true.`；设为 `.false.` 时用标准算法，官方文档说明"much slower"。

### 2.3 双层自洽循环

pw.x 不直接对密度矩阵自洽，而是：

1. 先做普通（非杂化）自洽，得到初始轨道；
2. 用当前轨道构造 Fock 算符 V_x[ψ]，固定它做内层 SCF；
3. 计算 dexx（一个正定的误差估计）；若大于阈值，用新轨道重建 V_x，回到第 2 步。

所以输出中会看到多段 SCF，每段之后打印 `est. exchange err (dexx)`。外层循环次数由 `exx_maxstep`（默认 100）控制。

### 2.4 为什么杂化泛函不能直接用 `calculation='bands'`

GGA 的哈密顿量只依赖电荷密度；杂化泛函的 Fock 算符依赖**所有 k+q 处的占据轨道**。对能带路径上的任意 k 点，pw.x 需要 k+q（q 在 EXX 网格上）处的轨道，而这些点不在自洽网格中。QE 官方 EXX_example README 写明：在不同于自洽计算的 k 点上做非自洽计算"at present it is not possible"。实际可行的做法见 §6：Wannier 插值（推荐）或零权重 k 点。

### 2.5 meta-GGA

SCAN（2015）满足全部 17 个已知精确约束，但数值上对实空间网格敏感；r2SCAN（2020）做了正则化，数值稳定性更好。pw.x 中 meta-GGA 需要动能密度 τ，因此：

- **只支持模守恒赝势**：US/PAW 报错 `Meta-GGA not implemented with USPP/PAW`；
- **不支持非共线**：`Non-collinear Meta-GGA not implemented`；
- 赝势含 NLCC 时打印 `nonlinear core correction is not consistent with meta-GGA` 警告；QE 7.6 新增"NC 赝势 + meta-GGA 的 NLCC 支持"和"ONCV 赝势用于 r2SCAN（模型原子动能密度）"。

---

## 3. 通用准备

示例目录：[../examples/14-hybrid-metagga/README.md](../examples/14-hybrid-metagga/README.md)。

### 3.1 赝势：为什么用模守恒

| 功能 | NC | US | PAW | 依据 |
|---|---|---|---|---|
| 杂化泛函能量 | 可以 | 可以 | 可以 | 7.6 源码 `exx.f90` 含 US 与 PAW（`PAW_init_fock_kernel`）分支 |
| 杂化泛函 + relax/md（力） | 可以 | 报错 | 报错 | `setup.f90`：`forces for hybrid functionals + US/PAW not implemented` |
| `ecutfock` 加速 | 可以 | 被忽略 | 被忽略 | INPUT_PW：Not implemented for US-PP and PAW |
| 非共线杂化泛函 | 可以 | 报错 | 报错 | `Noncolinear hybrid XC for USPP not implemented` |
| 杂化泛函应力 | 7.6.1 修正 | 7.6.1 起不再计算 | 同 US | 开发版 release notes |
| meta-GGA | 可以 | 报错 | 报错 | `Meta-GGA not implemented with USPP/PAW` |

因此本章统一用 **SG15 ONCV**（PBE，模守恒，无 NLCC）：

```bash
mkdir -p pseudo && cd pseudo
wget http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2020-02-06.tar.gz
tar -xzf sg15_oncv_upf_2020-02-06.tar.gz Si_ONCV_PBE-1.2.upf Mg_ONCV_PBE-1.2.upf O_ONCV_PBE-1.2.upf
grep -m1 'core_correction=' *.upf          # "F" for all three
grep -m1 'z_valence=' Mg_ONCV_PBE-1.2.upf   # Mg includes 2s2p semicore: 10 electrons
cd ..
```

也可以用 PseudoDojo 的 NC 赝势（`nc-sr-04_pbe_standard_upf.tgz` 中的 `Si.upf` 等），它们大多含 NLCC，pw.x 会打印警告但可以运行。赝势的一般选择原则见 [第3章](03-pseudopotentials.md)。

> [!IMPORTANT]
> 目前没有为杂化泛函专门生成的赝势：通常用 PBE 赝势 + `input_dft` 覆盖（官方 EXX README 也是这样说明的）。这本身是一种近似，比较不同代码的杂化泛函带隙时要说明所用赝势。

### 3.2 先做 PBE 参考计算

```fortran
! Si (diamond), PBE reference run with SG15 ONCV norm-conserving pseudopotentials.
! Same cell, cutoff and k mesh as si.hse.scf.in, so that gaps can be compared directly.
&CONTROL
  calculation = 'scf'
  prefix      = 'si_pbe'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 40.0      ! start value for SG15 Si; converge it
  ecutrho   = 160.0     ! 4 x ecutwfc (norm-conserving)
  nbnd      = 8         ! 4 occupied + 4 empty bands -> HOMO/LUMO printed
/
&ELECTRONS
  conv_thr  = 1.0d-10
/
ATOMIC_SPECIES
  Si  28.0855  Si_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
  -2.7150   0.0000   2.7150
   0.0000   2.7150   2.7150
  -2.7150   2.7150   0.0000
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  4 4 4 0 0 0
```

```bash
mpirun -np 8 pw.x -nk 4 -in si.pbe.scf.in > si.pbe.scf.out
grep "highest occupied" si.pbe.scf.out
```

`nbnd` 必须大于占据带数（Si 原胞 8 个价电子 → 4 条占据带），否则只打印 `highest occupied level`，读不到导带。

---

## 4. HSE 自洽计算

### 4.1 输入

```fortran
! Si (diamond), HSE (QE 'hse': 25% short-range exact exchange, mu = 0.106 bohr^-1).
! pw.x first converges a plain PBE SCF, then switches on exact exchange (outer EXX loop).
&CONTROL
  calculation = 'scf'
  prefix      = 'si_hse'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 40.0
  ecutrho   = 160.0
  nbnd      = 8
  input_dft = 'hse'      ! overrides the PBE functional read from the UPF files
  nqx1 = 2, nqx2 = 2, nqx3 = 2   ! q mesh for the Fock operator (default = k mesh); must divide the k mesh
  ! ecutfock = 80.0      ! optional speed-up (NC only), must be >= ecutwfc and <= ecutrho
  ! exxdiv_treatment = 'gygi-baldereschi' (default), x_gamma_extrapolation = .true. (default)
  ! ace = .true. (default): Adaptively Compressed Exchange
/
&ELECTRONS
  conv_thr  = 1.0d-10
/
ATOMIC_SPECIES
  Si  28.0855  Si_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
  -2.7150   0.0000   2.7150
   0.0000   2.7150   2.7150
  -2.7150   2.7150   0.0000
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  4 4 4 0 0 0
```

### 4.2 关键词表

以下默认值据 INPUT_PW（QE 7.5/7.6 文档）与 7.6 源码核对：

| 关键词（namelist） | 默认值 | 推荐 / 说明 |
|---|---|---|
| `input_dft`（SYSTEM） | 从赝势读取 | `'hse'`、`'pbe0'`、`'b3lyp'` 等；官方提示"Use with care" |
| `nqx1, nqx2, nqx3`（SYSTEM） | 等于 k 网格 | EXX 的 q 网格，可以小于 k 网格；必须与 k 网格相容（k+q 在 k 列表中） |
| `exx_fraction`（SYSTEM） | 随泛函 | PBE0 为 0.25，B3LYP 为 0.20；做"调 α"的杂化泛函时修改 |
| `screening_parameter`（SYSTEM） | 0.106 | HSE 类泛函的屏蔽参数（bohr⁻¹） |
| `ecutfock`（SYSTEM） | 文档写 `ecutrho`；7.6 源码未设置时实际取 min(`ecutrho`, 4·`ecutwfc`) | Fock 算符的截断（Ry）；减小可显著加速；必须在 [`ecutwfc`, `ecutrho`] 之间；US/PAW 下被忽略；INPUT_PW 说明应力计算未实现该近似；金属中慎用 |
| `exxdiv_treatment`（SYSTEM） | `'gygi-baldereschi'` | q → 0 发散的处理；强各向异性晶胞可用 `'vcut_ws'`（配合 `ecutvcut`），`'vcut_spherical'` 适合立方/准立方 |
| `x_gamma_extrapolation`（SYSTEM） | `.true.` | G = 0 项外推；与 `'vcut_ws'`、`'vcut_spherical'` 不能同时使用（源码报 `cannot USE x_gamm_extrap and vcut_ws`），GAU-PBE 需设为 `.false.` |
| `ecutvcut`（SYSTEM） | 0.0 Ry | `vcut_ws` 的倒空间截断，官方 README 建议在 0.7–2.0 范围内测试收敛 |
| `ace`（SYSTEM） | `.true.` | ACE 算符；保持默认 |
| `localization_thr`（SYSTEM） | 0.0 | 大于 0 时启用 SCDM 局域化并只计算重叠超过阈值的轨道对；US/PAW 下不支持 |
| `exx_maxstep`（ELECTRONS） | 100 | 外层 EXX 循环最大次数 |
| `adaptive_thr`（ELECTRONS） | `.false.` | 内层 SCF 使用自适应阈值；配合 `conv_thr_init`（默认 1e-3）和 `conv_thr_multi`（默认 0.1） |
| `nbnd`（SYSTEM） | 绝缘体为占据带数 | 想看带隙必须加空带 |

### 4.3 运行与输出

```bash
mpirun -np 16 pw.x -nk 4 -in si.hse.scf.in > si.hse.scf.out
grep -e "EXX: q-point mesh" -e "EXX: setup a grid" si.hse.scf.out
grep -e "est. exchange err" -e "EXX self-consistency reached" si.hse.scf.out
grep -e "^!!" si.hse.scf.out
grep "highest occupied" si.hse.scf.out | tail -1
```

输出片段（示意输出，数值仅供参考）：

```text
     EXX: q-point mesh:     2    2    2
     EXX: setup a grid of NN q-points centered on each k-point
     ...
     EXX: now go back to refine exchange calculation

     total energy              =     -xx.xxxxxxxx Ry
     Harris-Foulkes estimate   =     -xx.xxxxxxxx Ry
     est. exchange err (dexx)  =       0.000xxxxx Ry
     - averaged Fock potential =       x.xxxxxxxx Ry
     + Fock energy (ACE)       =      -x.xxxxxxxx Ry
     ...
!!   total energy              =     -xx.xxxxxxxx Ry
     Harris-Foulkes estimate   =     -xx.xxxxxxxx Ry
     est. exchange err (dexx)  =          x.xE-xx Ry
     ...
     EXX self-consistency reached
```

读法与成功判据：

- `EXX: q-point mesh` 应回显你设置的 `nqx`。如果出现 `EXX: WARNING: q-point mesh has been updated!`，说明 q 网格与 k 网格不相容，pw.x 自动换成了别的 q 网格（源码会依次尝试更小或更大的网格，最坏退到 1×1×1）。
- 第一段 SCF 是普通 PBE；之后每段结束打印一次 `est. exchange err (dexx)`，应单调下降。
- 最终总能行以 `!!` 开头（而不是普通计算的 `!`），随后打印 `EXX self-consistency reached`。脚本里提取总能时要匹配 `^!!`。
- 如果到达 `exx_maxstep` 仍未收敛，打印 `EXX convergence NOT achieved after N iterations: stopping`。
- 最后一次 `highest occupied, lowest unoccupied level` 才是杂化泛函的结果；前面的是 PBE 阶段的。
- 这里的带隙是**网格点上**的最小间隙。Si 的导带底在 Γ–X 线上约 85% 处，不在 4×4×4 网格点上，真正的间接带隙要看 §6 的插值能带。

> [!TIP]
> 想更快得到 HSE 能量：先单独跑一次 PBE，再用 `startingwfc = 'file'`、`startingpot = 'file'` 并让 HSE 计算使用同一 `prefix` 与 `outdir` 读取它。注意 HSE 运行会覆盖 save 目录，必要时先复制一份。

### 4.4 q 网格与截断的收敛

| 测试 | 做法 | 判据 |
|---|---|---|
| `nqx` | 1 → 2 → 4（k 网格 4×4×4，q 必须整除 k） | 带隙变化 < 0.02 eV |
| k 网格 | 4³ → 6³ → 8³（同时调整 `nqx`） | 网格点上的能量与带隙稳定 |
| `ecutfock` | `ecutrho` → 0.75·`ecutrho` → 0.5·`ecutrho` | 总能差 < 1 meV/atom 才可接受 |
| `ecutwfc` | 与 PBE 相同方法 | 杂化泛函一般不需要比 PBE 更高的截断 |

```bash
for nq in 1 2 4; do
  sed "s/nqx1 = 2, nqx2 = 2, nqx3 = 2/nqx1 = $nq, nqx2 = $nq, nqx3 = $nq/; s/si_hse/si_hse_nq$nq/" \
      si.hse.scf.in > si.hse.nq$nq.in
  mpirun -np 16 pw.x -nk 4 -in si.hse.nq$nq.in > si.hse.nq$nq.out
  echo "nq=$nq  $(grep 'highest occupied' si.hse.nq$nq.out | tail -1)"
done
```

仓库脚本 `qe_conv_test.sh` 只内置 `ecut`/`kpt`/`smear` 三种扫描模式，但支持 `-s KEY=VALUE` 填充任意 `@KEY@` 占位符；在示例目录下用 `bash ../../scripts/qe_conv_test.sh -h` 查看用法，用 `python ../../scripts/qe_parse_out.py -h` 查看如何汇总总能与带隙。

---

## 5. PBE0：MgO

```fortran
! MgO (rocksalt), PBE0 with SG15 ONCV norm-conserving pseudopotentials.
&CONTROL
  calculation = 'scf'
  prefix      = 'mgo_pbe0'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tprnfor     = .true.
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 2
  ecutwfc   = 80.0      ! start value (O and Mg semicore are hard); converge it
  ecutrho   = 320.0
  nbnd      = 12        ! 8 occupied (16 valence electrons: Mg 10 + O 6) + 4 empty
  input_dft = 'pbe0'    ! 25% exact exchange (exx_fraction default for PBE0)
  nqx1 = 2, nqx2 = 2, nqx3 = 2
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.5
/
ATOMIC_SPECIES
  Mg  24.305   Mg_ONCV_PBE-1.2.upf
  O   15.999   O_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
   0.0000   2.1060   2.1060
   2.1060   0.0000   2.1060
   2.1060   2.1060   0.0000
ATOMIC_POSITIONS crystal
  Mg  0.00  0.00  0.00
  O   0.50  0.50  0.50
K_POINTS automatic
  4 4 4 0 0 0
```

```bash
mpirun -np 16 pw.x -nk 4 -in mgo.pbe0.scf.in > mgo.pbe0.scf.out
grep "highest occupied" mgo.pbe0.scf.out | tail -1
grep -A4 "Forces acting on atoms" mgo.pbe0.scf.out | tail -4
```

MgO 的价带顶和导带底都在 Γ 点（Γ 在 4×4×4 网格上），所以网格点给出的就是直接带隙。PBE0 没有屏蔽，长程交换对 q 网格收敛更慢，`nqx` 要单独测试。

杂化泛函 + NC 赝势可以做 `relax`（固定晶胞）。官方 EXX README 写明："Energy and forces ... So structural optimization is OK if the cell shape is kept fixed"。对 vc-relax 见 §8 的应力说明。

---

## 6. 杂化泛函能带

### 6.1 推荐路线：显式全网格自洽 + Wannier90

思路：在规则网格上做完 HSE 自洽，用 Wannier90 构造最大局域 Wannier 函数，再在任意 k 路径上插值。Wannier 插值不需要再次计算 Fock 算符。详细原理与参数见 [第15章 Wannier90 接口与插值能带](15-wannier90.md)，这里只列出与杂化泛函有关的部分。

与 PBE 流程（scf → nscf → pw2wannier90）的区别：杂化泛函不能在新 k 点上做 nscf，所以**直接在自洽这一步**就把全网格显式写出来（`K_POINTS crystal`，与 Wannier90 的 `kpoints` 块顺序一致），并加 `nosym`/`noinv`：

```fortran
! (excerpt of si.hse.fullgrid.in; the full 64-point list is in the example directory)
&SYSTEM
  ...
  nbnd      = 12          ! 4 occupied + 8 empty: enough for 8 Wannier functions (sp3 on both Si)
  input_dft = 'hse'
  nqx1 = 2, nqx2 = 2, nqx3 = 2
  nosym     = .true.      ! keep the k list exactly as given (same order as the .nnkp file)
  noinv     = .true.
/
...
K_POINTS crystal
64
  0.00000000  0.00000000  0.00000000  1.562500e-02
  0.00000000  0.00000000  0.25000000  1.562500e-02
  ...
```

> [!IMPORTANT]
> 显式 k 列表时，pw.x 内部没有"自动网格"信息（源码中 `nk1 = nk2 = nk3 = 0`），EXX q 网格的默认值随之失效：**必须**显式写 `nqx1/2/3`，否则 `exx_grid_init` 报 `wrong EXX q grid`。k 列表要用 Wannier90 的 `utility/kmesh.pl 4 4 4` 生成（它输出的就是 `K_POINTS crystal` 格式与等权重），或用示例目录中已生成好的文件。

然后按第 15 章的顺序：

```bash
mpirun -np 16 pw.x -nk 4 -in si.hse.fullgrid.in > si.hse.fullgrid.out
wannier90.x -pp si_hse
mpirun -np 16 pw2wannier90.x -in si_hse.pw2wan.in > si_hse.pw2wan.out
wannier90.x si_hse
grep -A6 "Final State" si_hse.wout
```

检查：

- 冻结窗口（`dis_froz_max`）内的插值能带应与 `si.hse.fullgrid.out`（`verbosity='high'` 时打印全部 k 点本征值）在网格点上一致；
- 用插值能带读取间接带隙（导带底在 Γ–X 线上），与 §4 网格点上的"带隙"比较，体会两者的差别；
- 同一流程用 PBE 再做一遍，可以得到同一 k 路径上的 PBE vs HSE 对比图。

### 6.2 零权重 k 点法（进阶）

另一种思路是把路径上的 k 点以**零权重**加入自洽 k 列表。QE 7.6 源码 `exx.f90` 中，权重为零的 k 点的占据数被设为 0，因此它们不参与构造 Fock 算符；但 pw.x 仍需要每个路径点 k_p 的所有 k_p + q（q 在 EXX 网格上）都出现在 k 列表中，否则 `exx_grid_init` 会自动缩小 q 网格（打印 `EXX: WARNING: q-point mesh has been updated!`），结果将不正确。

因此零权重法的 k 列表必须包含：规则网格（有权重）+ 每个路径点 k_p + 所有 k_p + q（均为零权重），并显式设置 `nqx1/2/3` 和 `nosym`。点数增长很快（路径点数 × q 点数），代价往往高于 Wannier 插值。使用时务必检查输出中的 `EXX: q-point mesh` 与你设定的一致。本章不提供该方法的完整输入，推荐使用 §6.1。

### 6.3 不要这样做

| 做法 | 问题 |
|---|---|
| HSE scf 后直接 `calculation='bands'` | Fock 算符需要路径点 k+q 处的轨道，不在 save 目录中；官方 README 说明暂不支持任意 k 点非自洽 |
| PBE 能带 + 剪刀算符冒充 HSE 能带 | 只适合粗略比较，色散与带宽并不相同 |
| 只报网格点上的"HSE 带隙" | 对间接带隙材料（Si）会高估 |

---

## 7. meta-GGA：SCAN 与 r2SCAN

### 7.1 编译时启用 Libxc

QE 用户手册（Libxc 一节）给出的方法：

```bash
# configure
./configure --with-libxc --with-libxc-prefix=/path/to/libxc
# CMake
cmake -B build -DQE_ENABLE_LIBXC=ON ...   # plus your usual options
```

- QE 与 Libxc 4.3.0 及以上兼容；
- 用户手册建议：Libxc 5.0.0 以后的版本默认强制 Fermi hole curvature，可能导致某些 meta-GGA 结果不准或不收敛，编译 Libxc 时加 `--disable-fhc`；
- 检查已安装的 pw.x 是否链接了 Libxc：`ldd $(which pw.x) | grep -i xc`（动态链接时）；最直接的办法是跑一个 r2SCAN 小算例，没有 Libxc 时会报 `libxc needed for this functional`。

conda-forge 等预编译包是否带 Libxc 取决于打包方式，请自行确认（见 [第1章](01-install.md)）。

### 7.2 r2SCAN 输入

```fortran
! Si (diamond), r2SCAN meta-GGA. REQUIRES pw.x compiled with Libxc
! (otherwise: "libxc needed for this functional"). Norm-conserving pseudopotentials only.
&CONTROL
  calculation = 'scf'
  prefix      = 'si_r2scan'
  outdir      = './tmp'
  pseudo_dir  = './pseudo'
  tstress     = .true.   ! meta-GGA stress was fixed in QE 7.6; check it before vc-relax
  tprnfor     = .true.
/
&SYSTEM
  ibrav     = 0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 50.0       ! meta-GGA usually needs a higher cutoff than PBE; converge it
  ecutrho   = 200.0
  nbnd      = 8
  input_dft = 'r2scan'
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.4
/
ATOMIC_SPECIES
  Si  28.0855  Si_ONCV_PBE-1.2.upf
CELL_PARAMETERS angstrom
  -2.7150   0.0000   2.7150
   0.0000   2.7150   2.7150
  -2.7150   2.7150   0.0000
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  6 6 6 0 0 0
```

`input_dft = 'r2scan'` 在源码中被映射为 Libxc 的 497（交换）+ 498（关联）。也可以用 QE 7.0 起要求的 Libxc 索引写法 `input_dft = 'XC-000i-000i-000i-000i-497L-498L'`（六个槽依次为 LDA 交换、LDA 关联、GGA 交换、GGA 关联、MGGA 交换、MGGA 关联，`i`/`L` 分别表示 QE 内置/Libxc，写法见 `XClib` 源码注释与用户手册），两者等价。

### 7.3 meta-GGA 的注意事项

| 问题 | 说明 |
|---|---|
| 赝势 | 只能用 NC。用 PBE 赝势跑 meta-GGA 是常见折中；7.6 起 ONCV 赝势可用于 r2SCAN（模型原子动能密度） |
| 截断与网格 | SCAN 对 FFT 网格敏感，截断要单独收敛；r2SCAN 通常更稳定 |
| 应力 | QE 7.6 release notes："Stress for meta-GGA functionals was broken"，已在 7.6 修复；开发版 7.6.1 又修正了 `ecutrho > 4·ecutwfc` 时 meta-GGA 应力错误与并行相关的问题。做 vc-relax 前先用有限差分（E–V 曲线）核对应力 |
| Libxc 数值阈值 | 7.6.1（开发版）修复了 Libxc meta-GGA 在梯度很小处把能量密度置零的问题；7.6 下若遇到 nspin=1 与 nspin=2 非磁结果不一致，可能与此有关 |
| 非共线 / SOC | 不支持 |
| ph.x | 不支持 meta-GGA（`The phonon code with meta-GGA functionals is not yet available`），声子用有限位移（phonopy） |
| 色散修正 | r2SCAN 常配 rVV10 或 D3/D4，具体可用组合见 [第9章](09-vdw-surface-workfunction.md) |

---

## 8. 力、应力与结构优化

| 计算 | NC + 杂化 | US/PAW + 杂化 | NC + meta-GGA |
|---|---|---|---|
| scf 能量 | 可以 | 可以 | 可以 |
| relax / md（力） | 可以 | 报错 `forces for hybrid functionals + US/PAW not implemented` | 可以 |
| vc-relax（应力） | 7.6 及以前应力有误；7.6.1 修正 | 7.6.1 起不计算应力，vc-relax 停止 | 7.6 修复，7.6.1 再修正部分情况 |

开发分支 release notes（7.6.1，截至 2026-11 未发布）："Stress for hybrid functionals: the contribution of the q=0 divergent term was missing and the normalization was wrong"；"With ultrasoft or PAW pseudopotentials the stress for hybrid functionals is not computed, variable-cell calculations stop"。另外，QE 7.6 的 `PW/src/setup.f90` 中，杂化泛函计算应力时不能用 k 点池并行：`tstress = .true.`（vc-relax 自动打开）且 `-nk` > 1 时报 `stress for hybrid functionals not available with pools`。

实用建议：

- 晶格常数用 PBE/PBEsol/r2SCAN 优化，再做 HSE 单点（记作 HSE@PBE），发表时写清楚；
- 一定要杂化泛函结构时，用 NC 赝势、固定晶胞 relax + 手动 E–V 曲线（[第5章](05-relax-eos.md)），而不是依赖 7.6 及以前的杂化应力。

---

## 9. 并行与计算量

| 选项 | 作用 | 杂化泛函中的建议 |
|---|---|---|
| `-nk`（pools） | 按 k 点分组 | 首选，效率高；注意内存随 pool 数复制 |
| `-nb`（band groups） | 按能带分组 | 用户手册："Especially useful for calculations with hybrid functionals"；`pw2wannier90.x` 不支持 `-nb` |
| `-nd`（线性代数） | 子空间对角化 | 大体系有用 |
| OpenMP（`OMP_NUM_THREADS`） | 线程并行 | 内存紧张时用更多线程、更少 MPI 进程 |
| `ecutfock` | 降低 EXX FFT 网格 | NC 才有效，需测试误差 |
| `nqx` | 降低 q 点数 | HSE 通常可以用比 k 网格粗的 q 网格 |

```bash
# 16 MPI ranks: 4 pools x 2 band groups x 2 ranks for plane waves
mpirun -np 16 pw.x -nk 4 -nb 2 -in si.hse.scf.in > si.hse.scf.out
```

作业脚本模板见仓库的 `QE/scripts/job_qe_slurm.sh`（提交前先阅读文件头部注释，修改 module、进程数和 `-nk`/`-nb`）和 [Common 第2章 作业调度](../../Common/tutorials/02-job-scheduler.md)；更系统的并行调优见 [第18章](18-automation-parallel-troubleshooting.md)。

---

## 10. 与 VASP 的对应关系

VASP 一侧的完整讲解见 [VASP 第6章 杂化泛函与 meta-GGA](../../VASP/tutorials/06-hybrid-metagga.md)。

| 概念 | VASP | Quantum ESPRESSO |
|---|---|---|
| 打开杂化泛函 | `LHFCALC = .TRUE.` + `GGA`/`HFSCREEN`/`AEXX` | `input_dft = 'hse'`/`'pbe0'`/`'b3lyp'` |
| 精确交换比例 | `AEXX` | `exx_fraction` |
| 屏蔽参数 | `HFSCREEN = 0.2`（Å⁻¹） | `screening_parameter = 0.106`（bohr⁻¹，≈ 0.2 Å⁻¹） |
| q 网格降采样 | `NKRED`、`NKREDX/Y/Z` | `nqx1/2/3` |
| Fock FFT 网格精度 | `PRECFOCK` | `ecutfock` |
| q → 0 发散处理 | `HFRCUT` 等 | `exxdiv_treatment`、`x_gamma_extrapolation`、`ecutvcut` |
| 压缩交换 | ACE（默认） | `ace = .true.`（默认） |
| 杂化能带 | 零权重 k 点、`KPOINTS_OPT` | Wannier 插值（推荐）、零权重 k 点（需加入所有 k+q） |
| r2SCAN | `METAGGA = R2SCAN` + `LASPH` | `input_dft = 'r2scan'`（需 Libxc），NC 赝势 |
| meta-GGA 赝势 | POTCAR 需含动能密度信息 | 只能用 NC；7.6 起 ONCV 可用于 r2SCAN |
| 总能行 | `F=`、`E0=` | `!!   total energy`（杂化）/ `!    total energy`（普通） |

---

## 参数选择与收敛建议

| 参数 | 建议起点 | 收敛判据 |
|---|---|---|
| 赝势 | NC（SG15、PseudoDojo） | 需要力/应力/meta-GGA 时必须 NC |
| `ecutwfc` | 与 PBE 收敛值相同；meta-GGA 加 10–20 Ry 测试 | ΔE < 1 meV/atom |
| k 网格 | 与 PBE 相同或略粗 | 带隙变化 < 0.02 eV |
| `nqx` | k 网格的 1/2（必须整除） | 与 `nqx` = k 网格比较，带隙差 < 0.02 eV |
| `ecutfock` | 默认；需要加速时降到 0.5–0.75 倍并测试 | 总能差 < 1 meV/atom |
| `nbnd` | 占据带 + 4（看带隙）；Wannier 需要更多 | 能读出 LUMO；Wannier 窗口内能带完整 |
| `conv_thr` | 1e-8 ~ 1e-10 | 出现 `EXX self-consistency reached` |
| `exx_maxstep` | 默认 100 | 一般 5–15 次外层循环收敛；不收敛先查结构与 SCF 设置 |
| `exxdiv_treatment` | 体材料默认；强各向异性晶胞（slab、长链）试 `'vcut_ws'` | 结果对盒子尺寸与 k 网格不敏感 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 处理 |
|---|---|---|
| `forces for hybrid functionals + US/PAW not implemented` | US/PAW 赝势 + relax/md | 换 NC 赝势，或只做单点 |
| `Noncolinear hybrid XC for USPP not implemented` | 非共线 + 杂化 + US | 换 NC 赝势 |
| `Warning: ecutfock not valid for US/PAW, ignored` | US/PAW 下设置了 `ecutfock` | 正常提示；需要加速时换 NC |
| `ecutfock can not be < ecutwfc or > ecutrho!` | `ecutfock` 越界 | 设在 `ecutwfc` 与 `ecutrho` 之间 |
| `cannot USE x_gamm_extrap and vcut_ws`（或 `vcut_spherical`） | 两种发散处理不能同时用 | 设 `x_gamma_extrapolation = .false.` |
| `invalid exxdiv_treatment: ...` | 拼写错误 | 只能是 `gygi-baldereschi`、`vcut_spherical`、`vcut_ws`、`none` |
| `EXX: WARNING: q-point mesh has been updated!` | `nqx` 与 k 网格不相容，或显式 k 列表缺少 k+q | 让 `nqx` 整除 k 网格；显式列表时补齐 k+q 点 |
| `EXX convergence NOT achieved after N iterations: stopping` | 外层循环不收敛 | 检查结构；收紧 `conv_thr`；增大 `exx_maxstep`；金属体系慎用杂化泛函 |
| `ACE is turned off because number of occupied orbitals <= 1` | 体系总价电子数 ≤ 1（如单个 H 原子；源码判据是 `nelec <= 1`） | 正常提示，使用标准算法 |
| `stress for hybrid functionals not available with pools` | 杂化泛函 + `tstress`（或 vc-relax）+ `-nk` > 1 | 去掉 `-nk`（改用 `-nb`/平面波并行），或不算应力 |
| `wrong EXX q grid` | 显式 k 列表（`K_POINTS crystal`/`tpiba`）却没写 `nqx1/2/3` | 显式写 `nqx1/2/3` |
| `localization without ACE not implemented` | `localization_thr > 0` 且 `ace = .false.` | 保持 `ace = .true.` |
| `hybrid XC and electric fields untested` | 杂化泛函 + `lelfield` | 不支持，改用 PBE 做有限电场 |
| BEWARE: `nonlinear core correction is not consistent with hybrid XC`（或 meta-GGA） | 赝势含 NLCC | 提示性；精确计算换无 NLCC 的赝势（如 SG15） |
| `Meta-GGA not implemented with USPP/PAW` | meta-GGA + US/PAW | 换 NC 赝势 |
| `Non-collinear Meta-GGA not implemented` | meta-GGA + `noncolin` | 不支持，改用共线计算 |
| `libxc needed for this functional` | pw.x 未链接 Libxc，却用了 `scan`/`r2scan`/`tb09` 等 | 重新编译（`--with-libxc` 或 `-DQE_ENABLE_LIBXC=ON`） |
| `The phonon code with hybrid functionals is not yet available` | 在杂化泛函基态上跑 ph.x | 用有限位移（phonopy）或回到 PBE |
| 杂化 vc-relax 的晶格常数明显异常 | 7.6 及以前杂化应力有误 | 用 E–V 曲线；等待 7.6.1 或用 NC |
| 带隙与文献差很多 | 只看了网格点；`nqx` 太小；结构不同 | 用 Wannier 插值能带；收敛 `nqx`；统一结构 |
| `pw2wannier90`：`bands (-nb) not implemented` | 在 pw2wannier90.x 上用了 `-nb` | 去掉 `-nb` |

---

## 练习

1. 对 Si 做 PBE 与 HSE 的 4×4×4 计算，记录网格点上的"带隙"；再用 §6.1 的 Wannier 流程得到插值能带，读出真正的间接带隙，比较两者。
2. 固定 k 网格 4×4×4，扫描 `nqx = 1, 2, 4`，画 HSE 带隙随 `nqx` 的变化，并记录墙钟时间。
3. 把 `ecutfock` 设为 `ecutrho` 的 0.5 倍，比较总能、带隙和耗时。
4. 对 MgO 比较 PBE、PBE0、HSE 的 Γ 点直接带隙；把 `exx_fraction` 改为 0.20、0.30，观察带隙随 α 的线性变化趋势。
5. 在启用 Libxc 的 pw.x 上用 r2SCAN 做 Si 的 E–V 曲线（[第5章](05-relax-eos.md) 的流程），与 PBE 晶格常数比较；同时检查 `tstress` 输出的压强与 E–V 曲线斜率是否一致。
6. 用 `-nk 4`、`-nk 4 -nb 2`、`-nk 2 -nb 4` 三种并行方式运行同一个 HSE 计算，比较时间与内存。

---

## 参考资料

**官方文档**

- [INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html) — `input_dft`、`ace`、`exx_fraction`、`screening_parameter`、`exxdiv_treatment`、`x_gamma_extrapolation`、`ecutvcut`、`nqx1/2/3`、`ecutfock`、`localization_thr`、`exx_maxstep`、`adaptive_thr`
- [PW/examples/EXX_example](https://github.com/QEF/q-e/tree/develop/PW/examples/EXX_example) — 官方 EXX README（算法、代价、q → 0 发散、限制）与 Si、分子算例
- [QE 用户手册](https://www.quantum-espresso.org/Doc/user_guide/) — Libxc 编译与 `input_dft` 索引写法、并行层级（`-nb`）
- [QE 7.6 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.6) / [release-notes（develop，含 7.6.1/8.0 预告）](https://gitlab.com/QEF/q-e/-/blob/develop/Doc/release-notes)
- [Libxc](https://libxc.gitlab.io/) — 泛函列表与 ID
- [Wannier90 文档](https://wannier90.readthedocs.io/) — 插值能带
- 赝势：[SG15 ONCV](http://www.quantum-simulation.org/potentials/sg15_oncv/)、[PseudoDojo](https://www.pseudo-dojo.org/)

**课程与教程（只链接）**

- [QuantumESPRESSO-school-2023](https://github.com/materialscloud-org/QuantumESPRESSO-school-2023) — EPFL/Materials Cloud 高级学校，含杂化泛函与 Wannier90 上机，最贴近本章
- [MaX e-School 2021](https://indico.ictp.it/event/9616/) — ICTP/MaX 两周课程视频，包括杂化泛函与 GPU 并行
- [ENCCS: Efficient materials modelling on HPC](https://enccs.github.io/efficient-materials-modelling-on-hpc/) — `-nk`/`-nd`/`-nt` 与 GPU 实践
- [MSE404-MM（Imperial College London）](https://imperialcollegelondon.github.io/MSE404-MM/) — QE 基础实验课，适合作为本章前置
- [PARADIM Computational Tutorials（Cornell）](https://www.paradim.org/toolbox/computation/tutorials) — QE 基础与 UV-Vis/Yambo 系列，可作为带隙问题的背景

**相关 GitHub 仓库（只链接）**

- [QEF/q-e](https://github.com/QEF/q-e) — 官方源码镜像
- [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90) — Wannier90 源码与 tutorials
- [aiidateam/aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso) — 自动化工作流

**关键文献**

- J. Heyd, G. E. Scuseria, M. Ernzerhof, *J. Chem. Phys.* 118, 8207 (2003)（HSE）. DOI: [10.1063/1.1564060](https://doi.org/10.1063/1.1564060)
- A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, G. E. Scuseria, *J. Chem. Phys.* 125, 224106 (2006)（HSE06 屏蔽参数）. DOI: [10.1063/1.2404663](https://doi.org/10.1063/1.2404663)
- C. Adamo, V. Barone, *J. Chem. Phys.* 110, 6158 (1999)（PBE0）. DOI: [10.1063/1.478522](https://doi.org/10.1063/1.478522)
- L. Lin, *J. Chem. Theory Comput.* 12, 2242 (2016)（ACE）. DOI: [10.1021/acs.jctc.6b00092](https://doi.org/10.1021/acs.jctc.6b00092)
- F. Gygi, A. Baldereschi, *Phys. Rev. B* 34, 4405 (1986)（q → 0 发散处理）. DOI: [10.1103/PhysRevB.34.4405](https://doi.org/10.1103/PhysRevB.34.4405)
- J. Sun, A. Ruzsinszky, J. P. Perdew, *Phys. Rev. Lett.* 115, 036402 (2015)（SCAN）. DOI: [10.1103/PhysRevLett.115.036402](https://doi.org/10.1103/PhysRevLett.115.036402)
- J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, J. Sun, *J. Phys. Chem. Lett.* 11, 8208 (2020)（r2SCAN）. DOI: [10.1021/acs.jpclett.0c02405](https://doi.org/10.1021/acs.jpclett.0c02405)
- I. Carnimeo et al., *J. Chem. Theory Comput.* 19, 6992 (2023)（QE 的 GPU 与性能可移植性）. DOI: [10.1021/acs.jctc.3c00249](https://doi.org/10.1021/acs.jctc.3c00249)
- P. Giannozzi et al., *J. Chem. Phys.* 152, 154105 (2020)（QE 面向百亿亿次计算）. DOI: [10.1063/5.0005082](https://doi.org/10.1063/5.0005082)

---

[← 上一章](13-dielectric-ir-raman.md) | [返回目录](../README.md) | [下一章 →](15-wannier90.md)
